#!/usr/bin/env python
# python script for underlying event studies

import argparse
import yaml
import IPython
import ROOT
import os
import DMesonJetCompare
import DMesonJetUtils
import subprocess
import numpy
import itertools

globalList = []

class Error(Exception):
    pass

class OverwriteError(Error):
    def __init__(self, mydict, field):
        self.field = field
        self.mydict = mydict

    def __str__(self):
        return "Trying to overwrite {}.\n Keys already present: {}".format(self.field, self.mydict.keys())

def Projections(hist, bins, rebinX, xmin=0, xmax=0):
    if hasattr(bins, '__iter__') and len(bins) > 1:
        res = []
        for min, max in zip(bins[:-1], bins[1:]):
            p = hist.ProjectionY("{}_{}_{}".format(hist.GetName(), int(min), int(max)), hist.GetXaxis().FindBin(min), hist.GetXaxis().FindBin(max) - 1)
            p.SetTitle("{} - {}%".format(int(min), int(max)))
            if rebinX > 1: p.Rebin(rebinX)
            if xmax > xmin: p.GetXaxis().SetRangeUser(xmin, xmax)
            p.Sumw2()
            res.append(p)
    else:
        res = hist.ProjectionY("{}_Projection".format(hist.GetName(), 0, -1))
        if rebinX > 1: res.Rebin(rebinX)
        if xmax > xmin: res.GetXaxis().SetRangeUser(xmin, xmax)
        res.Sumw2()
    return res

def CleanUpHistogramName(orig):
    return orig.replace("AliAnalysisTask", "").replace("_histos", "").replace("fHist", "").replace('/', "_").replace("__", "_")

def GenerateStdDevProfile(h, name, title, yaxisTitle):
    prof_stddev = h.ProfileX("{}_profSigma".format(name), 1, -1, "s")
    if h.GetXaxis().GetXbins().GetSize() == h.GetNbinsX() + 1:
        std_dev = ROOT.TH1F(name, title, h.GetNbinsX(), h.GetXaxis().GetXbins().GetArray())
    else:
        std_dev = ROOT.TH1F(name, title, h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    std_dev.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
    std_dev.GetYaxis().SetTitle(yaxisTitle)
    std_dev.GetYaxis().SetTitleOffset(1)
    for i in range(0, prof_stddev.GetNbinsX() + 2):
        if prof_stddev.GetBinContent(i) == 0: continue
        std_dev.SetBinContent(i, prof_stddev.GetBinError(i) / prof_stddev.GetBinContent(i))
    return std_dev

class RhoDefinition:
    def __init__(self, rho_name, rho_label, rho_title, jet_name, occ_corr):
        self.fRhoObjectName = rho_name[rho_name.find("_") + 1:rho_name.rfind("_")]
        self.fRhoName = rho_name
        self.fRhoLabel = rho_label
        self.fRhoTitle = rho_title
        self.fJetName = jet_name
        self.fHasOccCorrFactor = occ_corr
        self.GenerateHash()
        self.GenerateName()

    def GenerateName(self):
        self.fShortName = CleanUpHistogramName("{}_{}_{}".format(self.fRhoName, self.fRhoLabel, self.fJetName))

    def GenerateHash(self):
        rho_name = self.fRhoName.replace("Gen", "")
        self.fHash = hash((rho_name, self.fRhoLabel, self.fJetName))

    def __hash__(self):
        return self.fHash

    def __eq__(self, other):
        if self.fHash == other.fHash:
            return True
        else:
            return False

    def __gt__(self, other):
        if self.fHash > other.fHash:
            return True
        else:
            return False

    def __ne__(self, other):
        return (not self == other)

    def __ge__(self, other):
        return (self > other or self == other)

    def __lt__(self, other):
        return (not self >= other)

    def __le__(self, other):
        return (not self > other)

    def GetRhoVsCentName(self):
        return "{}/fHist{}RhoVsCent".format(self.fRhoName, self.fRhoLabel)

    def GetRhoVsLeadJetPtName(self):
        return "{}/{}_fHist{}RhoVsLeadJetPt".format(self.fRhoName, self.fJetName, self.fRhoLabel)

    def GetRhoVsTrackPtName(self):
        return "{}/fHist{}RhoVsLeadTrackPt".format(self.fRhoName, self.fRhoLabel)

    def GetLeadJetPtVsCentName(self):
        return "{}/{}_fHistLeadJetPtVsCent".format(self.fRhoName, self.fJetName)

    def GetLeadTrackPtVsCentName(self):
        return "{}/fHistLeadTrackPtVsCent".format(self.fRhoName)

    def GetOccCorrVsCentName(self):
        if self.fHasOccCorrFactor: return "{}/fHistOccCorrvsCent".format(self.fRhoName)
        else: return ""

    def GetRCDeltaPtVsCentName(self, jet_type, jet_radius):
        hist_name = "RCDeltaPt"
        return self.GetDeltaPtVsCentName(jet_type, jet_radius, hist_name)

    def GetRCPerpDeltaPtVsCentName(self, jet_type, jet_radius):
        hist_name = "RCPerpDeltaPt"
        return self.GetDeltaPtVsCentName(jet_type, jet_radius, hist_name)

    def GetRCExclLeadJetDeltaPtVsCentName(self, jet_type, jet_radius):
        hist_name = "RCExclLeadJetDeltaPt"
        return self.GetDeltaPtVsCentName(jet_type, jet_radius, hist_name)

    def GetDeltaPtVsCentName(self, jet_type, jet_radius, hist_name):
        if "Gen" in self.fRhoName: gen = True
        else: gen = False
        if gen:
            task_name = "AliAnalysisTaskJetUEStudies_gen_histos"
            jet_name = "Jet_AKT{jet_type}{jet_radius}_mcparticles_pT0000_pt_scheme".format(jet_type=jet_type, jet_radius=jet_radius)
        else:
            task_name = "AliAnalysisTaskJetUEStudies_histos"
            if jet_type == "Full":
                jet_name = "Jet_AKTFull{jet_radius}_tracks_pT0150_caloClusters_ET0300_pt_scheme".format(jet_radius=jet_radius)
            else:
                jet_name = "Jet_AKTCharged{jet_radius}_tracks_pT0150_pt_scheme".format(jet_radius=jet_radius)
        return "{task_name}/{jet_name}/{rho_object_name}/fHist{rho_label}{hist_name}VsCent".format(task_name=task_name, jet_name=jet_name, rho_object_name=self.fRhoObjectName, rho_label=self.fRhoLabel, hist_name=hist_name)

    def Print(self):
        print(self.fHash, self.fJetName, self.fRhoLabel, self.fRhoName, self.fRhoTitle)

class UEHistograms:
    def __init__(self, config, meson_name, jet_type, jet_radius, gen, cent):
        self.fYamlConfig = config
        self.fMesonName = meson_name
        self.fJetType = jet_type
        self.fJetRadius = jet_radius
        self.fGeneratorLevel = gen
        self.fCentrality = cent
        self.fPtBins = numpy.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 80, 100], dtype=numpy.float32)
        self.fRhoBins = numpy.array(list(DMesonJetUtils.frange(0, 15, 0.5, True)), dtype=numpy.float32)
        self.fCentBins = numpy.array(list(DMesonJetUtils.frange(0, 100, 5, True)), dtype=numpy.float32)
        self.fOccCorrBins = numpy.array(list(DMesonJetUtils.frange(0, 1.01, 0.005, True)), dtype=numpy.float32)
        self.fCoarseCentBins = numpy.array([0, 10, 30, 50, 90], dtype=numpy.float32)
        self.fFiles = []
        self.fTitle = config["collision_system"]
        self.fFileNames = []
        self.fRhoGenLevDefinitions = []
        self.fRhoDetLevDefinitions = []
        self.fJetNames = []
        self.fCanvases = []
        self.fHistograms = dict()

        if not self.fYamlConfig["monte_carlo"]: self.fGeneratorLevel = False

        self.GenerateStandardRhoDefinitions()
        self.GenerateDMesonTaskName()

    def GenerateStandardRhoDefinitions(self):
        self.fRhoDetLevDefinitions.append(RhoDefinition("AliAnalysisTaskRhoDev_Rho_histos", "", "CMS Method", "Signal", True))
        # self.fRhoDetLevDefinitions.append(RhoDefinition("AliAnalysisTaskRhoDev_RhoExclLeadJets_histos", "", "CMS Method, excl. lead. jets", "Signal", True))
        self.fRhoDetLevDefinitions.append(RhoDefinition("AliAnalysisTaskRhoTransDev_RhoTrans_histos", "", "Trans Plane", "Signal", False))
        # self.fRhoDetLevDefinitions.append(RhoDefinition("AliAnalysisTaskRhoTransDev_RhoTrans_histos", "B2B", "Trans Plane, back-to-back", "Signal", False))

        self.fDefaultRhoDetLevDefinition = RhoDefinition("AliAnalysisTaskRhoDev_Rho_histos", "", "CMS Method", "Signal", True)

        self.fRhoGenLevDefinitions.append(RhoDefinition("AliAnalysisTaskRhoDev_RhoGen_histos", "", "CMS Method (gen. lev.)", "Signal", True))
        # self.fRhoGenLevDefinitions.append(RhoDefinition("AliAnalysisTaskRhoDev_RhoExclLeadJetsGen_histos", "", "CMS Method, excl. lead. jets", "Signal", True))
        self.fRhoGenLevDefinitions.append(RhoDefinition("AliAnalysisTaskRhoTransDev_RhoTransGen_histos", "", "Trans Plane (gen. lev.)", "Signal", False))
        # self.fRhoGenLevDefinitions.append(RhoDefinition("AliAnalysisTaskRhoTransDev_RhoTransGen_histos", "B2B", "Trans Plane, back-to-back (gen. lev.)", "Signal", False))

        self.fDefaultRhoGenLevDefinition = RhoDefinition("AliAnalysisTaskRhoDev_RhoGen_histos", "", "CMS Method", "Background", True)

    def GenerateDMesonTaskName(self):
        self.fTrigger = self.fYamlConfig["analysis"][0]["trigger"][0]
        if self.fTrigger:
            self.fDMesonTaskName = "{}_{}".format(self.fYamlConfig["task_name"], self.fTrigger)
        else:
            self.fDMesonTaskName = self.fYamlConfig["task_name"]

    def OpenFiles(self):
        path = "{0}/{1}".format(self.fYamlConfig["input_path"], self.fYamlConfig["train"])
        print("Looking for file {0} in path {1}".format(self.fYamlConfig["file_name"], path))
        self.fFileNames.extend(DMesonJetUtils.find_file(path, self.fYamlConfig["file_name"]))
        for fname in self.fFileNames:
            f = ROOT.TFile(fname, "read")
            self.fFiles.append(f)

    def LoadHistograms(self, hname, title, varname, units, binsX, binsY):
        histograms = dict()

        hist_2d = DMesonJetUtils.GetObjectAndMerge(self.fFiles, hname)
        if hist_2d.GetXaxis().GetTitle():
            hist_2d.GetXaxis().SetTitle(hist_2d.GetXaxis().GetTitle().replace("Centrality", "Multiplicity"))
        else:
            hist_2d.GetXaxis().SetTitle("Multiplicity (%)")
        hist_2d.GetZaxis().SetTitle("counts")
        hist_2d.SetName(CleanUpHistogramName(hname))
        histograms[hist_2d.GetName()] = hist_2d

        h_rebin_name = "{}_Rebinned".format(hist_2d.GetName())
        hist_2d_rebin = DMesonJetUtils.Rebin2D_fromBins(hist_2d, h_rebin_name, len(binsX) - 1, binsX, len(binsY) - 1, binsY)
        histograms[hist_2d_rebin.GetName()] = hist_2d_rebin

        prof = hist_2d_rebin.ProfileX("{}_Profile".format(hist_2d.GetName()), 1, -1, "i")
        if units: prof.GetYaxis().SetTitle("<{}> ({})".format(varname, units))
        else: prof.GetYaxis().SetTitle("<{}>".format(varname))
        prof.GetYaxis().SetTitleOffset(1)
        prof.SetTitle(title)
        histograms[prof.GetName()] = prof

        std_dev = GenerateStdDevProfile(hist_2d_rebin, "{}_StdDev".format(hist_2d.GetName()), "Std Dev {}".format(title), "#sigma({varname}) / <{varname}>".format(varname=varname))
        histograms[std_dev.GetName()] = std_dev

        return histograms

    def LoadDetLevHistograms(self):
        self.LoadRhoVsDPt(self.fMesonName)
        self.LoadRhoVsJetPt(self.fMesonName)
        for rhoDef in self.fRhoDetLevDefinitions:
            if rhoDef.fHasOccCorrFactor: self.LoadOccCorrVsCent(rhoDef)
            self.LoadDeltaPt(rhoDef)

            h_new_name = self.LoadRhoVsTrackPt(rhoDef)
            self.LoadRhoVsLeadJetPt(rhoDef)

            rho_vs_trackpt_2d = self.fHistograms[h_new_name][h_new_name]
            rho = Projections(rho_vs_trackpt_2d, None, 5, 0, 25)
            rho.SetTitle(rhoDef.fRhoTitle)
            rho.GetYaxis().SetTitle("counts")
            hrhoname = "RhoDistribution_{}".format(rhoDef.fShortName)
            rho.SetName(hrhoname)
            self.fHistograms[hrhoname] = rho
        if self.fCentrality: self.LoadCentDetLevHistograms()

    def LoadCentDetLevHistograms(self):
        self.LoadLeadTrackPtVsCent(self.fDefaultRhoDetLevDefinition)
        self.LoadLeadJetPtVsCent(self.fDefaultRhoDetLevDefinition)
        for rhoDef in self.fRhoDetLevDefinitions: self.LoadRhoVsCent(rhoDef)

    def LoadGenLevHistograms(self):
        self.LoadRhoVsDPt("{}_MCTruth".format(self.fMesonName))
        self.LoadRhoVsJetPt("{}_MCTruth".format(self.fMesonName))
        for rhoDef in self.fRhoGenLevDefinitions:
            if rhoDef.fHasOccCorrFactor: self.LoadOccCorrVsCent(rhoDef)
            self.LoadDeltaPt(rhoDef)

            h_new_name = self.LoadRhoVsTrackPt(rhoDef)
            self.LoadRhoVsLeadJetPt(rhoDef)

            rho_vs_trackpt_2d = self.fHistograms[h_new_name][h_new_name]
            rho = Projections(rho_vs_trackpt_2d, None, 5, 0, 25)
            rho.GetYaxis().SetTitle("counts")
            rho.SetTitle(rhoDef.fRhoTitle)
            hrhoname = "RhoDistribution_{}".format(rhoDef.fShortName)
            rho.SetName(hrhoname)
            self.fHistograms[hrhoname] = rho
        if self.fCentrality: self.LoadCentGenLevHistograms()

    def LoadCentGenLevHistograms(self):
        self.LoadLeadTrackPtVsCent(self.fDefaultRhoGenLevDefinition)
        self.LoadLeadJetPtVsCent(self.fDefaultRhoGenLevDefinition)
        for rhoDef in self.fRhoGenLevDefinitions: self.LoadRhoVsCent(rhoDef)

    def LoadAllHistograms(self):
        self.LoadDetLevHistograms()
        if self.fGeneratorLevel: self.LoadGenLevHistograms()

    def LoadDeltaPt(self, rho_definition):
        RC_names = ["RCDeltaPt", "RCPerpDeltaPt", "RCExclLeadJetDeltaPt"]
        RC_titles = ["Random Cones", "RC perp. lead. jet", "RC excl. lead. jet"]
        for RCname, RCtitle in zip(RC_names, RC_titles):
            deltaPt_list = dict()
            hname = rho_definition.GetDeltaPtVsCentName(self.fJetType, self.fJetRadius, RCname)
            h_new_name = CleanUpHistogramName(hname)

            h_2d = DMesonJetUtils.GetObjectAndMerge(self.fFiles, hname)
            h_2d.GetXaxis().SetTitle(h_2d.GetXaxis().GetTitle().replace("Centrality", "Multiplicity"))
            deltaPt_list[h_new_name] = h_2d

            deltaPt = Projections(h_2d, None, 0)
            deltaPt.Scale(1. / deltaPt.Integral(), "width")
            deltaPt.GetYaxis().SetTitle("probability density")
            deltaPt.SetTitle(RCtitle)
            h_deltapt_name = "{}Distribution_{}".format(RCname, rho_definition.fShortName)
            deltaPt.SetName(h_deltapt_name)
            deltaPt_list[h_deltapt_name] = deltaPt

            self.fHistograms[h_new_name] = deltaPt_list

    def LoadOccCorrVsCent(self, rho_definition):
        hname = rho_definition.GetOccCorrVsCentName()
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.LoadHistograms(hname, rho_definition.fRhoTitle, "C", "", self.fCentBins, self.fOccCorrBins)

        occ_corr_vs_cent_2d = histograms[h_new_name]
        occ_corr_vs_cent_2d.GetXaxis().SetTitle("Multiplicity (%)")
        occ_corr_vs_cent_2d.GetYaxis().SetTitle("#it{C}")
        occ_corr = Projections(occ_corr_vs_cent_2d, None, 0)
        h_occ_corr_name = "OccCorrFactDistribution_{}".format(rho_definition.fShortName)
        occ_corr.SetName(h_occ_corr_name)
        occ_corr.GetXaxis().SetTitle("#it{C}")
        occ_corr.GetYaxis().SetTitle("#counts")
        histograms[h_occ_corr_name] = occ_corr

        if h_new_name in self.fHistograms: raise OverwriteError(self.fHistograms, h_new_name)
        self.fHistograms[h_new_name] = histograms

        return h_new_name

    def LoadRhoVsCent(self, rho_definition):
        hname = rho_definition.GetRhoVsCentName()
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.LoadHistograms(hname, rho_definition.fRhoTitle, "#rho", "GeV/#it{c}", self.fCentBins, self.fRhoBins)

        rho_vs_cent_2d = histograms[h_new_name]
        histograms["RhoCentralityBins"] = Projections(rho_vs_cent_2d, self.fCoarseCentBins, 5, 0, 25)
        for h in histograms["RhoCentralityBins"]: h.GetYaxis().SetTitle("counts")

        if h_new_name in self.fHistograms: raise OverwriteError(self.fHistograms, h_new_name)
        self.fHistograms[h_new_name] = histograms

        return h_new_name

    def LoadRhoVsLeadJetPt(self, rho_definition):
        hname = rho_definition.GetRhoVsLeadJetPtName()
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.LoadHistograms(hname, rho_definition.fRhoTitle, "#rho", "GeV/#it{c}", self.fPtBins, self.fRhoBins)

        if h_new_name in self.fHistograms: raise OverwriteError(self.fHistograms, h_new_name)
        self.fHistograms[h_new_name] = histograms

        return h_new_name

    def LoadRhoVsTrackPt(self, rho_definition):
        hname = rho_definition.GetRhoVsTrackPtName()
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.LoadHistograms(hname, rho_definition.fRhoTitle, "#rho", "GeV/#it{c}", self.fPtBins, self.fRhoBins)

        if h_new_name in self.fHistograms: raise OverwriteError(self.fHistograms, h_new_name)
        self.fHistograms[h_new_name] = histograms

        return h_new_name

    def LoadLeadJetPtVsCent(self, rho_definition):
        hname = rho_definition.GetLeadJetPtVsCentName()
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.LoadHistograms(hname, rho_definition.fRhoTitle, "#it{p}_{T,jet}^{lead}", "GeV/#it{c}", self.fCentBins, self.fPtBins)

        if h_new_name in self.fHistograms: raise OverwriteError(self.fHistograms, h_new_name)
        self.fHistograms[h_new_name] = histograms

        return h_new_name

    def LoadLeadTrackPtVsCent(self, rho_definition):
        hname = rho_definition.GetLeadTrackPtVsCentName()
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.LoadHistograms(hname, rho_definition.fRhoTitle, "#it{p}_{T,track}^{lead}", "GeV/#it{c}", self.fCentBins, self.fPtBins)

        if h_new_name in self.fHistograms: raise OverwriteError(self.fHistograms, h_new_name)
        self.fHistograms[h_new_name] = histograms

        return h_new_name

    def LoadRhoVsDPt(self, meson_name):
        hname = "{taskName}_histos/histos{taskName}/{meson_name}/Jet_AKT{jet_type}{jet_radius}_pt_scheme/fHistRhoVsLeadDPt".format(taskName=self.fDMesonTaskName, meson_name=meson_name, jet_type=self.fJetType, jet_radius=self.fJetRadius)
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.LoadHistograms(hname, "#rho vs. #it{p}_{T,D}", "#rho", "GeV/#it{c}", self.fPtBins, self.fRhoBins)

        if h_new_name in self.fHistograms: raise OverwriteError(self.fHistograms, h_new_name)
        self.fHistograms[h_new_name] = histograms

        return h_new_name

    def LoadRhoVsJetPt(self, meson_name):
        hname = "{taskName}_histos/histos{taskName}/{meson_name}/Jet_AKT{jet_type}{jet_radius}_pt_scheme/fHistRhoVsLeadJetPt".format(taskName=self.fDMesonTaskName, meson_name=meson_name, jet_type=self.fJetType, jet_radius=self.fJetRadius)
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.LoadHistograms(hname, "#rho vs. #it{p}_{T,jet}", "#rho", "GeV/#it{c}", self.fPtBins, self.fRhoBins)

        if h_new_name in self.fHistograms: raise OverwriteError(self.fHistograms, h_new_name)
        self.fHistograms[h_new_name] = histograms

        return h_new_name

    def PlotDetLevHistograms(self):
        self.PlotDeltaPt(self.fRhoDetLevDefinitions)
        self.PlotOccCorrVsCent(self.fRhoDetLevDefinitions)
        self.PlotRhoVsDPt(self.fMesonName)
        self.PlotRhoVsJetPt(self.fMesonName)
        self.PlotRhoVsTrackPt(self.fRhoDetLevDefinitions, "DetLev")
        self.PlotRhoVsLeadJetPt(self.fRhoDetLevDefinitions, "DetLev")
        self.PlotRhoDistributions(self.fRhoDetLevDefinitions, "DetLev")
        if self.fCentrality:
            self.PlotLeadTrackPtVsCent(self.fDefaultRhoDetLevDefinition)
            self.PlotLeadJetPtVsCent(self.fDefaultRhoDetLevDefinition)
            self.PlotRhoVsCent(self.fRhoDetLevDefinitions, "DetLev")

    def PlotGenLevHistograms(self):
        self.PlotDeltaPt(self.fRhoGenLevDefinitions)
        self.PlotOccCorrVsCent(self.fRhoGenLevDefinitions)
        self.PlotRhoVsDPt("{}_MCTruth".format(self.fMesonName))
        self.PlotRhoVsJetPt("{}_MCTruth".format(self.fMesonName))
        self.PlotRhoVsTrackPt(self.fRhoGenLevDefinitions, "GenLev")
        self.PlotRhoVsLeadJetPt(self.fRhoGenLevDefinitions, "GenLev")
        self.PlotRhoDistributions(self.fRhoGenLevDefinitions, "GenLev")
        if self.fCentrality:
            self.PlotLeadTrackPtVsCent(self.fDefaultRhoGenLevDefinition)
            self.PlotLeadJetPtVsCent(self.fDefaultRhoGenLevDefinition)
            self.PlotRhoVsCent(self.fRhoGenLevDefinitions, "GenLev")

    def PlotAllHistograms(self):
        self.PlotDetLevHistograms()
        if self.fGeneratorLevel: self.PlotGenLevHistograms()

    def PlotRhoDistributions(self, rho_definitions, suffix):
        hListRho = []
        for rho_def in rho_definitions:
            hname = "RhoDistribution_{}".format(rho_def.fShortName)
            h = self.fHistograms[hname]
            hListRho.append(h)
        self.PlotMultiple("RhoDistribution{}".format(suffix), hListRho, True, True, "stat")

    def PlotOccCorrVsCent(self, rho_definitions):
        for rho_definition in rho_definitions:
            if not rho_definition.fHasOccCorrFactor: continue
            hname = rho_definition.GetOccCorrVsCentName()
            h_new_name = CleanUpHistogramName(hname)
            histograms = self.fHistograms[h_new_name]

            h = histograms[h_new_name]
            self.PlotSingle(h)
            h = histograms["OccCorrFactDistribution_{}".format(rho_definition.fShortName)]
            self.PlotSingle(h)
            h = histograms["{}_Profile".format(h_new_name)]
            self.PlotSingle(h)
            h = histograms["{}_StdDev".format(h_new_name)]
            self.PlotSingle(h)

    def PlotRhoVsCent(self, rho_definitions, suffix):
        hListMean = []
        hListStdDev = []
        for rho_def in rho_definitions:
            hname = rho_def.GetRhoVsCentName()
            h_new_name = CleanUpHistogramName(hname)
            histograms = self.fHistograms[h_new_name]
            h = histograms["{}_Profile".format(h_new_name)]
            hListMean.append(h)
            h = histograms["{}_StdDev".format(h_new_name)]
            hListStdDev.append(h)

            self.PlotMultiple("{}_DistributionCentralityBins".format(h_new_name), histograms["RhoCentralityBins"], True, True, "stat")

            h_2d = histograms[h_new_name]
            self.PlotSingle(h_2d)

        self.PlotMultiple("MeanRhoVsCent{}".format(suffix), hListMean, True)
        self.PlotMultiple("StdDevRhoVsCent{}".format(suffix), hListStdDev, False)

    def PlotRhoVsLeadJetPt(self, rho_definitions, suffix):
        hListMean = []
        hListStdDev = []
        for rho_def in rho_definitions:
            hname = rho_def.GetRhoVsLeadJetPtName()
            h_new_name = CleanUpHistogramName(hname)
            histograms = self.fHistograms[h_new_name]
            h = histograms["{}_Profile".format(h_new_name)]
            hListMean.append(h)
            h = histograms["{}_StdDev".format(h_new_name)]
            hListStdDev.append(h)

            h_2d = histograms[h_new_name]
            self.PlotSingle(h_2d)
        self.PlotMultiple("MeanRhoVsLeadJetPt{}".format(suffix), hListMean, True)
        self.PlotMultiple("StdDevRhoVsLeadJetPt{}".format(suffix), hListStdDev, False)

    def PlotRhoVsTrackPt(self, rho_definitions, suffix):
        hListMean = []
        hListStdDev = []
        for rho_def in rho_definitions:
            hname = rho_def.GetRhoVsTrackPtName()
            h_new_name = CleanUpHistogramName(hname)
            histograms = self.fHistograms[h_new_name]
            h = histograms["{}_Profile".format(h_new_name)]
            hListMean.append(h)
            h = histograms["{}_StdDev".format(h_new_name)]
            hListStdDev.append(h)

            h_2d = histograms[h_new_name]
            self.PlotSingle(h_2d)
        self.PlotMultiple("MeanRhoVsLeadTrackPt{}".format(suffix), hListMean, True)
        self.PlotMultiple("StdDevRhoVsLeadTrackPt{}".format(suffix), hListStdDev, False)

    def PlotLeadJetPtVsCent(self, rho_def):
        hname = rho_def.GetLeadJetPtVsCentName()
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.fHistograms[h_new_name]
        h = histograms["{}_Profile".format(h_new_name)]
        self.PlotSingle(h)
        h = histograms["{}_StdDev".format(h_new_name)]
        self.PlotSingle(h)

        h_2d = histograms[h_new_name]
        self.PlotSingle(h_2d)

    def PlotLeadTrackPtVsCent(self, rho_def):
        hname = rho_def.GetLeadTrackPtVsCentName()
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.fHistograms[h_new_name]
        h = histograms["{}_Profile".format(h_new_name)]
        self.PlotSingle(h)
        h = histograms["{}_StdDev".format(h_new_name)]
        self.PlotSingle(h)

        h_2d = histograms[h_new_name]
        self.PlotSingle(h_2d)

    def PlotRhoVsDPt(self, meson_name):
        hname = "{taskName}_histos/histos{taskName}/{meson_name}/Jet_AKT{jet_type}{jet_radius}_pt_scheme/fHistRhoVsLeadDPt".format(taskName=self.fDMesonTaskName, meson_name=meson_name, jet_type=self.fJetType, jet_radius=self.fJetRadius)
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.fHistograms[h_new_name]
        h = histograms["{}_Profile".format(h_new_name)]
        self.PlotSingle(h)
        h = histograms["{}_StdDev".format(h_new_name)]
        self.PlotSingle(h)

        h_2d = histograms[h_new_name]
        self.PlotSingle(h_2d)

    def PlotRhoVsJetPt(self, meson_name):
        hname = "{taskName}_histos/histos{taskName}/{meson_name}/Jet_AKT{jet_type}{jet_radius}_pt_scheme/fHistRhoVsLeadJetPt".format(taskName=self.fDMesonTaskName, meson_name=meson_name, jet_type=self.fJetType, jet_radius=self.fJetRadius)
        h_new_name = CleanUpHistogramName(hname)
        histograms = self.fHistograms[h_new_name]
        h = histograms["{}_Profile".format(h_new_name)]
        self.PlotSingle(h)
        h = histograms["{}_StdDev".format(h_new_name)]
        self.PlotSingle(h)

        h_2d = histograms[h_new_name]
        self.PlotSingle(h_2d)

    def PlotDeltaPt(self, rho_definitions):
        for rho_def in rho_definitions:
            RC_names = ["RCDeltaPt", "RCPerpDeltaPt", "RCExclLeadJetDeltaPt"]
            hlist = []
            for RCname in RC_names:
                hname = rho_def.GetDeltaPtVsCentName(self.fJetType, self.fJetRadius, RCname)
                h_new_name = CleanUpHistogramName(hname)
                h_deltapt_name = "{}Distribution_{}".format(RCname, rho_def.fShortName)
                deltaPt = self.fHistograms[h_new_name][h_deltapt_name]
                hlist.append(deltaPt)
            comp = self.PlotMultiple("RCDeltaPt_{}".format(rho_def.fShortName), hlist, True, True, "stat")
            comp.fLegendSpectra.SetX1(0.08)
            comp.fLegendSpectra.SetX2(0.47)
            for r in comp.fResults: globalList.append(r)

    def PlotSingle(self, h):
        cname = h.GetName()
        c = ROOT.TCanvas(cname, cname)
        globalList.append(c)
        c.cd()
        if isinstance(h, ROOT.TH2):
            c.SetLogz()
            h.Draw("colz")
        else:
            h.Draw("")
        self.fCanvases.append(c)
        return c

    def PlotMultiple(self, name, hlist, errors, logY=False, leg=True):
        if len(hlist) < 2:
            print("Not enough histograms for {}".format(name))
            print(hlist)
            return None
        else:
            print("Plotting {}".format(name))
        comp = DMesonJetCompare.DMesonJetCompare(name)
        comp.fDoSpectrumLegend = leg
        if leg == "stat":
            xaxis = hlist[0].GetXaxis().GetTitle()
            comp.fUnits = xaxis[xaxis.rfind("(") + 1:xaxis.rfind(")")]
        if logY:
            comp.fDoSpectraPlot = "logy"
        else:
            comp.fDoSpectraPlot = "lineary"
        comp.fDoRatioPlot = False
        comp.fX1LegSpectrum = 0.25
        if not errors:
            comp.fOptSpectrumBaseline = "hist"
            comp.fOptSpectrum = "hist"
        r = comp.CompareSpectra(hlist[0], hlist[1:])
        globalList.extend(r)
        self.fCanvases.append(comp.fCanvasSpectra)
        return comp

    def SaveAllPlots(self):
        output_path = "{0}/{1}/ue_studies".format(self.fYamlConfig["input_path"], self.fYamlConfig["train"])
        if not os.path.isdir(output_path): os.makedirs(output_path)
        for c in self.fCanvases: c.SaveAs("{}/{}.pdf".format(output_path, c.GetName()))

    def CompareWith(self, ue2):
        self.fCompareUE = ue2

        self.fCompareRhoDetLevDefinitions1 = set(self.fRhoDetLevDefinitions)
        if self.fGeneratorLevel: self.fCompareRhoGenLevDefinitions1 = set(self.fRhoGenLevDefinitions)
        else: self.fCompareRhoGenLevDefinitions1 = set()

        self.fCompareRhoDetLevDefinitions2 = set(ue2.fRhoDetLevDefinitions)
        if ue2.fGeneratorLevel: self.fCompareRhoGenLevDefinitions2 = set(ue2.fRhoGenLevDefinitions)
        else: self.fCompareRhoGenLevDefinitions2 = set()

        rho_definitions1 = self.fCompareRhoDetLevDefinitions1 | self.fCompareRhoGenLevDefinitions1
        rho_definitions2 = self.fCompareRhoDetLevDefinitions2 | self.fCompareRhoGenLevDefinitions2

        self.fCompareRhoDefinitions = rho_definitions1 & rho_definitions2

        self.CompareGeneric("RCDeltaPtVsCent", "", True)
        self.CompareGeneric("RhoVsLeadDPt", self.fMesonName)
        self.CompareGeneric("RhoVsLeadJetPt", self.fMesonName)
        self.CompareGeneric("RhoVsTrackPt")
        self.CompareGeneric("RhoVsLeadJetPt")
        self.CompareGeneric("RhoDistr", "", True)
        if self.fCentrality and ue2.fCentrality:
            self.CompareGeneric("LeadTrackPtVsCent")
            self.CompareGeneric("LeadJetPtVsCent")
            self.CompareGeneric("RhoVsCent")

    def CompareGeneric(self, func, meson_name="", logY=False):
        if meson_name: meson_name_genlev = "{}_MCTruth".format(meson_name)
        else: meson_name_genlev = ""
        for rho_def_comp in self.fCompareRhoDefinitions:
            if meson_name and rho_def_comp != self.fDefaultRhoDetLevDefinition: continue
            hListMean = []
            hListStdDev = []

            def AddRhoDefinition(rho_def, title, ue, meson_name_internal):
                if not rho_def_comp == rho_def: return
                if func == "RhoDistr":
                    hname = "RhoDistribution_{}".format(rho_def.fShortName)
                    h = ue.fHistograms[hname]
                    h_copy = h.Clone()
                    h_copy.SetTitle(title)
                    globalList.append(h_copy)
                    hListMean.append(h_copy)
                else:
                    if meson_name_internal:
                        hname = "{taskName}_histos/histos{taskName}/{meson_name}/Jet_AKT{jet_type}{jet_radius}_pt_scheme/fHist{hname}".format(taskName=ue.fDMesonTaskName, meson_name=meson_name_internal, jet_type=ue.fJetType, jet_radius=ue.fJetRadius, hname=func)
                    else:
                        if "DeltaPt" in func:
                            hname = getattr(rho_def, "Get{}Name".format(func))(self.fJetType, self.fJetRadius)
                        else:
                            hname = getattr(rho_def, "Get{}Name".format(func))()
                    if "DeltaPt" in func:
                        h_new_name = CleanUpHistogramName(hname)
                        histograms = ue.fHistograms[h_new_name]
                        RCname = "RCDeltaPt"
                        h = histograms["{}Distribution_{}".format(RCname, rho_def.fShortName)]
                        h_copy = h.Clone()
                        h_copy.SetTitle(title)
                        globalList.append(h_copy)
                        hListMean.append(h_copy)
                    else:
                        h_new_name = CleanUpHistogramName(hname)
                        histograms = ue.fHistograms[h_new_name]
                        h = histograms["{}_Profile".format(h_new_name)]
                        h_copy = h.Clone()
                        h_copy.SetTitle(title)
                        globalList.append(h_copy)
                        hListMean.append(h_copy)
                        h = histograms["{}_StdDev".format(h_new_name)]
                        h_copy = h.Clone()
                        h_copy.SetTitle(title)
                        globalList.append(h_copy)
                        hListStdDev.append(h_copy)

            for rho_def in self.fCompareRhoDetLevDefinitions1:
                AddRhoDefinition(rho_def, self.fTitle, self, meson_name)

            for rho_def in self.fCompareRhoGenLevDefinitions1:
                AddRhoDefinition(rho_def, "{} (gen. lev.)".format(self.fTitle), self, meson_name_genlev)

            for rho_def in self.fCompareRhoDetLevDefinitions2:
                AddRhoDefinition(rho_def, self.fCompareUE.fTitle, self.fCompareUE, meson_name)

            for rho_def in self.fCompareRhoGenLevDefinitions2:
                AddRhoDefinition(rho_def, "{} (gen. lev.)".format(self.fCompareUE.fTitle), self.fCompareUE, meson_name_genlev)

            if len(hListMean) > 1: self.PlotMultiple("{}_{}{}_Mean{}_{}".format(self.fYamlConfig["name"], self.fCompareUE.fYamlConfig["name"], func, meson_name, rho_def_comp.fShortName), hListMean, True, logY)
            if len(hListStdDev) > 1: self.PlotMultiple("{}_{}{}_StdDev{}_{}".format(self.fYamlConfig["name"], self.fCompareUE.fYamlConfig["name"], func, meson_name, rho_def_comp.fShortName), hListStdDev, False, logY)

def main(config, meson_name, jet_type, jet_radius, gen, cent, config2):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if config2:
        ue_studies1 = UEHistograms(config, meson_name, jet_type, jet_radius, gen, cent)
        ue_studies1.OpenFiles()
        ue_studies1.LoadAllHistograms()

        ue_studies2 = UEHistograms(config2, meson_name, jet_type, jet_radius, gen, cent)
        ue_studies2.OpenFiles()
        ue_studies2.LoadAllHistograms()

        ue_studies1.CompareWith(ue_studies2)
        ue_studies1.SaveAllPlots()

        globalList.append(ue_studies1)
        globalList.append(ue_studies2)
    else:
        ue_studies = UEHistograms(config, meson_name, jet_type, jet_radius, gen, cent)
        ue_studies.OpenFiles()
        ue_studies.LoadAllHistograms()
        ue_studies.PlotAllHistograms()
        ue_studies.SaveAllPlots()
        globalList.append(ue_studies)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Underlying event studies.')
    parser.add_argument('yaml', metavar='config.yaml')
    parser.add_argument('--meson', metavar='MESON',
                        default="D0")
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    parser.add_argument('--cent', action='store_const',
                        default=False, const=True,
                        help='Centrality available.')
    parser.add_argument('--gen', action='store_const',
                        default=False, const=True,
                        help='Generator level.')
    parser.add_argument('--compare', metavar='config2.yaml',
                        default=None)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    if args.compare:
        f = open(args.compare, 'r')
        config2 = yaml.load(f)
        f.close()
    else:
        config2 = None

    main(config, args.meson, args.jet_type, args.jet_radius, args.gen, args.cent, config2)

    IPython.embed()
