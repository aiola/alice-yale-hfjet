#!/usr/bin/env python
# python script to do extract B feed down correction factors

import argparse
import yaml
import IPython
import ROOT
import DMesonJetCompare
import DMesonJetUtils

globalList = []

class SpectrumDef:
    def __init__(self, _name, _title, _styles):
        self.name = _name
        self.title = _title
        self.styles = _styles

class SpectraSet:
    def __init__(self, _name, _title):
        self.name = _name
        self.title = _title
        self.spectra = []

    def add(self, s):
        self.spectra.append(s)

class QuarkSetting:
    def __init__(self, _name):
        self.name = _name
        self.histos = dict()

def main(charm_ts, beauty_ts, jet_type, jet_radius, data, pythia):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    generators = ["powheg"]
    for gen in generators:
        if data:
            data_comparison_for_generator(gen, charm_ts, beauty_ts, jet_type, jet_radius, data)
        elif pythia:
            pythia_comparison_for_generator(gen, charm_ts, beauty_ts, jet_type, jet_radius, pythia)
        else:
            feed_down_analysis_for_generator(gen, charm_ts, beauty_ts, jet_type, jet_radius)

def generate_correction_factors(quarks, jet_type, jet_radius):
    for quark in quarks.itervalues():
        quark.histogram = GetSpectrum(quark.file, "D0_MCTruth", jet_type, jet_radius, "JetPtDPtSpectrum")
        quark.histogram.SetName("{0}_{1}".format(quark.name, quark.histogram.GetName()))
    denominator = quarks["beauty"].histogram.Clone("corr_fact_den")
    denominator.Add(quarks["charm"].histogram)
    numerator = quarks["charm"].histogram.Clone("corr_fact_num")
    corr_fact = numerator.Clone("corr_fact")
    corr_fact.Divide(denominator)
    corr_fact.GetZaxis().SetTitle("(c #rightarrow D^{0}) / (c,b #rightarrow D^{0})")
    corr_fact_err = ROOT.TH2D("corr_fact_err", "corr_fact_err", corr_fact.GetNbinsX(), corr_fact.GetXaxis().GetXbins().GetArray(), corr_fact.GetNbinsY(), corr_fact.GetYaxis().GetXbins().GetArray())
    corr_fact_err.GetXaxis().SetTitle(corr_fact.GetXaxis().GetTitle())
    corr_fact_err.GetYaxis().SetTitle(corr_fact.GetYaxis().GetTitle())
    corr_fact_err.GetZaxis().SetTitle("rel. stat. unc. of FD corr. factors")
    for xbin in range(0, corr_fact.GetNbinsX() + 2):
        for ybin in range(0, corr_fact.GetNbinsY() + 2):
            if corr_fact.GetBinContent(xbin, ybin) == 0: continue
            corr_fact_err.SetBinContent(xbin, ybin, corr_fact.GetBinError(xbin, ybin) / corr_fact.GetBinContent(xbin, ybin))
    return corr_fact, corr_fact_err

class DataSpectrumDef:
    def __init__(self, spectrumName, unfoldingMethod, unfoldingReg, unfoldingPrior):
        self.fSpectrumName = spectrumName
        self.fUnfoldingMethod = unfoldingMethod
        self.fUnfoldingReg = unfoldingReg
        self.fUnfoldingPrior = unfoldingPrior
        self.fCrossSection = 62.3  # mb CINT1
        self.fBranchingRatio = 0.0393  # D0->Kpi
        self.fAntiParticleNorm = 2.0

    def LoadSpectrum(self, file):
        spectrumList = file.Get(self.fSpectrumName)
        if not spectrumList:
            print("Could not get list {0} from file {1}".format(self.fSpectrumName, file.GetName()))
            exit(1)
        else:
            print("List {0} loaded from file {1}".format(self.fSpectrumName, file.GetName()))
        methodList = spectrumList.FindObject(self.fUnfoldingMethod)
        if not methodList:
            print("Could not get list {0} from list {1}".format(self.fUnfoldingMethod, self.fSpectrumName))
            exit(1)
        else:
            print("List {0} loaded from list {1}".format(self.fUnfoldingMethod, self.fSpectrumName))
        hname = "_".join([self.fSpectrumName, "UnfoldedSpectrum", self.fUnfoldingMethod, self.fUnfoldingReg, self.fUnfoldingPrior])
        h = methodList.FindObject(hname)
        if not h:
            print("Could not get histogram {0} from list {1}".format(hname, self.fUnfoldingMethod))
            exit(1)
        else:
            print("Histogram {0} loaded from list {1}".format(hname, self.fUnfoldingMethod))
        self.fHistogram = h.Clone("{0}_copy".format(h.GetName()))
        self.fEvents = spectrumList.FindObject("Events")
        if not h:
            print("Could not get histogram {0} from list {1}".format("Events", spectrumList.GetName()))
            exit(1)
        else:
            print("Histogram {0} loaded from list {1}".format("Events", spectrumList.GetName()))
        self.fNumberOfEvents = self.fEvents.GetBinContent(1)
        self.fNormalizedHistogram = self.fHistogram.Clone("{0}_Normalized".format(hname))
        self.fNormalizedHistogram.Scale(self.fCrossSection / (self.fNumberOfEvents * self.fBranchingRatio * self.fAntiParticleNorm), "width")

def GetTotalMCSpectrum(quarks, jet_type, jet_radius, spectrumName, title):
    res = None
    hname = ""
    for quark in quarks.itervalues():
        h = GetSpectrum(quark.file, "D0_MCTruth", jet_type, jet_radius, spectrumName)
        hname = h.GetName()
        h.SetName("{0}_{1}".format(quark.name, h.GetName()))
        h.SetTitle("{0} #rightarrow D^{{0}}".format(quark.name[0]))
        quark.histos[spectrumName] = h
        if res:
            res.Add(h)
        else:
            res = h.Clone(hname)
            res.SetTitle(title)
    res.Scale(1. / 2)  # particle/antiparticle normalization
    if res.GetDimension() == 1:
        res.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T} d#eta} [mb (GeV/#it{c})^{-1}]")
    else:
        res.GetZaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T} d#eta} [mb (GeV/#it{c})^{-1}]")
    return res

def data_comparison_for_generator(gen, charm_ts, beauty_ts, jet_type, jet_radius, data):
    rootPath = "/Volumes/DATA/ALICE/JetResults"

    quarks = dict()

    if charm_ts:
        charmQuark = QuarkSetting("charm")
        (charmQuark.ts, charmQuark.stage) = get_ts_stage(charm_ts)
        quarks["charm"] = charmQuark

    if beauty_ts:
        beautyQuark = QuarkSetting("beauty")
        (beautyQuark.ts, beautyQuark.stage) = get_ts_stage(beauty_ts)
        quarks["beauty"] = beautyQuark

    for quark in quarks.itervalues():
        if quark.stage >= 0:
            quark.path = "{0}/FastSim_{1}_{2}_{3}/stage_{4}/output".format(rootPath, gen, quark.name, quark.ts, quark.stage)
            quark.filename = "{0}/FastSimAnalysis_{1}_{2}_{3}.root".format(quark.path, gen, quark.name, quark.ts)
        else:
            quark.path = "{0}/FastSim_{1}_{2}_{3}/output".format(rootPath, gen, quark.name, quark.ts)
            quark.filename = "{0}/FastSimAnalysis_{1}_{2}.root".format(quark.path, gen, quark.name)
        quark.file = ROOT.TFile(quark.filename)
        if not quark.file or quark.file.IsZombie():
            print("Could not open file {0}".format(quark.filename))
            exit(1)

    dataFileName = "{0}/{1}/{1}.root".format(rootPath, data)
    dataFile = ROOT.TFile(dataFileName)
    if not dataFile or dataFile.IsZombie():
        print("Could not open file {0}".format(dataFileName))
        exit(1)
    else:
        print("File {0} successfully open".format(dataFileName))

    if gen == "powheg":
        genTitle = "POWHEG+PYTHIA6"
    else:
        genTitle = "MC"
    MCspectrum = GetTotalMCSpectrum(quarks, jet_type, jet_radius, "JetPtSpectrum_DPt_3", genTitle)

    spectra = []
    spectra.append(DataSpectrumDef("SideBand_DPt_30", "Bayes", "Reg4", "PriorResponseTruth"))

    cname = "_".join(["TheoryComparison", gen, jet_type, jet_radius])
    ratioAxis = "data / theory"
    histos = []
    for spectrum in spectra:
        spectrum.LoadSpectrum(dataFile)
        spectrum.fNormalizedHistogram.SetTitle("ALICE")
        histos.append(spectrum.fNormalizedHistogram)

    MCspectrumRebinned = MCspectrum.Rebin(histos[0].GetNbinsX(), "{0}_Rebinned".format(MCspectrum.GetName()), histos[0].GetXaxis().GetXbins().GetArray())
    MCspectrumRebinned.Scale(1., "width")
    globalList.append(MCspectrumRebinned)
    globalList.extend(histos)
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fYaxisRatio = ratioAxis
    r = comp.CompareSpectra(MCspectrumRebinned, histos)
    for obj in r:
        globalList.append(obj)
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{0}/{1}.pdf".format(rootPath, obj.GetName()))

def GetPYTHIASpectrum(pythiaFile, jet_type, jet_radius, spectrumName):
    hpath = "D0_Jet_AKT{jet_type}{jet_radius}_pt_scheme_{spectrum}/D0_Jet_AKT{jet_type}{jet_radius}_pt_scheme_{spectrum}_Truth".format(jet_type=jet_type, jet_radius=jet_radius, spectrum=spectrumName)
    return DMesonJetUtils.GetObject(pythiaFile, hpath)

def pythia_comparison_for_generator(gen, charm_ts, beauty_ts, jet_type, jet_radius, pythia):
    f = open("{0}.yaml".format(pythia), 'r')
    pythiaConfig = yaml.load(f)
    f.close()

    spectrumName = "JetPtDPtSpectrum"

    pythiaFileName = "{0}/{1}/{2}.root".format(pythiaConfig["input_path"], pythiaConfig["train"], pythiaConfig["name"])
    pythiaFile = ROOT.TFile(pythiaFileName)
    if not pythiaFile or pythiaFile.IsZombie():
        print("Could not open file {0}".format(pythiaFileName))
        exit(1)
    else:
        print("File {0} successfully open".format(pythiaFileName))
    pythiaSpectrum2D = GetPYTHIASpectrum(pythiaFile, jet_type, jet_radius, spectrumName)
    pythiaSpectrum2D.GetYaxis().SetRangeUser(3, 30)
    pythiaSpectrumDPt = pythiaSpectrum2D.ProjectionY("{0}_DPt_PYTHIA6".format(spectrumName), pythiaSpectrum2D.GetXaxis().FindBin(5), pythiaSpectrum2D.GetXaxis().FindBin(30 - 1e-3))
    pythiaSpectrumDPt.SetTitle("PYTHIA6 (LHC15i2bcde)")
    pythiaSpectrum2D.GetXaxis().SetRangeUser(5, 30)
    pythiaSpectrumJetPt = pythiaSpectrum2D.ProjectionX("{0}_JetPt_PYTHIA6".format(spectrumName), pythiaSpectrum2D.GetYaxis().FindBin(3), pythiaSpectrum2D.GetYaxis().FindBin(30 - 1e-3))
    pythiaSpectrumJetPt.SetTitle("PYTHIA6 (LHC15i2bcde)")
    # a factor of 1/4 comes from the fact that this is the sum of 4 cross sections
    # (the periods; in principle it should be 8 because there is a MB cross section as well
    # but at high pt the MB has vanishing statistics); a factor 1/2 comes from the antiparticles
    # so there is a missing factor 1/2
    pythiaSpectrumDPt.Scale(0.25, "width")
    pythiaSpectrumJetPt.Scale(0.25, "width")

    input_path = pythiaConfig["input_path"]

    quarks = dict()

    if charm_ts:
        charmQuark = QuarkSetting("charm")
        (charmQuark.ts, charmQuark.stage) = get_ts_stage(charm_ts)
        quarks["charm"] = charmQuark

    if beauty_ts:
        beautyQuark = QuarkSetting("beauty")
        (beautyQuark.ts, beautyQuark.stage) = get_ts_stage(beauty_ts)
        quarks["beauty"] = beautyQuark


    cnameDPt = "_".join(["PYTHIA_POWHEG_DPtSpectrumComparison", jet_type, jet_radius, pythia])
    cnameJetPt = "_".join(["PYTHIA_POWHEG_JetPtSpectrumComparison", jet_type, jet_radius, pythia])
    fname = "_".join(["PYTHIA_POWHEG_Comparison", jet_type, jet_radius, pythia])
    for quark in quarks.itervalues():
        cnameDPt += "_{0}_{1}".format(quark.name, quark.ts)
        cnameJetPt += "_{0}_{1}".format(quark.name, quark.ts)
        fname += "_{0}_{1}".format(quark.name, quark.ts)
        if quark.stage >= 0:
            quark.path = "{0}/FastSim_{1}_{2}_{3}/stage_{4}/output".format(input_path, gen, quark.name, quark.ts, quark.stage)
            quark.filename = "{0}/FastSimAnalysis_{1}_{2}_{3}.root".format(quark.path, gen, quark.name, quark.ts)
        else:
            quark.path = "{0}/FastSim_{1}_{2}_{3}/output".format(input_path, gen, quark.name, quark.ts)
            quark.filename = "{0}/FastSimAnalysis_{1}_{2}.root".format(quark.path, gen, quark.name)
        quark.file = ROOT.TFile(quark.filename)
        if not quark.file or quark.file.IsZombie():
            print("Could not open file {0}".format(quark.filename))
            exit(1)

    if gen == "powheg":
        genTitle = "POWHEG+PYTHIA6"
    else:
        genTitle = "MC"
    powhegSpectrum2D = GetTotalMCSpectrum(quarks, jet_type, jet_radius, spectrumName, genTitle)
    powhegSpectrum2D.GetYaxis().SetRangeUser(3, 30)
    powhegSpectrumDPt = powhegSpectrum2D.ProjectionY("{0}_DPt_POWHEG_PYTHIA6".format(spectrumName), powhegSpectrum2D.GetXaxis().FindBin(5), powhegSpectrum2D.GetXaxis().FindBin(30 - 1e-3))
    powhegSpectrum2D.GetXaxis().SetRangeUser(5, 30)
    powhegSpectrumJetPt = powhegSpectrum2D.ProjectionX("{0}_JetPt_POWHEG_PYTHIA6".format(spectrumName), powhegSpectrum2D.GetYaxis().FindBin(3), powhegSpectrum2D.GetYaxis().FindBin(30 - 1e-3))

    ratioAxis = "PYTHIA6 / POWHEG+PYTHIA6"

    powhegSpectrumDPt.Scale(1., "width")
    powhegSpectrumJetPt.Scale(1., "width")

    globalList.append(powhegSpectrumDPt)
    globalList.append(pythiaSpectrumDPt)
    globalList.append(powhegSpectrumJetPt)
    globalList.append(pythiaSpectrumJetPt)

    comp = DMesonJetCompare.DMesonJetCompare(cnameJetPt)
    comp.fYaxisRatio = ratioAxis
    rJetPt = comp.CompareSpectra(powhegSpectrumJetPt, [pythiaSpectrumJetPt])

    comp = DMesonJetCompare.DMesonJetCompare(cnameDPt)
    comp.fYaxisRatio = ratioAxis
    rDPt = comp.CompareSpectra(powhegSpectrumDPt, [pythiaSpectrumDPt])

    file_out = ROOT.TFile("{0}/{1}.root".format(input_path, fname), "recreate")
    file_out.cd()
    powhegSpectrumDPt.Write()
    pythiaSpectrumDPt.Write()
    powhegSpectrumJetPt.Write()
    pythiaSpectrumJetPt.Write()
    for obj in rJetPt + rDPt:
        globalList.append(obj)
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{0}/{1}.pdf".format(pythiaConfig["input_path"], obj.GetName()))
        if isinstance(obj, ROOT.TH1):
            file_out.cd()
            obj.Write()

def get_ts_stage(ts_stage):
    b = ts_stage.split(":")
    ts = int(b[0])
    if len(b) > 1:
        stage = int(b[1])
    else:
        stage = -1
    return ts, stage

def feed_down_analysis_for_generator(gen, charm_ts, beauty_ts, jet_type, jet_radius):
    rootPath = "/Volumes/DATA/ALICE/JetResults"
    charmQuark = QuarkSetting("charm")
    beautyQuark = QuarkSetting("beauty")
    quarks = dict()
    quarks["charm"] = charmQuark
    quarks["beauty"] = beautyQuark

    (charmQuark.ts, charmQuark.stage) = get_ts_stage(charm_ts)
    (beautyQuark.ts, beautyQuark.stage) = get_ts_stage(beauty_ts)

    for quark in quarks.itervalues():
        if quark.stage >= 0:
            quark.path = "{0}/FastSim_{1}_{2}_{3}/stage_{4}/output".format(rootPath, gen, quark.name, quark.ts, quark.stage)
            quark.filename = "{0}/FastSimAnalysis_{1}_{2}_{3}.root".format(quark.path, gen, quark.name, quark.ts)
        else:
            quark.path = "{0}/FastSim_{1}_{2}_{3}/output".format(rootPath, gen, quark.name, quark.ts)
            quark.filename = "{0}/FastSimAnalysis_{1}_{2}.root".format(quark.path, gen, quark.name)
        quark.file = ROOT.TFile(quark.filename)
        if not quark.file or quark.file.IsZombie():
            print("Could not open file {0}".format(quark.filename))
            exit(1)

    ptD = SpectraSet("BFeedDownVsPtD_{0}_{1}_{2}".format(gen, jet_type, jet_radius), "B feed-down vs #it{p}_{T,D}")
    ptD.add(SpectrumDef("DPtSpectrum_JetPt_0", "#it{p}_{T,ch jet} > 0", dict(colors=[ROOT.kBlue + 2, ROOT.kGreen + 2], markers=[ROOT.kFullCircle, ROOT.kFullCircle], lines=[None, None])))
    ptD.add(SpectrumDef("DPtSpectrum_JetPt_2", "#it{p}_{T,ch jet} > 2 GeV/#it{c}", dict(colors=[ROOT.kRed + 2, ROOT.kOrange + 2], markers=[ROOT.kOpenSquare, ROOT.kOpenSquare], lines=[None, None])))
    ptD.add(SpectrumDef("DPtSpectrum_JetPt_5", "#it{p}_{T,ch jet} > 5 GeV/#it{c}", dict(colors=[ROOT.kAzure + 2, ROOT.kCyan + 2], markers=[ROOT.kOpenDiamond, ROOT.kOpenDiamond], lines=[None, None])))
    ptD.add(SpectrumDef("DPtSpectrum_JetPt_8", "#it{p}_{T,ch jet} > 8 GeV/#it{c}", dict(colors=[ROOT.kMagenta + 2, ROOT.kPink + 2], markers=[ROOT.kOpenStar, ROOT.kOpenStar], lines=[None, None])))
    ptD.add(SpectrumDef("DPtSpectrum_JetPt_15", "#it{p}_{T,ch jet} > 15 GeV/#it{c}", dict(colors=[ROOT.kTeal + 2, ROOT.kSpring + 2], markers=[ROOT.kOpenCross, ROOT.kOpenCross], lines=[None, None])))

    ptJet = SpectraSet("BFeedDownVsPtJet_{0}_{1}_{2}".format(gen, jet_type, jet_radius), "B feed-down vs #it{p}_{T,ch jet}")
    ptJet.add(SpectrumDef("JetPtSpectrum_DPt_0", "#it{p}_{T,D} > 0", dict(colors=[ROOT.kBlue + 2, ROOT.kGreen + 2], markers=[ROOT.kFullCircle, ROOT.kFullCircle], lines=[None, None])))
    ptJet.add(SpectrumDef("JetPtSpectrum_DPt_2", "#it{p}_{T,D} > 2 GeV/#it{c}", dict(colors=[ROOT.kRed + 2, ROOT.kOrange + 2], markers=[ROOT.kOpenSquare, ROOT.kOpenSquare], lines=[None, None])))

    zJet = SpectraSet("BFeedDownVsZ_{0}_{1}_{2}".format(gen, jet_type, jet_radius), "B feed-down vs #it{z}_{||,D}^{ch jet}")
    zJet.add(SpectrumDef("JetZSpectrum_JetPt_0_5", "0 < #it{p}_{T,ch jet} < 5 GeV/#it{c}", dict(colors=[ROOT.kBlue + 2, ROOT.kGreen + 2], markers=[ROOT.kFullCircle, ROOT.kFullCircle], lines=[None, None])))
    zJet.add(SpectrumDef("JetZSpectrum_JetPt_5_10", "5 < #it{p}_{T,ch jet} < 10 GeV/#it{c}", dict(colors=[ROOT.kRed + 2, ROOT.kOrange + 2], markers=[ROOT.kOpenSquare, ROOT.kOpenSquare], lines=[None, None])))
    zJet.add(SpectrumDef("JetZSpectrum_JetPt_10_15", "10 < #it{p}_{T,ch jet} < 15 GeV/#it{c}", dict(colors=[ROOT.kAzure + 2, ROOT.kCyan + 2], markers=[ROOT.kOpenDiamond, ROOT.kOpenDiamond], lines=[None, None])))
    zJet.add(SpectrumDef("JetZSpectrum_JetPt_15_20", "15 < #it{p}_{T,ch jet} < 20 GeV/#it{c}", dict(colors=[ROOT.kMagenta + 2, ROOT.kPink + 2], markers=[ROOT.kOpenStar, ROOT.kOpenStar], lines=[None, None])))
    zJet.add(SpectrumDef("JetZSpectrum_JetPt_20_25", "20 < #it{p}_{T,ch jet} < 25 GeV/#it{c}", dict(colors=[ROOT.kTeal + 2, ROOT.kSpring + 2], markers=[ROOT.kStar, ROOT.kStar], lines=[None, None])))
    zJet.add(SpectrumDef("JetZSpectrum_JetPt_25", "#it{p}_{T,ch jet} > 25 GeV/#it{c}", dict(colors=[ROOT.kViolet + 2, ROOT.kYellow + 2], markers=[ROOT.kOpenCross, ROOT.kOpenCross], lines=[None, None])))

    spectraSets = [ptD, ptJet, zJet]

    ratioAxis = "({0} #rightarrow D^{{0}}) / ({1} #rightarrow D^{{0}})".format(quarks.values()[1].name[0], quarks.values()[0].name[0])

    for spectraSet in spectraSets:
        comp = DMesonJetCompare.DMesonJetCompare(spectraSet.name)
        comp.fOptSpectrum = ""
        comp.fOptRatio = ""
        comp.fYaxisRatio = ratioAxis
        comp.fX1LegSpectrum = 0.10
        comp.fX1LegRatio = 0.4
        comp.fNColsLegSpectrum = 2
        comp.fLogUpperSpace = 500
        comp.fLegTextSize = 18
        for spectrumDef in spectraSet.spectra:
            comp.fColors = spectrumDef.styles["colors"]
            comp.fMarkers = spectrumDef.styles["markers"]
            comp.fLines = spectrumDef.styles["lines"]
            histos = []
            for quark in quarks.itervalues():
                h = GetSpectrum(quark.file, "D0_MCTruth", jet_type, jet_radius, spectrumDef.name)
                h.Scale(1., "width")
                h.SetName("{0}_{1}".format(quark.name, h.GetName()))
                h.SetTitle("{0} #rightarrow D^{{0}}, {1}".format(quark.name[0], spectrumDef.title))
                histos.append(h)
                globalList.append(h)
                quark.histos[spectrumDef] = h

            r = comp.CompareSpectra(histos[0], histos[1:])
            comp.fOptSpectrum = "same"
            comp.fOptRatio = "same"
            for obj in r:
                globalList.append(obj)
                if isinstance(obj, ROOT.TLegend):
                    if "Ratio" in obj.GetName():
                        legRatio = obj
                    else:
                        leg = obj
                elif isinstance(obj, ROOT.TCanvas):
                    if "Ratio" in obj.GetName():
                        cRatio = obj
                    else:
                        c = obj
        c.SaveAs("{0}/{1}.pdf".format(rootPath, c.GetName()))
        cRatio.SaveAs("{0}/{1}.pdf".format(rootPath, cRatio.GetName()))

    corr_fact, corr_fact_err = generate_correction_factors(quarks, jet_type, jet_radius)
    globalList.append(corr_fact)
    globalList.append(corr_fact_err)

    cname = "BFeedDown_CorrFact_{0}_{1}_{2}_{3}_{4}".format(gen, jet_type, jet_radius, charmQuark.ts, beautyQuark.ts)
    c = ROOT.TCanvas(cname, cname)
    c.cd()
    globalList.append(c)
    corr_fact.SetMaximum(1)
    corr_fact.Draw("colz")
    c.SetRightMargin(0.15)
    c.SaveAs("{0}/{1}.pdf".format(rootPath, c.GetName()))

    cname = "BFeedDown_CorrFactUnc_{0}_{1}_{2}_{3}_{4}".format(gen, jet_type, jet_radius, charmQuark.ts, beautyQuark.ts)
    c = ROOT.TCanvas(cname, cname)
    c.cd()
    globalList.append(c)
    corr_fact_err.Draw("colz")
    c.SetRightMargin(0.15)
    c.SaveAs("{0}/{1}.pdf".format(rootPath, c.GetName()))

def GetSpectrum(file, meson_name, jet_type, jet_radius, spectrum):
    mesonlistname = meson_name
    mesonlist = file.Get(mesonlistname)
    if not mesonlist:
        print("Could not get list {0} from file {1}".format(mesonlistname, file.GetName()))
        file.ls()
        exit(1)
    jetlistname = "{0}_{1}".format(jet_type, jet_radius)
    jetlist = mesonlist.FindObject(jetlistname)
    if not jetlist:
        print("Could not get list {0} from list {1} in file {2}".format(jetlistname, mesonlistname, file.GetName()))
        mesonlist.Print()
        exit(1)
    spectrumlistname = "_".join([meson_name, jet_type, jet_radius, spectrum])
    spectrumlist = jetlist.FindObject(spectrumlistname)
    if not spectrumlist:
        print("Could not get list {0} from list {1} in list {2} in file {3}".format(spectrumlistname, jetlistname, mesonlistname, file.GetName()))
        jetlist.Print()
        exit(1)
    spectrumname = "_".join([meson_name, jet_type, jet_radius, spectrum])
    h = spectrumlist.FindObject(spectrumname)
    if not h:
        print("Could not find object {0} in list {1}/{2}/{3} in file {4}".format(spectrumname, mesonlistname, jetlistname, spectrumlistname, file.GetName()))
        jetlist.Print()
        exit(1)
    h_copy = h.Clone("{0}_copy".format(spectrum))
    return h_copy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='B feed-down.')
    parser.add_argument('--charm', metavar='timestamp',
                        default=None)
    parser.add_argument('--beauty', metavar='timestamp',
                        default=None)
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    parser.add_argument('--data', metavar='LHC10_Train823_LHC15i2_Train961_efficiency',
                        default=None)
    parser.add_argument('--pythia', metavar='LHC15i2analysis_Train961',
                        default=None)
    args = parser.parse_args()

    main(args.charm, args.beauty, args.jet_type, args.jet_radius, args.data, args.pythia)

    IPython.embed()
