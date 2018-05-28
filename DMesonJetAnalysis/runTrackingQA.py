#!/usr/bin/env python
# python script to produce tracking QA plots
# use it with:
# train Jets_EMC_pp_MC_1522_1508_1524_1525 (MC)
# train Jets_EMC_pp_1151_1152_1153_1154 (data)

import argparse
import IPython
import ROOT
import os
import DMesonJetCompare
import DMesonJetUtils
import numpy
import yaml
from enum import Enum

globalList = []


def ResetTHnFilters(hn):
    for iaxis in range(0, hn.GetNdimensions()):
        hn.GetAxis(iaxis).SetRange(0, -1)

class Variable2D:
    def __init__(self, varx, vary, plotsigma):
        self.fVariableX = varx
        self.fVariableY = vary
        self.fPlotSigma = plotsigma
    
    def GetName(self):
        return "{}vs{}".format(self.fVariableY.fName, self.fVariableX.fName)

class Variable:

    def __init__(self, name, title, units, logy, tracks_axis, particles_axis, matched_tracks_axis, matched_particles_axis):
        self.fName = name
        self.fTitle = title
        self.fUnits = units
        self.fTracksAxis = tracks_axis
        self.fParticlesAxis = particles_axis
        self.fMatchedTracksAxis = matched_tracks_axis
        self.fMatchedParticlesAxis = matched_particles_axis
        self.fLogy = logy
        self.fLogx = False
        self.fBins = None
    
    def GetName(self):
        return self.fName

    def GetAxisTitle(self):
        axis_title = self.fTitle
        if self.fUnits: axis_title += " " + self.fUnits
        return axis_title

class HistogramType(Enum):
    Track = 1
    Particle = 2
    MatchedTrack = 3
    MatchedParticle = 4

class Histogram:

    def __init__(self, variable, htype):
        self.fVariable = variable
        self.fHistogramType = htype
        self.fHistograms = dict()
        self.fNormHistograms = dict()
        self.fTrackTypes = { 1 : "Global tracks", 2 : "Constrained tracks w/ ITS refit", 3 : "Constrained tracks w/o ITS refit"}

    def MergeHistograms(self):
        if len(self.fHistograms) < 2: return
        htot = None
        for h in self.fHistograms.itervalues():
            if htot:
                htot.Add(h)
            else:
                htot = h.Clone(h.GetName()[:-2])
        htot.SetTitle("Hybrid tracks")
        self.fHistograms["all"] = htot

    def NormalizeHistogram(self, events):
        for tt, h in self.fHistograms.iteritems():
            hnorm = h.Clone("{}_Normalized".format(h.GetName()))
            hnorm.Scale(1.0 / events, "width")
            if self.fVariable.fUnits:
                yaxis_title = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d#it{{N}}}}{{d{var}}} {units}^{{-1}}".format(var=self.fVariable.fTitle, units=self.fVariable.fUnits)
            else:
                yaxis_title = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d#it{{N}}}}{{d{var}}}".format(var=self.fVariable.fTitle)
            hnorm.GetYaxis().SetTitle(yaxis_title)
            self.fNormHistograms[tt] = hnorm

    def GetFullHistogram(self, normalized=True):
        if normalized:
            return self.fNormHistograms["all"]
        else:
            return self.fHistograms["all"]

    def GetPartialHistograms(self, normalized=True):
        if normalized:
            l = self.fNormHistograms
        else:
            l = self.fHistograms
        return [obj for k, obj in l.iteritems() if not k == "all"]

class Histogram2D(Histogram):
    def RebinIfNeeded(self, h):
        if not self.fVariable.fVariableX.fBins and not self.fVariable.fVariableY.fBins: 
            return h
        else:
            if self.fVariable.fVariableX.fBins:
                nxbins = len(self.fVariable.fVariableX.fBins) - 1
                xbins = numpy.array(self.fVariable.fVariableX.fBins, dtype=numpy.float64)
            else:
                nxbins = h.GetXaxis().GetNbins()
                xbins = h.GetXaxis().GetXbins().GetArray()
            if self.fVariable.fVariableY.fBins:
                nybins = len(self.fVariable.fVariableY.fBins) - 1
                ybins = numpy.array(self.fVariable.fVariableY.fBins, dtype=numpy.float64)
            else:
                nybins = h.GetYaxis().GetNbins()
                ybins = h.GetYaxis().GetXbins().GetArray()

            hrebinned = DMesonJetUtils.Rebin2D_fromBins(h, h.GetName(), nxbins, xbins, nybins, ybins, False)
        return hrebinned

    def GenerateProfiles(self):
        self.fProfiles = dict()
        for k, h in self.fHistograms.iteritems():
            if self.fVariable.fPlotSigma:
                p = h.ProfileX("{}_Profile".format(h.GetName()), 1, -1, "s")
                hname = "{}_Sigma".format(h.GetName())
                if h.GetXaxis().GetXbins():
                    s = ROOT.TH1D(hname, h.GetTitle(), h.GetXaxis().GetNbins(), h.GetXaxis().GetXbins().GetArray())
                else:
                    s = ROOT.TH1D(hname, h.GetTitle(), h.GetXaxis().GetNbins(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
                for ibin in range(0, p.GetNbinsX() + 2):
                    s.SetBinContent(ibin, p.GetBinError(ibin))
                    s.SetBinError(ibin, 0)
                self.fProfiles[k] = s
            else:
                p = h.ProfileX("{}_Profile".format(h.GetName()), 1, -1, "i")
                self.fProfiles[k] = p

    def GetFullProfile(self):
        return self.fProfiles["all"]

    def GetPartialProfiles(self, normalized=True):
        return [obj for k, obj in self.fProfiles.iteritems() if not k == "all"]

class Histogram1D(Histogram):
    def RebinIfNeeded(self, h):
        if not self.fVariable.fBins: return h
        hrebinned = DMesonJetUtils.Rebin1D_fromBins(h, h.GetName(), len(self.fVariable.fBins) - 1, numpy.array(self.fVariable.fBins, dtype=numpy.float64), False)
        return hrebinned

class TrackHistogram(Histogram1D):

    def Add(self, hn, track_type_axis):
        print("Adding histogram '{}' (variable '{}')...".format(hn.GetName(), self.fVariable.fName))
        for track_type, track_type_title in self.fTrackTypes.iteritems():
            hn.GetAxis(track_type_axis).SetRange(track_type, track_type)
            if self.fHistogramType == HistogramType.MatchedTrack:
                h = self.RebinIfNeeded(hn.Projection(self.fVariable.fMatchedTracksAxis, "A"))
                hname = "MatchedTracks{}_{}".format(self.fVariable.fName, track_type)
            else:
                h = self.RebinIfNeeded(hn.Projection(self.fVariable.fTracksAxis, "A"))
                hname = "Tracks{}_{}".format(self.fVariable.fName, track_type)
            hn.GetAxis(track_type_axis).SetRange(0, -1)
            xaxis_title = self.fVariable.GetAxisTitle()
            h.SetName(hname)
            h.SetTitle(track_type_title)
            h.GetXaxis().SetTitle(xaxis_title)
            h.GetYaxis().SetTitle("counts")
            h.Sumw2()
            if track_type in self.fHistograms:
                self.fHistograms[track_type].Add(h)
            else:
                self.fHistograms[track_type] = h


class ParticleHistogram(Histogram1D):

    def Add(self, hn, track_type_axis):
        if self.fHistogramType == HistogramType.MatchedParticle:
            track_types = self.fTrackTypes
        else:
            track_types = { "all" : "Hybrid tracks" }
        for track_type, track_type_title in track_types.iteritems():
            if self.fHistogramType == HistogramType.MatchedParticle:
                hn.GetAxis(track_type_axis).SetRange(track_type, track_type)
                h = self.RebinIfNeeded(hn.Projection(self.fVariable.fMatchedParticlesAxis, "A"))
                hn.GetAxis(track_type_axis).SetRange(0, -1)
                hname = "MatchedParticles{}_{}".format(self.fVariable.fName, track_type)
            else:
                hn.GetAxis(track_type_axis).SetRange(2, 2)  # findable
                h = self.RebinIfNeeded(hn.Projection(self.fVariable.fParticlesAxis, "A"))
                hn.GetAxis(track_type_axis).SetRange(0, -1)
                hname = "Particles{}".format(self.fVariable.fName)
            xaxis_title = self.fVariable.GetAxisTitle()
            h.SetName(hname)
            h.SetTitle(track_type_title)
            h.GetXaxis().SetTitle(xaxis_title)
            h.GetYaxis().SetTitle("counts")
            h.Sumw2()
            if track_type in self.fHistograms:
                self.fHistograms[track_type].Add(h)
            else:
                self.fHistograms[track_type] = h

class TrackHistogram2D(Histogram2D):

    def Add(self, hn, track_type_axis):
        print("Adding histogram '{}' (variable '{}':'{}')...".format(hn.GetName(), self.fVariable.fVariableX.fName, self.fVariable.fVariableY.fName))
        for track_type, track_type_title in self.fTrackTypes.iteritems():
            hn.GetAxis(track_type_axis).SetRange(track_type, track_type)
            if self.fHistogramType == HistogramType.MatchedTrack:
                h = self.RebinIfNeeded(hn.Projection(self.fVariable.fVariableY.fMatchedTracksAxis, self.fVariable.fVariableX.fMatchedTracksAxis, "A"))
                hname = "MatchedTracks{}_{}".format(self.fVariable.GetName(), track_type)
            else:
                h = self.RebinIfNeeded(hn.Projection(self.fVariable.fVariableY.fTracksAxis, self.fVariable.fVariableX.fTracksAxis, "A"))
                hname = "Tracks{}_{}".format(self.fVariable.GetName(), track_type)
            hn.GetAxis(track_type_axis).SetRange(0, -1)
            xaxis_title = self.fVariable.fVariableX.GetAxisTitle()
            yaxis_title = self.fVariable.fVariableY.GetAxisTitle()
            h.SetName(hname)
            h.SetTitle(track_type_title)
            h.GetXaxis().SetTitle(xaxis_title)
            h.GetYaxis().SetTitle(yaxis_title)
            h.GetZaxis().SetTitle("counts")
            h.Sumw2()
            if track_type in self.fHistograms:
                self.fHistograms[track_type].Add(h)
            else:
                self.fHistograms[track_type] = h

class EfficiencyHistogram(Histogram):

    def Generate(self, all_histograms, matched_histograms):
        hall = all_histograms.fHistograms["all"]
        for track_type, hmatched in matched_histograms.fHistograms.iteritems():
            eff = hmatched.Clone("{}_Efficiency".format(hmatched.GetName()))
            eff.Divide(hall)
            eff.GetYaxis().SetTitle("Efficiency")
            self.fHistograms[track_type] = eff
        self.fNormHistograms = self.fHistograms


class TrackingQA:

    def __init__(self, input_path, train):
        if "_MC_" in train:
            print("Train '{}' is a MC train.".format(train))
            self.fIsMC = True
            self.fHnParticlesName = "fParticlesPhysPrim"
            self.fHnMatchedName = "fParticlesMatched"
        else:
            print("Train '{}' is NOT a MC train.".format(train))
            self.fIsMC = False
            self.fHnParticlesName = ""
            self.fHnMatchedName = ""
        self.fInputPath = input_path
        self.fTrain = train
        self.fTaskName = "AliEmcalTrackingQATask_histos"
        self.fHnTracksName = "fTracks"
        self.fVariables = []
        pt = Variable("Pt", "#it{p}_{T}", "(GeV/#it{c})", True, 0, 0, 3, 0)
        pt.fBins = [0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
                    1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,
                    3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
                    15, 20, 30, 40]
        pt.fLogx = True
        eta = Variable("Eta", "#eta", "", True, 1, 1, 4, 1)
        phi = Variable("Phi", "#phi", "", True, 2, 2, 5, 2)
        sigmapt = Variable("SigmaPt", "#sigma(#it{p}_{T}) / #it{p}_{T}", "", True, 5, -1, -1, -1)
        reldiffpt = Variable("RelDiffPt", "(#it{p}_{T}^{gen} - #it{p}_{T}^{det}) / #it{p}_{T}^{det}", "", True, -1, -1, 7, -1)
        self.fVariables.append(pt)
        self.fVariables.append(eta)
        self.fVariables.append(phi)
        self.fVariables.append(Variable2D(pt, sigmapt, False))
        self.fVariables.append(Variable2D(eta, sigmapt, False))
        self.fVariables.append(Variable2D(phi, sigmapt, False))
        self.fVariables.append(Variable2D(pt, reldiffpt, True))
        self.fVariables.append(Variable2D(eta, reldiffpt, True))
        self.fVariables.append(Variable2D(phi, reldiffpt, True))
        self.fTrackTypeAxis = 4
        self.fParticleFindableAxis = 4
        self.fMatchedTrackTypeAxis = 6
        self.fTrackHistograms = dict()
        self.fParticleHistograms = dict()
        self.fMatchedTracksHistograms = dict()
        self.fMatchedParticleHistograms = dict()
        self.fEfficiencyHistograms = dict()
        self.fFakeRateHistograms = dict()
        self.fHnTracks = []
        self.fHnParticles = []
        self.fHnMatched = []
        self.fNEvents = []
        self.fTotEvents = 0
        self.fCanvases = []

    def ExtractPeriod(self, fname):
        train_index = fname.find(self.fTrain)
        index1 = fname.find("/", train_index)
        index2 = fname.find("/", index1 + 1)
        period = fname[index1 + 1: index2]
        return period

    def LoadHns(self):
        file_name = "AnalysisResults.root"
        path = "{0}/{1}".format(self.fInputPath, self.fTrain)
        print("Looking for file {0} in path {1}".format(file_name, path))
        file_names = DMesonJetUtils.find_file(path, file_name)

        for fname in file_names:
            print("Opening file '{}'".format(fname))
            period = self.ExtractPeriod(fname)

            file = ROOT.TFile(fname, "read")

            h = DMesonJetUtils.GetObject(file, "{}/{}".format(self.fTaskName, self.fHnTracksName))
            h.SetName("{}_{}".format(h.GetName(), period))
            self.fHnTracks.append(h)

            if self.fIsMC:
                h = DMesonJetUtils.GetObject(file, "{}/{}".format(self.fTaskName, self.fHnParticlesName))
                h.SetName("{}_{}".format(h.GetName(), period))
                self.fHnParticles.append(h)

                h = DMesonJetUtils.GetObject(file, "{}/{}".format(self.fTaskName, self.fHnMatchedName))
                h.SetName("{}_{}".format(h.GetName(), period))
                self.fHnMatched.append(h)

            h = DMesonJetUtils.GetObject(file, "{}/{}".format(self.fTaskName, "fHistEventCount"))
            events = h.GetBinContent(1)
            self.fNEvents.append(events)
            self.fTotEvents += events

    def ProjetHistogram1D(self, var):
        track_histogram = TrackHistogram(var, HistogramType.Track)
        for hn in self.fHnTracks:
            track_histogram.Add(hn, self.fTrackTypeAxis)
        track_histogram.MergeHistograms()
        track_histogram.NormalizeHistogram(self.fTotEvents)
        self.fTrackHistograms[var.fName] = track_histogram

        if self.fIsMC:
            particle_histogram = ParticleHistogram(var, HistogramType.Particle)
            for hn in self.fHnParticles:
                particle_histogram.Add(hn, self.fParticleFindableAxis)
            particle_histogram.MergeHistograms()
            particle_histogram.NormalizeHistogram(self.fTotEvents)
            self.fParticleHistograms[var.fName] = particle_histogram

            matched_track_histogram = TrackHistogram(var, HistogramType.MatchedTrack)
            matched_particle_histogram = ParticleHistogram(var, HistogramType.MatchedParticle)
            for hn in self.fHnMatched:
                matched_track_histogram.Add(hn, self.fMatchedTrackTypeAxis)
                matched_particle_histogram.Add(hn, self.fMatchedTrackTypeAxis)
            matched_track_histogram.MergeHistograms()
            matched_particle_histogram.MergeHistograms()
            matched_track_histogram.NormalizeHistogram(self.fTotEvents)
            matched_particle_histogram.NormalizeHistogram(self.fTotEvents)
            self.fMatchedTracksHistograms[var.fName] = matched_track_histogram
            self.fMatchedParticleHistograms[var.fName] = matched_particle_histogram

    def DoProjetHistogram2D(self, histogram, hns, track_type_axis):
        for hn in hns: histogram.Add(hn, track_type_axis)
        histogram.MergeHistograms()
        histogram.GenerateProfiles()
        return histogram

    def ProjetHistogram2D(self, var):
        print("Working on variable '{}'".format(var.GetName()))
        if var.fVariableY.fTracksAxis >= 0:
            print("Histogram type is track.")
            track_histogram = TrackHistogram2D(var, HistogramType.Track)
            self.fTrackHistograms[var.GetName()] = self.DoProjetHistogram2D(track_histogram, self.fHnTracks, self.fTrackTypeAxis)
        
        if self.fIsMC:
            if var.fVariableY.fMatchedTracksAxis >= 0:
                print("Histogram type is matched track.")
                matched_track_histogram = TrackHistogram2D(var, HistogramType.MatchedTrack)
                self.fMatchedTracksHistograms[var.GetName()] = self.DoProjetHistogram2D(matched_track_histogram, self.fHnMatched, self.fMatchedTrackTypeAxis)

    def ProjectHistograms(self):
        for var in self.fVariables:
            if isinstance(var, Variable): 
                self.ProjetHistogram1D(var)
            elif isinstance(var, Variable2D): 
                self.ProjetHistogram2D(var)

    def Plot(self):
        self.PlotSpectra("Tracks", self.fTrackHistograms)
        self.Plot2DHistograms("Tracks", self.fTrackHistograms)
        self.Plot2DHistograms("MatchedTracks", self.fMatchedTracksHistograms)
        self.PlotProfiles()
        if self.fIsMC:
            self.PlotSpectra("MatchedTracks", self.fMatchedTracksHistograms)
            self.PlotSpectra("MatchedParticles", self.fMatchedParticleHistograms)
            self.PlotEfficiency()

    def Plot2DHistograms(self, name, histograms):
        for var in self.fVariables:
            if not isinstance(var, Variable2D): continue
            if not var.GetName() in histograms: continue
            cname = "{}_{}".format(name, var.GetName())
            h = histograms[var.GetName()].GetFullHistogram(False)
            c = ROOT.TCanvas(cname, cname)
            h.Draw("colz")
            globalList.append(c)

    def PlotSpectra(self, name, histograms):
        for var in self.fVariables:
            if not isinstance(var, Variable): continue
            if not var.fName in histograms: continue
            comp_name = "{}_{}".format(name, var.fName)
            comp = DMesonJetCompare.DMesonJetCompare(comp_name)
            comp.fDoRatioPlot = "logy"
            comp.fX1LegRatio = 0.15
            comp.fX2LegRatio = 0.50
            if var.fLogy:
                comp.fDoSpectraPlot = "logy"
            else:
                comp.fDoSpectraPlot = "lineary"
            baseline = histograms[var.fName].GetFullHistogram()
            histos = histograms[var.fName].GetPartialHistograms()
            tot = baseline.Integral()
            for h in histos:
                partial = h.Integral() / tot
                print("Fraction of '{}': {}".format(h.GetTitle(), partial))
            results = comp.CompareSpectra(baseline, histos)
            for r in results:
                if isinstance(r, ROOT.TCanvas):
                    if var.fLogx: r.SetLogx()
                    self.fCanvases.append(r)
                globalList.append(r)

    def PlotEfficiency(self):
        for var in self.fVariables:
            if not isinstance(var, Variable): continue
            comp_name = "Efficiency_{}".format(var.fName)
            comp = DMesonJetCompare.DMesonJetCompare(comp_name)
            comp.fDoRatioPlot = False
            comp.fDoSpectraPlot = "lineary"
            comp.fX1LegSpectrum = 0.15
            comp.fX2LegSpectrum = 0.50
            comp.fLinUpperSpace = 0.4
            comp.fOptSpectrum = "hist"
            comp.fOptSpectrumBaseline = "hist"
            comp.fGridySpectrum = True
            baseline = self.fEfficiencyHistograms[var.GetName()].GetFullHistogram()
            histos = self.fEfficiencyHistograms[var.GetName()].GetPartialHistograms()
            results = comp.CompareSpectra(baseline, histos)
            for r in results:
                if isinstance(r, ROOT.TCanvas):
                    if var.fLogx: r.SetLogx()
                    self.fCanvases.append(r)
                globalList.append(r)

    def PlotProfiles(self):
        for var in self.fVariables:
            if not isinstance(var, Variable2D): continue
            if var.fVariableY.fTracksAxis >= 0:
                histograms = self.fTrackHistograms
            elif var.fVariableY.fMatchedTracksAxis >= 0:
                histograms = self.fMatchedTracksHistograms
            else:
                continue
            if not var.GetName() in histograms: continue
            comp_name = "Profile_{}".format(var.GetName())
            comp = DMesonJetCompare.DMesonJetCompare(comp_name)
            comp.fDoRatioPlot = False
            comp.fDoSpectraPlot = "lineary"
            comp.fX1LegSpectrum = 0.15
            comp.fX2LegSpectrum = 0.50
            comp.fLinUpperSpace = 0.4
            comp.fOptSpectrum = "hist"
            comp.fOptSpectrumBaseline = "hist"
            comp.fGridySpectrum = True

            baseline = histograms[var.GetName()].GetFullProfile()
            histos = histograms[var.GetName()].GetPartialProfiles()
            results = comp.CompareSpectra(baseline, histos)
            if comp.fMainHistogram.GetMaximum() > 0.3:
                comp.fMainHistogram.GetYaxis().SetRangeUser(0, 0.3)
            for r in results:
                if isinstance(r, ROOT.TCanvas):
                    self.fCanvases.append(r)
                globalList.append(r)

    def CalculateEfficiency(self):
        if not self.fIsMC: return
        for var in self.fVariables:
            if not isinstance(var, Variable): continue
            efficiencies = EfficiencyHistogram(var, True)
            efficiencies.Generate(self.fParticleHistograms[var.fName], self.fMatchedParticleHistograms[var.fName])
            self.fEfficiencyHistograms[var.fName] = efficiencies

    def SaveCanvases(self):
        for c in self.fCanvases:
            fname = "{}/{}/{}.pdf".format(self.fInputPath, self.fTrain, c.GetName())
            c.SaveAs(fname)


def main(input_path, train):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    qa = TrackingQA(input_path, train)
    globalList.append(qa)
    qa.LoadHns()
    qa.ProjectHistograms()
    qa.CalculateEfficiency()
    qa.Plot()
    qa.SaveCanvases()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tracking QA.')
    parser.add_argument('train', metavar='train')
    parser.add_argument('--input-path', metavar='input-path',
                        default="/Volumes/DATA/ALICE/JetResults")
    args = parser.parse_args()

    main(args.input_path, args.train)

    IPython.embed()
