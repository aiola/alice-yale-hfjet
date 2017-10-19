#!/usr/bin/env python
# python program to perform a D meson jet analysis

import math
import copy
import collections
import os

import ROOT

import DMesonJetProjectors
import DMesonJetCompare
import DMesonJetUtils
from DMesonJetBase import AnalysisType
import BinSet
import Axis

globalList = []

class DMesonJetTopoContainer:
    def __init__(self, jet_pt_bins, jet_type, jet_radius):
        self.fJetPtBins = jet_pt_bins
        self.fHistograms = dict()
        self.fJetType = jet_type
        self.fJetRadius = jet_radius
        for jetPtBin in self.fJetPtBins:
            self.fHistograms[jetPtBin] = self.GenerateHistograms()

    def AddHistogram(self, histograms, h):
        histograms[h.GetName()] = h

    def GenerateHistograms(self):
        histograms = dict()

        h = ROOT.TH1F("fDCA", "fDCA", 250, 0, 2.5)
        h.Sumw2()
        self.AddHistogram(histograms, h)

        h = ROOT.TH1F("fCosThetaStar", "fCosThetaStar", 200, -1, 1)
        h.Sumw2()
        self.AddHistogram(histograms, h)

        h = ROOT.TH1F("fd0d0", "fd0d0", 1000, -5e-4, 5e-4)
        h.Sumw2()
        self.AddHistogram(histograms, h)

        h = ROOT.TH1F("fMaxNormd0", "fMaxNormd0", 2000, -100, 100)
        h.Sumw2()
        self.AddHistogram(histograms, h)

        h = ROOT.TH1F("fCosPointing", "fCosPointing", 200, -1, 1)
        h.Sumw2()
        self.AddHistogram(histograms, h)

        return histograms

    def Fill(self, event, eventWeight):
        dmeson = event.DmesonJet

        if self.fJetType or self.fJetRadius:
            jetName = "Jet_AKT{0}{1}_pt_scheme".format(self.fJetType, self.fJetRadius)
            jet = getattr(event, jetName)
        else:
            jet = None

        for (jetPtMin, jetPtMax) in self.fJetPtBins:
            if jet.fPt < jetPtMin or jet.fPt >= jetPtMax: continue
            for h in self.fHistograms[(jetPtMin, jetPtMax)].itervalues():
                v = getattr(dmeson, h.GetName())
                h.Fill(v, eventWeight)

class DMesonJetTopoAnalysisManager:
    def __init__(self, trigger, dmeson, projector):
        self.fBkgAnalysis = DMesonJetTopoAnalysis("Background", trigger, "{}_kBackgroundOnly_D0toKpiCuts_loosest_nopid".format(dmeson), projector)
        self.fSigAnalysis = DMesonJetTopoAnalysis("Signal", trigger, "{}_kSignalOnly_D0toKpiCuts_loosest_nopid".format(dmeson), projector)

    def StartAnalysis(self):
        self.fBkgAnalysis.Start()
        self.fSigAnalysis.Start()
        self.Plot()

    def SaveRootFile(self, path):
        pass

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            if c: c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))

    def Plot(self):
        self.fCanvases = []
        self.fCompareObjects = []
        for (jetPtMin, jetPtMax), sigEfficiencies, bkgEfficiencies in zip(self.fSigAnalysis.fHistograms.iterkeys(), self.fSigAnalysis.fCutEfficiency.itervalues(), self.fBkgAnalysis.fCutEfficiency.itervalues()):
            for sigEfficiency, bkgEfficiency in zip(sigEfficiencies.itervalues(), bkgEfficiencies.itervalues()):
                cname = "{}_{}_{}".format(sigEfficiency.GetName(), jetPtMin, jetPtMax)
                comp = DMesonJetCompare.DMesonJetCompare(cname)
                comp.fDoSpectraPlot = "lineary"
                comp.fDoRatioPlot = False
                sigEfficiency.SetTitle("Signal")
                bkgEfficiency.SetTitle("Background")
                comp.CompareSpectra(sigEfficiency, [bkgEfficiency])
                self.fCompareObjects.append(comp)
                self.fCanvases.append(comp.fCanvasSpectra)

class DMesonJetTopoAnalysis:
    def __init__(self, name, trigger, dmeson, projector):
        self.fName = name
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fProjector = projector

    def SaveRootFile(self, file):
        rlist = ROOT.TList()
        if self.fTrigger:
            rlist.SetName("{}_{}".format(self.fTrigger, self.fDMeson))
        else:
            rlist.SetName(self.fDMeson)

        for (jetPtMin, jetPtMax), histograms in self.fHistograms.iteritems():
            hListName = "Histograms_{}_{}".format(jetPtMin, jetPtMax)
            hList = ROOT.TList()
            hList.SetName(jetName)
            for h in histograms:
                hList.Add(h)
            rlist.Add(hList)

        for (jetPtMin, jetPtMax), histograms in self.fCutEfficiency.iteritems():
            hListName = "CutEfficiency_{}_{}".format(jetPtMin, jetPtMax)
            hList = ROOT.TList()
            hList.SetName(jetName)
            for h in histograms:
                hList.Add(h)
            rlist.Add(hList)

        if rlist.GetEntries() > 0:
            file.cd()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)

    def DoProjections(self):
        data_container = DMesonJetTopoContainer([(5, 15), (15, 30)], "Charged", "R040")
        self.fProjector.StartProjection(self.fTrigger, self.fDMeson, None, data_container)
        self.fHistograms = data_container.fHistograms

    def Start(self):
        self.DoProjections()
        self.Analyze()

    def Analyze(self):
        self.fCutEfficiency = dict()
        for (jetPtMin, jetPtMax), histograms in self.fHistograms.iteritems():
            self.fCutEfficiency[(jetPtMin, jetPtMax)] = dict()
            for h in histograms.itervalues():
                integral = h.Integral(0, -1)
                if integral == 0: continue
                h.Scale(1. / integral)
                cutEff = h.Clone("{}_Efficiency".format(h.GetName()))
                cutEff.Reset()
                self.fCutEfficiency[(jetPtMin, jetPtMax)][cutEff.GetName()] = cutEff
                partialIntErr = 0.
                partialInt = 0.
                for ibin in xrange(0, h.GetNbinsX() + 2):
                    partialIntErr = math.sqrt(partialIntErr ** 2 + h.GetBinError(ibin) ** 2)
                    partialInt += h.GetBinContent(ibin)
                    cutEff.SetBinContent(ibin, partialInt)
                    cutEff.SetBinError(ibin, partialIntErr)
