#!/usr/bin/env python
# python program to perform a D meson jet analysis

import math
import collections
import os
import numpy

import ROOT

import DMesonJetProjectors
import DMesonJetCompare

globalList = []


class DMesonJetPIDContainer:

    def __init__(self, mc, jet_pt_bins, d_pt_bins, jet_type, jet_radius):
        self.fJetPtBins = jet_pt_bins
        self.fDPtBins = d_pt_bins
        self.fJetType = jet_type
        self.fJetRadius = jet_radius
        self.fMC = mc

        self.fMinDPt = 2
        self.fMaxDPt = 30

        self.fMinJetPt = 5
        self.fMaxJetPt = 30

        self.fMaxJetEta = 0.5

        self.fNoPID = dict()
        self.fD0 = dict()
        self.fD0bar = dict()
        self.fMCSignal = dict()
        self.fMCWrongPID = dict()

        if self.fJetType and self.fJetRadius:
            self.fJetName = "Jet_AKT{0}{1}_pt_scheme".format(self.fJetType, self.fJetRadius)
        else:
            self.fJetName = None

        self.GenerateAllHistograms()

    def LoopAllHistos(self):
        for h in self.fNoPID.itervalues(): yield h
        for h in self.fD0.itervalues(): yield h
        for h in self.fD0bar.itervalues(): yield h
        for h in self.fMCSignal.itervalues(): yield h
        for h in self.fMCWrongPID.itervalues(): yield h

    def DoSumw2(self):
        for h in self.LoopAllHistos():
            h.Scale(1.0, "width")
            h.Sumw2()

    def SetMC(self, switch):
        if self.fMC and switch:
            if switch == "Signal":
                self.Fill = self.FillMCSignal
            elif switch == "WrongPID":
                self.Fill = self.FillMCWrongPID
            else:
                print("Switch {} not recognized (DMesonJetPIDContainer.SetMC)".format(switch))
                exit(1)
        else:
            self.Fill = self.FillData

    def GetHistogramCollection(self, selection_type):
        if selection_type == 1:
            return self.fD0
        elif selection_type == 2:
            return self.fD0bar
        elif selection_type == 3:
            return self.fNoPID
        else:
            return None

    def GetCollectionName(self, selection_type):
        if selection_type == 1:
            return "D0"
        elif selection_type == 2:
            return "D0bar"
        elif selection_type == 3:
            return "NoPID"
        else:
            return None

    def CreateDPtHistogram(self, name):
        hname = "{}DPtHistogram".format(name)
        htilte = "{};#it{{p}}_{{T,D}} (GeV/#it{{c}});yield".format(hname)
        h = ROOT.TH1F(hname, htilte, len(self.fDPtBins) - 1, numpy.array(self.fDPtBins, dtype=numpy.float64))
        return h

    def CreateJetPtHistogram(self, name):
        hname = "{}JetPtHistogram".format(name)
        htilte = "{};#it{{p}}_{{T,ch jet}} (GeV/#it{{c}});yield".format(hname)
        h = ROOT.TH1F(hname, htilte, len(self.fJetPtBins) - 1, numpy.array(self.fJetPtBins, dtype=numpy.float64))
        return h

    def GenerateDataHistograms(self, selection_type):
        histos = self.GetHistogramCollection(selection_type)

        h = self.CreateDPtHistogram(self.GetCollectionName(selection_type))
        histos["DPtHistogram"] = h

        if self.fJetName:
            h = self.CreateJetPtHistogram(self.GetCollectionName(selection_type))
            histos["JetPtHistogram"] = h

    def GenerateMCHistograms(self):
        h = self.CreateDPtHistogram("MCSignal")
        self.fMCSignal["DPtHistogram"] = h

        h = self.CreateDPtHistogram("MCWrongPID")
        self.fMCWrongPID["DPtHistogram"] = h

        if self.fJetName:
            h = self.CreateJetPtHistogram("MCSignal")
            self.fMCSignal["JetPtHistogram"] = h

            h = self.CreateJetPtHistogram("MCWrongPID")
            self.fMCWrongPID["JetPtHistogram"] = h

    def GenerateAllHistograms(self):
        for selection_type in range(1, 4):
            self.GenerateDataHistograms(selection_type)

        if self.fMC:
            self.GenerateMCHistograms()

    def FillData(self, event, eventWeight):
        dmeson = event.DmesonJet
        if dmeson.fPt < self.fMinDPt or dmeson.fPt >= self.fMaxDPt: return

        if self.fJetName:
            jet = getattr(event, self.fJetName)
            if jet.fPt < self.fMinJetPt or jet.fPt >= self.fMaxJetPt: return
            if math.fabs(jet.fEta) > self.fMaxJetEta: return
        else:
            jet = None

        histos = self.GetHistogramCollection(dmeson.fSelectionType)

        hname = "DPtHistogram"
        histos[hname].Fill(dmeson.fPt, eventWeight)

        if jet:
            hname = "JetPtHistogram"
            histos[hname].Fill(jet.fPt, eventWeight)

    def FillMC(self, histos, event, eventWeight):
        dmeson = event.DmesonJet
        if dmeson.fPt < 2 or dmeson.fPt >= 30: return

        if self.fJetName:
            jet = getattr(event, self.fJetName)
            if jet.fPt < 5 or jet.fPt >= 30: return
            if math.fabs(jet.fEta) > 0.5: return
        else:
            jet = None

        hname = "DPtHistogram"
        histos[hname].Fill(dmeson.fPt, eventWeight)

        if jet:
            hname = "JetPtHistogram"
            histos[hname].Fill(jet.fPt, eventWeight)

    def FillMCSignal(self, event, eventWeight):
        self.FillMC(self.fMCSignal, event, eventWeight)

    def FillMCWrongPID(self, event, eventWeight):
        self.FillMC(self.fMCWrongPID, event, eventWeight)


class DMesonJetPIDAnalysis:

    def __init__(self, name, title, mc, trigger, dmeson, dmeson_suffix, jet_pt_bins, d_pt_bins, projector):
        self.fName = name
        self.fTitle = title
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fDMesonSuffix = dmeson_suffix
        self.fProjector = projector
        self.fJetPtBins = jet_pt_bins
        self.fDPtBins = d_pt_bins
        self.fMC = mc
        self.fCanvases = []

    def StartAnalysis(self):
        self.DoProjections()
        self.DoAnalysis()
        self.Plot()

    def SaveRootFile(self, path):
        pass

    def PlotNoPID(self):
        vars = ["JetPtHistogram", "DPtHistogram"]
        for var in vars:
            PIDOK = self.fDataContainer.fD0[var].Clone("PIDOK")
            PIDOK.Add(self.fDataContainer.fD0bar[var])
            PIDOK.SetTitle("PID")
            globalList.append(PIDOK)
            PIDNO = self.fDataContainer.fNoPID[var].Clone("PIDNO")
            PIDNO.SetTitle("No PID")
            globalList.append(PIDNO)
            comp = DMesonJetCompare.DMesonJetCompare("NoPID{}".format(var))
            comp.CompareSpectra(PIDOK, [PIDNO])
            for r in comp.fResults:
                globalList.append(r)
            self.fCanvases.append(comp.fCanvasSpectra)
            self.fCanvases.append(comp.fCanvasRatio)

    def PlotWrongPID(self):
        vars = ["JetPtHistogram", "DPtHistogram"]
        for var in vars:
            signal = self.fDataContainer.fMCSignal[var].Clone("Signal")
            signal.SetTitle("Signal")
            globalList.append(signal)
            wrongPID = self.fDataContainer.fMCWrongPID[var].Clone("WrongPID")
            wrongPID.SetTitle("Wrong PID")
            globalList.append(wrongPID)
            comp = DMesonJetCompare.DMesonJetCompare("WrongPID{}".format(var))
            comp.CompareSpectra(signal, [wrongPID])
            for r in comp.fResults:
                globalList.append(r)
            self.fCanvases.append(comp.fCanvasSpectra)
            self.fCanvases.append(comp.fCanvasRatio)

    def Plot(self):
        self.PlotNoPID()
        self.PlotWrongPID()

    def SavePlots(self, path, format):
        pass

    def DoAnalysis(self):
        self.fDataContainer.DoSumw2()

    def DoProjections(self):
        data_container = DMesonJetPIDContainer(self.fMC, self.fJetPtBins, self.fDPtBins, "Charged", "R040")
        data_container.SetMC(False)
        self.fProjector.StartProjection(self.fTrigger, self.fDMeson, self.fDMesonSuffix, data_container)

        if self.fMC:
            data_container.SetMC("Signal")
            self.fProjector.StartProjection(self.fTrigger, "{}_SignalOnly".format(self.fDMeson), self.fDMesonSuffix, data_container)

            data_container.SetMC("WrongPID")
            self.fProjector.StartProjection(self.fTrigger, "{}_OnlyWrongPIDAccepted".format(self.fDMeson), self.fDMesonSuffix, data_container)

        self.fDataContainer = data_container
