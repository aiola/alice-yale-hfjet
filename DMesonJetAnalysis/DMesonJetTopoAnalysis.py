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
import collections
import numpy

globalList = []

class DMesonJetVariable:
    def __init__(self, name, varname, vartitle, bins):
        self.fName = name
        self.fVariableName = varname
        self.fBins = bins
        self.fXaxisTitle = vartitle
        self.fHistogram = ROOT.TH1F(self.fName, "{};{};probability".format(self.fName, self.fXaxisTitle), len(bins) - 1, bins)
        self.fHistogram.Sumw2()
        self.fTargetEfficiency = None
        self.fTargetCutValue = None
        self.fEfficiencyAtTarget = None
        self.fCutValueAtTarget = None

    def FillStd(self, dmeson, w):
        v = getattr(dmeson, self.fVariableName)
        self.fHistogram.Fill(v, w)

    def FillAbs(self, dmeson, w):
        v = math.fabs(getattr(dmeson, self.fVariableName))
        self.fHistogram.Fill(v, w)

    def Analyze(self):
        self.fCutEfficiency = self.fHistogram.Clone("{}_CutEfficiency".format(self.fName))
        self.fCutEfficiency.Reset()
        if self.CalculateCutEfficiency == self.CalculateLeftCutEfficiency:
            self.fCutEfficiency.GetYaxis().SetTitle("Efficiency for cut < {}".format(self.fCutEfficiency.GetXaxis().GetTitle()))
        else:
            self.fCutEfficiency.GetYaxis().SetTitle("Efficiency for cut > {}".format(self.fCutEfficiency.GetXaxis().GetTitle()))
        integral = self.fHistogram.Integral(0, -1)
        if integral == 0: return
        self.fHistogram.Scale(1. / integral)
        self.CalculateCutEfficiency()

    def CalculateLeftCutEfficiency(self):
        partialIntErr = 0.
        partialInt = 0.
        for ibin in xrange(0, self.fHistogram.GetNbinsX() + 2):
            partialIntErr = math.sqrt(partialIntErr ** 2 + self.fHistogram.GetBinError(ibin) ** 2)
            partialInt += self.fHistogram.GetBinContent(ibin)
            self.fCutEfficiency.SetBinContent(ibin, partialInt)
            self.fCutEfficiency.SetBinError(ibin, partialIntErr)
            if self.fTargetEfficiency and self.fEfficiencyAtTarget is None:
                if partialInt >= self.fTargetEfficiency:
                    self.fEfficiencyAtTarget = partialInt
                    self.fCutValueAtTarget = self.fHistogram.GetXaxis().GetBinUpEdge(ibin)
            elif self.fTargetCutValue and self.fEfficiencyAtTarget is None:
                if self.fTargetCutValue > self.fHistogram.GetXaxis().GetBinLowEdge(ibin) and self.fTargetCutValue <= self.fHistogram.GetXaxis().GetBinUpEdge(ibin):
                    self.fEfficiencyAtTarget = partialInt
                    self.fCutValueAtTarget = self.fHistogram.GetXaxis().GetBinUpEdge(ibin)

    def CalculateRightCutEfficiency(self):
        partialIntErr = 0.
        partialInt = 0.
        for ibin in reversed(xrange(0, self.fHistogram.GetNbinsX() + 2)):
            partialIntErr = math.sqrt(partialIntErr ** 2 + self.fHistogram.GetBinError(ibin) ** 2)
            partialInt += self.fHistogram.GetBinContent(ibin)
            self.fCutEfficiency.SetBinContent(ibin, partialInt)
            self.fCutEfficiency.SetBinError(ibin, partialIntErr)
            if self.fTargetEfficiency and self.fEfficiencyAtTarget is None:
                if partialInt >= self.fTargetEfficiency:
                    self.fEfficiencyAtTarget = partialInt
                    self.fCutValueAtTarget = self.fHistogram.GetXaxis().GetBinLowEdge(ibin)
            elif self.fTargetCutValue and self.fEfficiencyAtTarget is None:
                if self.fTargetCutValue >= self.fHistogram.GetXaxis().GetBinLowEdge(ibin) and self.fTargetCutValue < self.fHistogram.GetXaxis().GetBinUpEdge(ibin):
                    self.fEfficiencyAtTarget = partialInt
                    self.fCutValueAtTarget = self.fHistogram.GetXaxis().GetBinLowEdge(ibin)

    @classmethod
    def DCA(cls):
        bins = numpy.linspace(0.0, 0.3, 61, True)
        obj = cls("DCA", "fDCA", "DCA (cm)", bins)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def CosThetaStar(cls):
        bins = numpy.linspace(0.0, 1.0, 51, True)
        obj = cls("CosThetaStar", "fCosThetaStar", "|cos(#theta*)|", bins)
        obj.Fill = obj.FillAbs
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def d0d0(cls):
        bins = numpy.linspace(-500e-6, 500e-6, 251, True)
        obj = cls("d0d0", "fd0d0", "#it{d}_{0,#pi}#it{d}_{0,K} (cm^{2})", bins)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def MaxNormd0(cls):
        bins = numpy.linspace(0.0, 10.0, 21, True)
        obj = cls("MaxNormd0", "fMaxNormd0", "max(|#it{d}_{0,#pi}|, |#it{d}_{0,K}|) / #sigma(#it{d}_{0})", bins)
        obj.Fill = obj.FillAbs
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def CosPointing(cls):
        bins = numpy.linspace(0.0, 1.0, 51, True)
        obj = cls("CosPointing", "fCosPointing", "cos(#theta_{p})", bins)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateRightCutEfficiency
        return obj

    @classmethod
    def PtK(cls):
        bins = numpy.concatenate((numpy.linspace(0.0, 5.0, 50), numpy.linspace(5.0, 20.0, 30), numpy.linspace(20.0, 40.0, 5, True)))
        obj = cls("PtK", "fPtK", "#it{p}_{T,K}", bins)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateRightCutEfficiency
        return obj

    @classmethod
    def PtPi(cls):
        bins = numpy.concatenate((numpy.linspace(0.0, 5.0, 50), numpy.linspace(5.0, 20.0, 30), numpy.linspace(20.0, 40.0, 5, True)))
        obj = cls("PtPi", "fPtPi", "#it{p}_{T,#pi}", bins)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateRightCutEfficiency
        return obj

DMesonJetVariable.fgVariableList = [DMesonJetVariable.DCA, DMesonJetVariable.CosThetaStar, DMesonJetVariable.d0d0, DMesonJetVariable.MaxNormd0,
                                    DMesonJetVariable.CosPointing, DMesonJetVariable.PtK, DMesonJetVariable.PtPi]

class DMesonJetTopoContainer:
    def __init__(self, jet_pt_bins, d_pt_bins, jet_type, jet_radius):
        self.fJetPtBins = jet_pt_bins
        self.fDPtBins = d_pt_bins
        self.fVariables = collections.OrderedDict()
        self.fJetType = jet_type
        self.fJetRadius = jet_radius
        for jetPtBin in self.fJetPtBins:
            self.fVariables[jetPtBin] = collections.OrderedDict()
            for dPtBin in self.fDPtBins:
                self.fVariables[jetPtBin][dPtBin] = self.GenerateVariables()

    def AddVariable(self, variables, v):
        variables[v.fName] = v

    def GenerateVariables(self):
        variables = collections.OrderedDict()

        for var in DMesonJetVariable.fgVariableList:
            self.AddVariable(variables, var())

        return variables

    def Fill(self, event, eventWeight):
        dmeson = event.DmesonJet

        if self.fJetType or self.fJetRadius:
            jetName = "Jet_AKT{0}{1}_pt_scheme".format(self.fJetType, self.fJetRadius)
            jet = getattr(event, jetName)
        else:
            print("Error no jet found!")
            exit(1)

        for (jetPtMin, jetPtMax) in self.fJetPtBins:
            if jet.fPt < jetPtMin or jet.fPt >= jetPtMax: continue
            for (dPtMin, dPtMax) in self.fDPtBins:
                if dmeson.fPt < dPtMin or dmeson.fPt >= dPtMax: continue
                for var in self.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)].itervalues():
                    var.Fill(dmeson, eventWeight)

class DMesonJetTopoAnalysis:
    def __init__(self, name, title, trigger, dmeson, jet_pt_bins, d_pt_bins, projector):
        self.fName = name
        self.fTitle = title
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fProjector = projector
        self.fJetPtBins = jet_pt_bins
        self.fDPtBins = d_pt_bins
        self.fTargetEfficiency = None
        self.fTargetCutValueFromAnalysis = None

    def SaveRootFile(self, file):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)

        for (jetPtMin, jetPtMax), jetPtBin in self.fVariables.iteritems():
            for (dPtMin, dPtMax), dPtBin in jetPtBin.iteritems():
                hListName = "Histograms_JetPt{}_{}_DPt_{}_{}".format(jetPtMin, jetPtMax, dPtMin, dPtMax)
                hList = ROOT.TList()
                hList.SetName(hListName)

                eListName = "CutEfficiency_JetPt{}_{}_DPt_{}_{}".format(jetPtMin, jetPtMax, dPtMin, dPtMax)
                eList = ROOT.TList()
                eList.SetName(eListName)

                for v in dPtBin.itervalues():
                    hList.Add(v.fHistogram)
                    eList.Add(v.fCutEfficiency)

                rlist.Add(hList)
                rlist.Add(eList)

        if rlist.GetEntries() > 0:
            file.cd()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)

    def DoProjections(self):
        data_container = DMesonJetTopoContainer(self.fJetPtBins, self.fDPtBins, "Charged", "R040")
        self.fProjector.StartProjection(self.fTrigger, self.fDMeson, None, data_container)
        self.fVariables = data_container.fVariables

    def Analyze(self):
        for (jetPtMin, jetPtMax), jetPtBin in self.fVariables.iteritems():
            for (dPtMin, dPtMax), dPtBin in jetPtBin.iteritems():
                for vname, v in dPtBin.iteritems():
                    if self.fTargetEfficiency:
                        v.fTargetEfficiency = self.fTargetEfficiency
                    elif self.fTargetCutValueFromAnalysis:
                        v.fTargetCutValue = self.fTargetCutValueFromAnalysis.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][vname].fCutValueAtTarget
                    v.Analyze()

class DMesonJetTopoAnalysisManager:
    def __init__(self, dmeson, jet_pt_bin_limits, d_pt_bin_limits):
        self.fDMeson = dmeson
        self.fAnalysisList = collections.OrderedDict()
        self.fJetPtBins = [(minpt, maxpt) for minpt, maxpt in zip(jet_pt_bin_limits[:-1], jet_pt_bin_limits[1:])]
        self.fDPtBins = [(minpt, maxpt) for minpt, maxpt in zip(d_pt_bin_limits[:-1], d_pt_bin_limits[1:])]

    def AddAnalysis(self, name, title, trigger, projector, dmesonSuffix):
        self.fAnalysisList[name] = DMesonJetTopoAnalysis(name, title, trigger, "{}_{}".format(self.fDMeson, dmesonSuffix), self.fJetPtBins, self.fDPtBins, projector)

    def Analyze(self):
        first = None
        for ana in self.fAnalysisList.itervalues():
            if not first:
                ana.fTargetEfficiency = 0.9
                first = ana
            else:
                ana.fTargetCutValueFromAnalysis = first
            ana.Analyze()

    def DoProjections(self):
        for ana in self.fAnalysisList.itervalues():
            ana.DoProjections()

    def StartAnalysis(self):
        self.DoProjections()
        self.Analyze()
        self.Plot()

    def SaveRootFile(self, path):
        fname = "{}/TopoAnalysis.root".format(path)
        file = ROOT.TFile(fname, "recreate")
        if not file or file.IsZombie():
            print("Could not open file '{}'".format(fname))
            return
        for ana in self.fAnalysisList.itervalues():
            ana.SaveRootFile(file)
        print("Results stored in '{}'".format(fname))

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            if c: c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))

    def Plot(self):
        self.fCanvases = []
        self.fCompareObjects = []
        self.fKeepObjects = []
        summary = ""
        for (jetPtMin, jetPtMax) in self.fJetPtBins:
            for (dPtMin, dPtMax) in self.fDPtBins:
                for variableName in [varFunct().fName for varFunct in DMesonJetVariable.fgVariableList]:
                    summary += "******\nJet Pt bin [{}. {}]\nD Pt bin [{}. {}]\nVariable name {}\n******\n".format(jetPtMin, jetPtMax, dPtMin, dPtMax, variableName)
                    hVariable = []
                    hEfficiency = []
                    cutValue = None
                    for ana in self.fAnalysisList.itervalues():
                        cutValue = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fCutValueAtTarget
                        summary += "Name: {}\nCut value: {}\nEfficiency: {}\n\n".format(ana.fTitle, ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fCutValueAtTarget, ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fEfficiencyAtTarget)
                        h = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fHistogram
                        hVar = h.Clone("{}_copy".format(h.GetName()))
                        hVar.SetTitle(ana.fTitle)
                        hVariable.append(hVar)
                        h = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fCutEfficiency
                        hEff = h.Clone("{}_copy".format(h.GetName()))
                        hEff.SetTitle(ana.fTitle)
                        hEfficiency.append(hEff)

                    cname = "{}_JetPt{}_{}_DPt{}_{}".format(hEfficiency[0].GetName(), jetPtMin, jetPtMax, dPtMin, dPtMax)
                    comp = DMesonJetCompare.DMesonJetCompare(cname)
                    comp.fDoSpectraPlot = "lineary"
                    comp.fDoRatioPlot = False
                    comp.CompareSpectra(hEfficiency[0], hEfficiency[1:])
                    if not cutValue is None:
                        line = ROOT.TLine(cutValue, comp.fMainHistogram.GetMinimum(), cutValue, comp.fMainHistogram.GetMaximum())
                        line.SetLineColor(ROOT.kRed)
                        line.SetLineStyle(2)
                        line.SetLineWidth(2)
                        line.Draw()
                    self.fKeepObjects.append(line)
                    self.fCompareObjects.append(comp)
                    self.fCanvases.append(comp.fCanvasSpectra)

                    cname = "{}_JetPt{}_{}_DPt{}_{}".format(hVariable[0].GetName(), jetPtMin, jetPtMax, dPtMin, dPtMax)
                    comp = DMesonJetCompare.DMesonJetCompare(cname)
                    comp.fDoSpectraPlot = "logy"
                    comp.fDoRatioPlot = False
                    comp.CompareSpectra(hVariable[0], hVariable[1:])
                    if not cutValue is None:
                        line = ROOT.TLine(cutValue, comp.fMainHistogram.GetMinimum(), cutValue, comp.fMainHistogram.GetMaximum())
                        line.SetLineColor(ROOT.kRed)
                        line.SetLineStyle(2)
                        line.SetLineWidth(2)
                        line.Draw()
                    self.fKeepObjects.append(line)
                    self.fCompareObjects.append(comp)
                    self.fCanvases.append(comp.fCanvasSpectra)

        print(summary)

