#!/usr/bin/env python
#python base classes and utilities for D Meson jet analysis

import ROOT
import array
import os
import math
import copy
import DMesonJetUtils
from enum import Enum

class DetectorResponse:
    def __init__(self, name, jetName, axis, cuts, weightEff):
        self.fWeightEfficiency = weightEff
        self.fAxis = axis
        self.fCuts = DMesonJetCuts(cuts)
        self.fJetInfo = False
        for a in self.fAxis:
            if "jet" in a.fTruthAxis.fName or a.fTruthAxis.fName == "d_z":
                self.fJetInfo = True
                break
        self.fName = name
        self.fJetName = jetName
        self.fResponseMatrix = self.GenerateResponseMatrix(axis)
        self.fTruth = self.GenerateTruth(axis)
        self.fMeasured = self.GenerateMeasured(axis)
        self.fReconstructedTruth = self.GenerateTruth(axis, "RecontructedTruth")
        self.fEfficiency = None
        self.fResolution = None
        self.fEnergyScaleShift = None
        self.fEnergyScaleShiftMedian = None
        if self.fJetInfo and len(self.fAxis) == 1 and "pt" in self.fAxis[0].fTruthAxis.fName:
            varName = "({0}-{1}) / {1}".format(self.fAxis[0].fDetectorAxis.GetVariableName(), self.fAxis[0].fTruthAxis.GetVariableName())
            sname = "{0}_DetectorResponse".format(self.fName)
            self.fStatistics = DMesonJetUtils.StatisticMultiSet(sname, self.fAxis[0].fTruthAxis, varName)
        else:
            self.fStatistics = None
        if len(self.fAxis) == 2 and "pt" in self.fAxis[0].fTruthAxis.fName:
            self.fResponseMatrix1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fDetectorAxis, [self.fAxis[1].fDetectorAxis, self.fAxis[1].fTruthAxis], bin, "DetectorResponse") for bin in range(0, len(self.fAxis[0].fTruthAxis.fBins)+1)]
            self.fTruth1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fTruthAxis, [self.fAxis[1].fTruthAxis], bin, "Truth") for bin in range(0, len(self.fAxis[0].fTruthAxis.fBins)+1)]
            self.fMeasured1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fDetectorAxis, [self.fAxis[1].fDetectorAxis], bin, "Measured") for bin in range(0, len(self.fAxis[0].fDetectorAxis.fBins)+1)]
            self.fReconstructedTruth1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fDetectorAxis, [self.fAxis[1].fDetectorAxis], bin, "RecontructedTruth") for bin in range(0, len(self.fAxis[0].fDetectorAxis.fBins)+1)]
            self.fEfficiency1D = []
            self.fEfficiency1DRatios = []
        else:
            self.fResponseMatrix1D = None
            self.fTruth1D = None
            self.fMeasured1D = None
            self.fReconstructedTruth1D = None
            self.fEfficiency1D = None
            self.fEfficiency1DRatios = None

    def GenerateResolution(self):
        if not self.fStatistics:
            return
        self.fStatistics.PrintSummary(self.fName)
        self.fEnergyScaleShift = self.fStatistics.GenerateMeanHistogram("EnergyScaleShift")
        self.fResolution = self.fStatistics.GenerateStdDevHistogram("Resolution")
        self.fEnergyScaleShiftMedian = self.fStatistics.GenerateMedianHistogram("EnergyScaleShiftMedian")

    def GenerateLowerDimensionHistogram(self, axisProj, axis, bin, name):
        if bin >= 0 and bin < len(axisProj.fBins):
            binLimits = BinLimits()
            binLimits.AddFromAxis(axisProj, bin)
            binName = binLimits.GetName()
            binTitle = binLimits.GetTitle()
        else:
            binLimits = None
            binName = "NoJet"
            binTitle = "All, no #it{p}_{T,jet} requirement"

        hist = self.GenerateHistogram(axis, "{0}_{1}".format(name, binName))
        hist.SetTitle(binTitle)
        return hist

    def GenerateRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        rlist.Add(self.fResponseMatrix)
        rlist.Add(self.fEfficiency)
        rlist.Add(self.fReconstructedTruth)
        rlist.Add(self.fTruth)
        rlist.Add(self.fMeasured)
        if self.fJetInfo and len(self.fAxis) == 1 and "pt" in self.fAxis[0].fTruthAxis.fName:
            rlist.Add(self.fResolution)
            rlist.Add(self.fEnergyScaleShift)
            rlist.Add(self.fEnergyScaleShiftMedian)
        if len(self.fAxis) == 2 and "pt" in self.fAxis[0].fTruthAxis.fName:
            for eff in self.fEfficiency1D:
                rlist.Add(eff)
            for h in self.fTruth1D:
                rlist.Add(h)
            for h in self.fMeasured1D:
                rlist.Add(h)
            for h in self.fReconstructedTruth1D:
                rlist.Add(h)
            for h in self.fResponseMatrix1D:
                rlist.Add(h)
        return rlist

    def GenerateEfficiency(self):
        #if len(self.fAxis) == 1:
        #self.fEfficiency = ROOT.TGraphAsymmErrors(self.fReconstructedTruth, self.fTruth, "b(1,1) mode")
        #self.fEfficiency.SetName("{0}_Efficiency".format(self.fName))
        #self.fEfficiency.SetTitle("{0}_Efficiency".format(self.fName))
        #else:
        self.fEfficiency = self.GenerateTruth(self.fAxis, "Efficiency")
        self.fEfficiency.Divide(self.fReconstructedTruth, self.fTruth)
        self.fEfficiency.GetXaxis().SetTitle(self.fAxis[0].fTruthAxis.GetTitle())
        if len(self.fAxis) == 1:
            self.fEfficiency.GetYaxis().SetTitle("Efficiency")
        elif len(self.fAxis) == 2:
            self.fEfficiency.GetZaxis().SetTitle("Efficiency")

        if len(self.fAxis) == 2:
            self.GeneratePartialMultiEfficiency()

    def IterateResponseMatrix(self, coord=None, axis=-1, minAxis=0, maxAxis=-1):
        if maxAxis == -1:
            maxAxis = len(self.fAxis)-1
        if coord is None:
            coord = array.array('i',[-1]*(len(self.fAxis)*2))
        if axis == -1:
            axis = minAxis
        for ibin in range(0, self.fResponseMatrix.GetAxis(axis).GetNbins()+2):
            coord[axis] = ibin 
            if axis == maxAxis: 
                yield coord
            else:
                iterNextAxis = self.IterateResponseMatrix(coord, axis+1, minAxis, maxAxis)
                for coord in iterNextAxis:
                    yield coord

    def IterateResponseMatrixMeasured(self, coord=None, axis=-1):
        return self.IterateResponseMatrix(coord, axis, 0, len(self.fAxis)-1)

    def IterateResponseMatrixTruth(self, coord=None, axis=-1):
        return self.IterateResponseMatrix(coord, axis, len(self.fAxis), len(self.fAxis)*2-1)

    def FoldResponse(self, truth):
        meas = self.GenerateMeasured(self.fAxis)

        if len(self.fAxis) == 1:
            for xbin in range(0, self.fResponseMatrix.GetNbinsX()+2):
                binContent = 0
                binError = 0
                for ybin in range(0, self.fResponseMatrix.GetNbinsY()+2):
                    if self.fReconstructedTruth.GetBinContent(ybin) == 0:
                        continue
                    binContent += self.fResponseMatrix.GetBinContent(xbin, ybin) * truth.GetBinContent(ybin) / self.fReconstructedTruth.GetBinContent(ybin)
                    binError += ((self.fResponseMatrix.GetBinError(xbin, ybin) * truth.GetBinContent(ybin)) ** 2 + (truth.GetBinError(ybin) * self.fResponseMatrix.GetBinContent(xbin, ybin)) ** 2) / self.fReconstructedTruth.GetBinContent(ybin) / self.fReconstructedTruth.GetBinContent(ybin)
                binError = math.sqrt(binError)
                meas.SetBinContent(xbin, binContent)
                meas.SetBinError(xbin, binError)
        else:
            measAxisIter = self.IterateResponseMatrixMeasured()
            for coord in measAxisIter:
                binContent = 0
                binError = 0
                norm = 0
                truthAxisIter = self.IterateResponseMatrixTruth(coord)
                for coord in truthAxisIter:
                    normValue = self.GetHistogramBinContent(self.fReconstructedTruth, coord[len(self.fAxis):len(self.fAxis)*2])
                    if normValue == 0:
                        continue
                    respMatrixValue = self.GetResponseMatrixBinContent(coord)
                    respMatrixError = self.GetResponseMatrixBinError(coord)
                    truthValue = self.GetHistogramBinContent(truth, coord[len(self.fAxis):len(self.fAxis)*2])
                    truthError = self.GetHistogramBinError(truth, coord[len(self.fAxis):len(self.fAxis)*2])
                    binContent +=  respMatrixValue * truthValue / normValue
                    binError += (respMatrixError * truthValue) ** 2 + (truthError * respMatrixValue) ** 2 / normValue / normValue

                binError = math.sqrt(binError)
                self.SetHistogramBinContent(meas, coord[0:len(self.fAxis)], binContent)
                self.SetHistogramBinError(meas, coord[0:len(self.fAxis)], binError)

        return meas

    def GetResponseMatrixBinContent(self, coord):
        return self.GetHistogramBinContent(self.fResponseMatrix, coord)

    def GetResponseMatrixBinError(self, coord):
        return self.GetHistogramBinError(self.fResponseMatrix, coord)

    def GetHistogramBinContent(self, histo, coord):
        if len(coord) == 1:
            return histo.GetBinContent(coord[0])
        elif len(coord) == 2:
            return histo.GetBinContent(coord[0], coord[1])
        elif len(coord) == 3:
            return histo.GetBinContent(coord[0], coord[1], coord[3])
        else:
            return histo.GetBinContent(coord)

    def GetHistogramBinError(self, histo, coord):
        if len(coord) == 1:
            return histo.GetBinError(coord[0])
        elif len(coord) == 2:
            return histo.GetBinError(coord[0], coord[1])
        elif len(coord) == 3:
            return histo.GetBinError(coord[0], coord[1], coord[3])
        else:
            return histo.GetBinError(coord)

    def SetHistogramBinContent(self, histo, coord, v):
        if len(coord) == 1:
            return histo.SetBinContent(coord[0], v)
        elif len(coord) == 2:
            return histo.SetBinContent(coord[0], coord[1], v)
        elif len(coord) == 3:
            return histo.SetBinContent(coord[0], coord[1], coord[3], v)
        else:
            return histo.SetBinContent(coord, v)

    def SetHistogramBinError(self, histo, coord, v):
        if len(coord) == 1:
            return histo.SetBinError(coord[0], v)
        elif len(coord) == 2:
            return histo.SetBinError(coord[0], coord[1], v)
        elif len(coord) == 3:
            return histo.SetBinError(coord[0], coord[1], coord[3], v)
        else:
            return histo.SetBinError(coord, v)

    def GeneratePartialMultiEfficiency(self):
        if not self.fTruth1D or not self.fReconstructedTruth1D:
            return
        for bin,(truth,recoTruth) in enumerate(zip(self.fTruth1D, self.fReconstructedTruth1D)):
            eff = self.GeneratePartialMultiEfficiencyForBin(self.fAxis[1], self.fAxis[0].fTruthAxis, bin, truth, recoTruth)
            self.fEfficiency1D.append(eff)
            ratio = self.MakeComparisonRatio(eff, self.fEfficiency1D[0])
            self.fEfficiency1DRatios.append(ratio)

    def MakeComparisonRatio(self, num, den):
        hname = "{0}_Over_{1}".format(num.GetName(), den.GetName()) 
        ratio = num.Clone(hname)
        ratio.Divide(den)
        ratio.GetYaxis().SetTitle("ratio")
        return ratio

    def GeneratePartialMultiEfficiencyForBin(self, axis, axisProj, bin, truth, recoTruth):
        hname = truth.GetName().replace("Truth", "Efficiency")
        #eff = ROOT.TGraphAsymmErrors(recoTruthProj, truthProj, "b(1,1) mode")
        #eff.SetName("{0}_{1}".format(binName, name))

        eff = self.GenerateTruth([axis], "Efficiency")
        eff.SetName(hname)
        eff.SetTitle(truth.GetTitle())
        eff.GetYaxis().SetTitle("Efficiency")
        eff.Divide(recoTruth, truth)

        return eff

    def GenerateResponseMatrix(self, axis):
        haxis = []
        for a in axis:
            haxis.append(a.fDetectorAxis)
        for a in axis:
            haxis.append(a.fTruthAxis)

        return self.GenerateHistogram(haxis, "DetectorResponse")
    
    def GenerateTruth(self, axis, name="Truth"):
        truth = []
        for a in axis:
            truth.append(a.fTruthAxis)

        return self.GenerateHistogram(truth, name)

    def GenerateMeasured(self, axis, name="Measured"):
        detector = []
        for a in axis:
            detector.append(a.fDetectorAxis)

        return self.GenerateHistogram(detector, name)

    def GenerateHistogram(self, axis, name):
        hname = "{0}_{1}".format(self.fName, name)
        if len(axis) == 1:
            hist = ROOT.TH1D(hname, hname, len(axis[0].fBins)-1, array.array('d',axis[0].fBins))
            hist.GetXaxis().SetTitle(axis[0].GetTitle())
            hist.GetYaxis().SetTitle("counts")
        elif len(axis) == 2:
            hist = ROOT.TH2D(hname, hname, len(axis[0].fBins)-1, array.array('d',axis[0].fBins), len(axis[1].fBins)-1, array.array('d',axis[1].fBins))
            hist.GetXaxis().SetTitle(axis[0].GetTitle())
            hist.GetYaxis().SetTitle(axis[1].GetTitle())
            hist.GetZaxis().SetTitle("counts")
        elif len(axis) == 3:
            hist = ROOT.TH3D(hname, hname, len(axis[0].fBins)-1, array.array('d',axis[0].fBins), len(axis[1].fBins)-1, array.array('d',axis[1].fBins), len(axis[2].fBins)-1, array.array('d',axis[2].fBins))
            hist.GetXaxis().SetTitle(axis[0].GetTitle())
            hist.GetYaxis().SetTitle(axis[1].GetTitle())
            hist.GetZaxis().SetTitle(axis[2].GetTitle())
        else:
            nbins = array.array('i', [len(a.fBins) for a in axis])
            hist = ROOT.THnSparseD(hname, hname, len(axis), nbins)
            for i,a in enumerate(axis):
                hist.GetAxis(i).Set(len(a.fBins)-1, array.array('d',a.fBins))
                hist.GetAxis(i).SetTitle(axis[i].GetTitle())

        hist.Sumw2()

        return hist

    def FillResolution(self, recoDmeson, truthDmeson, recoJet, truthJet, w):
        if not self.fStatistics:
            return

        if (not (recoJet and truthJet)) or recoJet.fPt <= 0 or truthJet.fPt <= 0:
            return

        self.fStatistics.Fill(truthJet.fPt, (recoJet.fPt - truthJet.fPt) / truthJet.fPt, w)

    def FillResponseMatrix(self, axis, resp, recoDmeson, truthDmeson, recoJet, truthJet, w):
        naxis = len(axis)
        values = array.array('d', [0]*(naxis*2))
        for i,a in enumerate(axis):
            if a == "jet_pt":
                values[i] = recoJet.fPt
                values[i+naxis] = truthJet.fPt
            elif a == "jet_eta":
                values[i] = recoJet.fEta
                values[i+naxis] = truthJet.fEta
            elif a == "d_pt":
                values[i] = recoDmeson.fPt
                values[i+naxis] = truthDmeson.fPt
            elif a == "d_eta":
                values[i] = recoDmeson.fEta
                values[i+naxis] = truthDmeson.fEta
            elif a == "d_z":
                values[i] = recoJet.fZ
                values[i+naxis] = truthJet.fZ

        self.FillHistogram(resp, values, w)

    def FillSpectrum(self, axis, hist, dmeson, jet, w):
        values = array.array('d')
        for a in axis:
            if a == "jet_pt":
                values.append(jet.fPt)
            elif a == "d_pt":
                values.append(dmeson.fPt)
            elif a == "jet_eta":
                values.append(jet.fEta)
            elif a == "d_eta":
                values.append(dmeson.fEta)
            elif a == "d_z":
                values.append(jet.fZ)

        self.FillHistogram(hist, values, w)

    def FillHistogram(self, hist, values, w):
        if len(values) == 1:
            hist.Fill(values[0], w)
        elif len(values) == 2:
            hist.Fill(values[0], values[1], w)
        elif len(values) == 3:
            hist.Fill(values[0], values[1], values[2], w)
        else:
            hist.Fill(values, w)

    def FillDetectorResponse(self, recoDmeson, truthDmeson, recoJet, truthJet, w):
        if (not recoJet or recoJet.fPt > 0) and \
           (not truthJet or truthJet.fPt > 0):
            self.FillResponseMatrix([a.fDetectorAxis.fName for a in self.fAxis], self.fResponseMatrix, recoDmeson, truthDmeson, recoJet, truthJet, w)

        if self.fResponseMatrix1D:
            self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[self.fResponseMatrix.GetAxis(0).GetNbins()+1], recoDmeson, truthDmeson, recoJet, truthJet, w)
            if self.fAxis[0].fTruthAxis.fName == "jet_pt":
                bin = self.fResponseMatrix.GetAxis(len(self.fAxis)).FindBin(truthJet.fPt)
            elif self.fAxis[0].fTruthAxis.fName == "d_pt":
                bin = self.fResponseMatrix.GetAxis(len(self.fAxis)).FindBin(truthDmeson.fPt)
            else:
                bin = 0
            if bin >= 1 and bin <= self.fResponseMatrix.GetAxis(0).GetNbins():
                self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[bin], recoDmeson, truthDmeson, recoJet, truthJet, w)
                self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[0], recoDmeson, truthDmeson, recoJet, truthJet, w)

    def FillMeasured(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fDetectorAxis.fName for a in self.fAxis], self.fMeasured, dmeson, jet, w)

        if self.fMeasured1D:
            self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[self.fMeasured.GetNbinsX()+1], dmeson, jet, w)
            if self.fAxis[0].fDetectorAxis.fName == "jet_pt":
                bin = self.fMeasured.GetXaxis().FindBin(jet.fPt)
            elif self.fAxis[0].fDetectorAxis.fName == "d_pt":
                bin = self.fMeasured.GetXaxis().FindBin(dmeson.fPt)
            else:
                bin = 0
            if bin >= 1 and bin <= self.fMeasured.GetNbinsX():
                self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[0], dmeson, jet, w)
                self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[bin], dmeson, jet, w)

    def FillTruth(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fTruthAxis.fName for a in self.fAxis], self.fTruth, dmeson, jet, w)

        if self.fTruth1D:
            self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[self.fTruth.GetNbinsX()+1], dmeson, jet, w)
            if self.fAxis[0].fTruthAxis.fName == "jet_pt":
                bin = self.fTruth.GetXaxis().FindBin(jet.fPt)
            elif self.fAxis[0].fTruthAxis.fName == "d_pt":
                bin = self.fTruth.GetXaxis().FindBin(dmeson.fPt)
            else:
                bin = 0
            if bin >= 1 and bin <= self.fTruth.GetNbinsX():
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[bin], dmeson, jet, w)
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[0], dmeson, jet, w)

    def FillRecoTruth(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fTruthAxis.fName for a in self.fAxis], self.fReconstructedTruth, dmeson, jet, w)

        if self.fReconstructedTruth1D:
            self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fReconstructedTruth1D[self.fReconstructedTruth.GetNbinsX()+1], dmeson, jet, w)
            if self.fAxis[0].fTruthAxis.fName == "jet_pt":
                bin = self.fReconstructedTruth.GetXaxis().FindBin(jet.fPt)
            elif self.fAxis[0].fTruthAxis.fName == "d_pt":
                bin = self.fReconstructedTruth.GetXaxis().FindBin(dmeson.fPt)
            else:
                bin = 0
            if bin >= 1 and bin <= self.fReconstructedTruth.GetNbinsX():
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fReconstructedTruth1D[bin], dmeson, jet, w)
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fReconstructedTruth1D[0], dmeson, jet, w)

    def Fill(self, dmeson, w):
        if self.fJetInfo:
            jetTruth = getattr(dmeson, "{0}_truth".format(self.fJetName))
            jetMeasured = getattr(dmeson, "{0}_reco".format(self.fJetName))
        else:
            jetTruth = None
            jetMeasured = None

        dMesonTruth = dmeson.DmesonJet.fGenerated
        dMesonMeasured = dmeson.DmesonJet.fReconstructed

        weff = self.fWeightEfficiency.GetEfficiencyWeight(dMesonMeasured, jetMeasured)

        if dMesonTruth.fPt > 0 and dMesonMeasured.fPt > 0:
            if self.fCuts.ApplyCuts(dMesonMeasured, jetMeasured):
                self.FillDetectorResponse(dMesonMeasured, dMesonTruth, jetMeasured, jetTruth, w*weff)
                self.FillResolution(dMesonMeasured, dMesonTruth, jetMeasured, jetTruth, w*weff)
                self.FillMeasured(dMesonMeasured, jetMeasured, w*weff)
                self.FillRecoTruth(dMesonTruth, jetTruth, w*weff)

        if dMesonTruth.fPt > 0:
            if self.fCuts.ApplyCuts(dMesonTruth, jetTruth):
                self.FillTruth(dMesonTruth, jetTruth, w)

class ResponseAxis:
    def __init__(self, detector, truth):
        self.fDetectorAxis = detector
        self.fTruthAxis = truth

class Axis:
    def __init__(self, name, bins, label = ""):
        self.fName = name
        self.fBins = bins
        self.fLabel = label
        
    def FindBin(self, x):
        for i,(min,max) in enumerate(zip(self.fBins[:-1], self.fBins[1:])):
            if x > min and x < max:
                return i
        return -1

    def GetTitle(self, label = ""):
        varName = self.GetVariableName(label)
        units = self.GetVariableUnits()

        if units:
            title = "{0} ({1})".format(varName, units)
        else:
            title = varName

        return title
        
    def GetVariableUnits(self):
        if "pt" in self.fName:
            return "GeV/#it{c}"
        else:
            return ""

    def GetVariableName(self, label = ""):
        if not label:
            label = self.fLabel

        if label == "nolabel":
            label = ""

        if self.fName == "jet_pt":
            if label:
                title = "#it{{p}}_{{T,ch jet}}^{{ch,{0}}}".format(label)
            else:
                title = "#it{p}_{T,jet}^{ch}"
        elif self.fName == "d_pt":
            if label:
                title = "#it{{p}}_{{T,D}}^{{{0}}}".format(label)
            else:
                title = "#it{p}_{T,D}"
        elif self.fName == "jet_eta":
            if label:
                title = "#it{{#eta}}_{{jet}}^{{{0}}}".format(label)
            else:
                title = "#it{#eta}_{jet}"
        elif self.fName == "d_eta":
            if label:
                title = "#it{{#eta}}_{{D}}^{{{0}}}".format(label)
            else:
                title = "#it{#eta}_{D}"
        elif self.fName == "d_z":
            if label:
                title = "#it{{z}}_{{||,D}}^{{ch,{0}}}".format(label)
            else:
                title = "#it{z}_{||,D}^{ch}"

        return title

class AnalysisType(Enum):
    InvMassFit = 1
    SideBand = 2
    LikeSign = 3

class Spectrum:
    def __init__(self, config, name, binSet=None):
        self.fName = name
        self.fBins = config["bins"]
        self.fHistogram = None
        self.fNormHistogram = None
        self.fUncertainty = None
        self.fMass = None
        self.fMassWidth = None
        self.fBackground = None
        self.fSideBandHistograms = None
        self.fSignalHistograms = None
        self.fAxis = []
        if "title" in config:
            self.fTitle = config["title"]
        else:
            self.fTitle = "" 
        if "side_band" in config:
            self.fAnalysisType = AnalysisType.SideBand
            self.fSideBandMinSigmas = config["side_band"]["min_sigmas"]
            self.fSideBandMaxSigmas = config["side_band"]["max_sigmas"]
            self.fBinCountSignalSigmas = config["side_band"]["max_signal_sigmas"]
        elif "like_sign" in config:
            self.fAnalysisType = AnalysisType.LikeSign
            self.fBinCountSignalSigmas = config["like_sign"]["sigmas"]
            self.fLikeSignTree = config["like_sign"]["name"]
        else:
            self.fAnalysisType = AnalysisType.InvMassFit
        if self.fAnalysisType ==  AnalysisType.SideBand and binSet:
            self.fAxis.append(binSet[self.fBins[0]].fBinCountAnalysisAxis)
        else:
            for axisName, axisBins in config["axis"].iteritems():
                self.fAxis.append(Axis(axisName, axisBins))

    def GenerateRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        if self.fHistogram:
            rlist.Add(self.fHistogram)
        if self.fUncertainty:    
            rlist.Add(self.fUncertainty)
        if self.fMass:    
            rlist.Add(self.fMass)
        if self.fMassWidth:    
            rlist.Add(self.fMassWidth)
        if self.fBackground:    
            rlist.Add(self.fBackground)
        if self.fAnalysisType == AnalysisType.SideBand and self.fSideBandHistograms and self.fSignalHistograms:
            SBlist = ROOT.TList()
            SBlist.SetName("SideBandAnalysis")
            for h in self.fSideBandHistograms:
                SBlist.Add(h)
            for h in self.fSignalHistograms:
                SBlist.Add(h)
            SBlist.Add(self.fSideBandWindowTotalHistogram)
            SBlist.Add(self.fSignalWindowTotalHistogram)
            rlist.Add(SBlist)
        return rlist

    def GenerateNormalizedSpectrum(self, events):
        if not self.fHistogram:
            return
        hname = "{0}_Normalized".format(self.fHistogram.GetName())
        self.fNormHistogram = self.fHistogram.Clone(hname)
        self.fNormHistogram.SetTitle(hname)
        if len(self.fAxis) == 1:
            axisTitle = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d#it{{N}}}}{{d{var}}}".format(var=self.fAxis[0].GetVariableName())
            if self.fAxis[0].GetVariableUnits():
                axisTitle += " ({0})^{{-1}}".format(self.fAxis[0].GetVariableUnits())
            self.fNormHistogram.GetYaxis().SetTitle(axisTitle)
        elif len(self.fAxis) == 2:
            axisTitle = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d^{{2}}#it{{N}}}}{{d{var1} d{var2}}}".format(var1=self.fAxis[0].GetVariableName(), var2=self.fAxis[1].GetVariableName())
            units1 = self.fAxis[0].GetVariableUnits()
            units2 = self.fAxis[1].GetVariableUnits()
            if units1 != units2:
                if units1:
                    axisTitle += " ({0})^{{-1}}".format(units1)
                if units2:
                    axisTitle += " ({0})^{{-1}}".format(units2)
            else:
                if units1:
                    axisTitle += " ({0})^{{-2}}".format(units1)
            self.fNormHistogram.GetZaxis().SetTitle(axisTitle)

        self.fNormHistogram.Scale(1. / events, "width")

class DMesonJetCuts:
    def __init__(self, cutList):
        self.fCuts = cutList

    def ApplyCuts(self, dmeson, jet):
        for cut in self.fCuts:
            if cut["object"] == "d" and dmeson:
                v = getattr(dmeson, cut["variable"])
                if "min" in cut and v < cut["min"]:
                    return False
                if "max" in cut and v > cut["max"]:
                    return False
            if cut["object"] == "jet" and jet:
                v = getattr(jet, cut["variable"])
                if "min" in cut and v < cut["min"]:
                    return False
                if "max" in cut and v > cut["max"]:
                    return False
        return True

class BinMultiSet:
    def __init__(self):
        self.fBinSets = dict()

    def GenerateInvMassRootLists(self):
        for binSet in self.fBinSets.itervalues():
            yield binSet.GenerateInvMassRootList()

    def AddBinSet(self, name, limitSetList, cutList=[], side_band=None, weight=None, fitOptions=""):
        self.fBinSets[name] = BinSet(name, limitSetList, cutList, side_band, weight, fitOptions)

    def FindBin(self, dmeson, jet):
        for binSet in self.fBinSets.itervalues():
            bins = binSet.FindBin(dmeson, jet)
            for bin, w in bins:
                yield bin, w

class BinSet:
    def __init__(self, name, limitSetList, cutList=[], bin_count_axis=None, weight=None, fitOptions=""):
        self.fName = name
        self.fBins = []
        self.fCuts = DMesonJetCuts(cutList)
        self.fFitOptions = fitOptions
        self.fWeightEfficiency = weight
        if bin_count_axis:
            self.fBinCountAnalysisAxis = Axis(bin_count_axis.keys()[0], bin_count_axis.values()[0])
        else:
            self.fBinCountAnalysisAxis = None
        limits = dict()
        self.AddBinsRecursive(limitSetList, limits)

    def GenerateInvMassRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        for bin in self.fBins:
            rlist.Add(bin.fInvMassHisto)
        return rlist

    def AddBinsRecursive(self, limitSetList, limits):
        if len(limitSetList) > 0:
            (limitSetName, limitSet) = limitSetList[0]
            for min,max in zip(limitSet[:-1], limitSet[1:]):
                limits[limitSetName] = min,max
                self.AddBinsRecursive(limitSetList[1:], limits)
        else:
            bin = BinLimits(limits)
            bin.fBinCountAnalysisAxis = self.fBinCountAnalysisAxis
            self.fBins.append(bin)

    def FindBin(self, dmeson, jet):
        w = self.fWeightEfficiency.GetEfficiencyWeight(dmeson, jet)
        if not self.fCuts.ApplyCuts(dmeson, jet):
            return
        for bin in self.fBins:
            if bin.IsInBinLimits(dmeson, jet):
                yield bin, w

class BinLimits:
    def __init__(self, limits=dict()):
        self.fLimits = copy.deepcopy(limits)
        self.fInvMassHisto = None
        self.fMassFitter = None
        self.fBinCountAnalysisAxis = None
        self.fBinCountAnalysisHisto = None

    def FillInvariantMass(self, dmeson, jet, w):
        self.fInvMassHisto.Fill(dmeson.fInvMass, w)
        if self.fBinCountAnalysisHisto:
            if self.fBinCountAnalysisAxis.fName == "jet_pt":
                obsVal = jet.fPt
            else:
                obsVal = 0
            self.fBinCountAnalysisHisto.Fill(dmeson.fInvMass, obsVal, w)
        
    def AddFromAxis(self, axis, binIndex):
        if binIndex == 0:
            min = axis.fBins[0]
            max = axis.fBins[-1]
        else:
            min = axis.fBins[binIndex-1]
            max = axis.fBins[binIndex]                
        self.fLimits[axis.fName] = min, max
      
    def SetMassFitter(self, fitter):
        self.fMassFitter = fitter
    
    def IsSameOf(self, bin):
        for name,(min, max) in self.fLimits.iteritems():
            if bin.fLimits.has_key(name):
                if not bin.fLimits[name] == (min, max):
                    return False
            else:
                return False
        return True

    def SetJetPtLimits(self, min, max):
        self.fLimits["jet_pt"] = min, max
        
    def SetDPtLimits(self, min, max):
        self.fLimits["d_pt"] = min, max
        
    def SetDZLimits(self, min, max):
        self.fLimits["d_z"] = min, max
        
    def SetJetEtaLimits(self, min, max):
        self.fLimits["jet_eta"] = min, max
        
    def SetDEtaLimits(self, min, max):
        self.fLimits["d_eta"] = min, max
        
    def IsInBinLimits(self, dmeson, jet):
        for name,(min,max) in self.fLimits.iteritems():
            if not min < max:
                continue
            if name == "d_pt" and (dmeson.fPt < min or dmeson.fPt >= max):
                return False
            elif name == "jet_pt" and (jet.fPt < min or jet.fPt >= max):
                return False
            elif name == "d_eta" and (dmeson.fEta < min or dmeson.fEta >= max):
                return False
            elif name == "jet_eta" and (jet.fEta < min or jet.fEta >= max):
                return False
            elif name == "d_z" and (jet.fZ < min or jet.fZ >= max):
                return False

        return True
    
    def GetBinCenter(self, axis):
        if axis in self.fLimits:
            (min, max) = self.fLimits[axis]
            return (min + max) / 2
        else:
            return -1

    def GetName(self):
        name = ""
        for varName,(min,max) in self.fLimits.iteritems():
            if varName == "d_pt":
                name += "DPt_{0}_{1}_".format(int(min*100), int(max*100))
            elif varName == "jet_pt":
                name += "JetPt_{0}_{1}_".format(int(min*100), int(max*100))
            elif varName == "d_eta":
                name += "DEta_{0}_{1}_".format(int(min*10), int(max*10))
            elif varName == "jet_eta":
                name += "JetEta_{0}_{1}_".format(int(min*10), int(max*10))
            elif varName == "d_z":
                name += "DZ_{0}_{1}_".format(int(min*100), int(max*100))
        
        #remove last "_"
        if name:
            name = name[:-1]
        return name
        
    def GetTitle(self):
        title = ""
        
        for varName,(min,max) in self.fLimits.iteritems():
            if varName == "d_pt":
                title += "{0:.1f} < #it{{p}}_{{T,D}} < {1:.1f} GeV/#it{{c}}, ".format(min, max)
            elif varName == "jet_pt":
                title += "{0:.0f} < #it{{p}}_{{T,jet}} < {1:.0f} GeV/#it{{c}}, ".format(min, max)
            elif varName == "d_eta":
                title += "{0:.1f} < #it{{#eta}}_{{D}} < {1:.1f}, ".format(min, max)
            elif varName == "jet_eta":
                title += "{0:.0f} < #it{{#eta}}_{{jet}} < {1:.0f}, ".format(min, max)
            elif varName == "d_z":
                title += "{0:.1f} < #it{{z}}_{{||, D}} < {1:.1f}, ".format(min, max)

        #remove last ", "
        if title:
            title = title[:-2]
        return title
    
    def Print(self):
        print(self.GetTitle())
    
    def CreateInvMassHisto(self, trigger, DMesonDef, xAxis, yAxis, nMassBins, minMass, maxMass):
        if trigger:
            hname = "InvMass_{0}_{1}_{2}".format(trigger, DMesonDef, self.GetName())
            htitle = "{0} - {1} Invariant Mass: {2};{3};{4}".format(trigger, DMesonDef, self.GetTitle(), xAxis, yAxis)
            if self.fBinCountAnalysisAxis:
                hnameSB = "InvMassBinCounting_{0}_{1}_{2}".format(trigger, DMesonDef, self.GetName())
                htitleSB = "{0} - {1} Invariant Mass: {2};{3};{4};{5}".format(trigger, DMesonDef, self.GetTitle(), xAxis, self.fBinCountAnalysisAxis.GetTitle(), yAxis)
        else:
            hname = "InvMass_{0}_{1}".format(DMesonDef, self.GetName())
            htitle = "{0} Invariant Mass: {1};{2};{3}".format(DMesonDef, self.GetTitle(), xAxis, yAxis)
            if self.fBinCountAnalysisAxis:
                hnameSB = "InvMassBinCounting_{0}_{1}".format(DMesonDef, self.GetName())
                htitleSB = "{0} Invariant Mass: {1};{2};{3};{4}".format(DMesonDef, self.GetTitle(), xAxis, self.fBinCountAnalysisAxis.GetTitle(), yAxis)
        
        self.fInvMassHisto = ROOT.TH1D(hname, htitle, nMassBins, minMass, maxMass)
        self.fInvMassHisto.Sumw2()
        self.fInvMassHisto.SetMarkerSize(0.9)
        self.fInvMassHisto.SetMarkerStyle(ROOT.kFullCircle)
        self.fInvMassHisto.SetMarkerColor(ROOT.kBlue+2)
        self.fInvMassHisto.SetLineColor(ROOT.kBlue+2)
        
        if self.fBinCountAnalysisAxis:
            self.fBinCountAnalysisHisto = ROOT.TH2D(hname, htitle, nMassBins, minMass, maxMass, len(self.fBinCountAnalysisAxis.fBins)-1, array.array('d', self.fBinCountAnalysisAxis.fBins))
            self.fBinCountAnalysisHisto.Sumw2()
