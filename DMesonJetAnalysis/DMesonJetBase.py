#!/usr/bin/env python
#python base classes and utilities for D Meson jet analysis

import ROOT
import array
import os
import math
import copy
import DMesonJetUtils
import DMesonJetProjectors
import numpy
from enum import Enum
import collections
import sys

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
        self.fResponseMatrix = None
        self.fTruth = None
        self.fMeasured = None
        self.fReconstructedTruth = None
        self.fEfficiency = None
        self.fResolution = None
        self.fEnergyScaleShift = None
        self.fEnergyScaleShiftMedian = None
        self.fStatistics = None
        self.fResponseMatrix1D = None
        self.fTruth1D = None
        self.fMeasured1D = None
        self.fReconstructedTruth1D = None
        self.fEfficiency1D = None
        self.fEfficiency1DRatios = None
        self.fResponseMatrixUncertainty = None
    
    def GenerateHistograms(self):
        self.fResponseMatrix = self.GenerateResponseMatrix(self.fAxis)
        if len(self.fAxis) == 1:
            self.fResponseMatrixUncertainty = self.GenerateResponseMatrix(self.fAxis)
            self.fResponseMatrixUncertainty.SetName(self.fResponseMatrixUncertainty.GetName().replace("DetectorResponse", "DetectorResponseUncertainty"))
            self.fResponseMatrixUncertainty.SetTitle(self.fResponseMatrixUncertainty.GetTitle().replace("DetectorResponse", "DetectorResponseUncertainty"))
            self.fResponseMatrixUncertainty.GetZaxis().SetTitle("relative statistical uncertainty")
        self.fTruth = self.GenerateTruth(self.fAxis)
        self.fMeasured = self.GenerateMeasured(self.fAxis)
        self.fReconstructedTruth = self.GenerateTruth(self.fAxis, "RecontructedTruth")
        if self.fJetInfo and len(self.fAxis) == 1 and "pt" in self.fAxis[0].fTruthAxis.fName and self.fAxis[0].fCoarseResponseAxis:
            self.SetupStatistics(self.fAxis[0].fCoarseResponseAxis)
        if len(self.fAxis) == 2 and "pt" in self.fAxis[0].fTruthAxis.fName and self.fAxis[0].fCoarseResponseAxis:
            self.fResponseMatrix1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fCoarseResponseAxis.fDetectorAxis, [self.fAxis[1].fDetectorAxis, self.fAxis[1].fTruthAxis], bin, "DetectorResponse") for bin in range(0, len(self.fAxis[0].fCoarseResponseAxis.fTruthAxis.fBins)+1)]
            self.fTruth1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fCoarseResponseAxis.fTruthAxis, [self.fAxis[1].fTruthAxis], bin, "Truth") for bin in range(0, len(self.fAxis[0].fCoarseResponseAxis.fTruthAxis.fBins)+1)]
            self.fMeasured1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fCoarseResponseAxis.fDetectorAxis, [self.fAxis[1].fDetectorAxis], bin, "Measured") for bin in range(0, len(self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.fBins)+1)]
            self.fReconstructedTruth1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fCoarseResponseAxis.fDetectorAxis, [self.fAxis[1].fDetectorAxis], bin, "RecontructedTruth") for bin in range(0, len(self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.fBins)+1)]
            self.fEfficiency1D = []
            self.fEfficiency1DRatios = []

    def SetupStatistics(self, axis):
        varName = "({0}-{1}) / {1}".format(axis.fDetectorAxis.GetVariableName(), axis.fTruthAxis.GetVariableName())
        sname = "{0}_DetectorResponse".format(self.fName)
        self.fStatistics = DMesonJetUtils.StatisticMultiSet(sname, axis.fTruthAxis, varName)

    @classmethod
    def fromfile(cls, name, jetName, axis, cuts, weightEff, file):
        resp = cls(name, jetName, axis, cuts, weightEff)
        resp.LoadFromRootFile(file)

    def GenerateResolution(self):
        if not self.fStatistics:
            return
        self.fStatistics.PrintSummary(self.fName)
        self.fEnergyScaleShift = self.fStatistics.GenerateMeanHistogram("EnergyScaleShift")
        self.fResolution = self.fStatistics.GenerateStdDevHistogram("Resolution")
        self.fEnergyScaleShiftMedian = self.fStatistics.GenerateMedianHistogram("EnergyScaleShiftMedian")

    def GenerateResponseUncertainty(self):
        if not self.fResponseMatrixUncertainty:
            return
        for x in range(0, self.fResponseMatrix.GetXaxis().GetNbins()+2):
            for y in range(0, self.fResponseMatrix.GetYaxis().GetNbins()+2):
                if self.fResponseMatrix.GetBinContent(x,y) == 0:
                    continue
                self.fResponseMatrixUncertainty.SetBinContent(x, y, self.fResponseMatrix.GetBinError(x,y) / self.fResponseMatrix.GetBinContent(x,y))

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

    def LoadFromRootFile(self, file):
        print("LoadFromRootFile: not yet implemented!")
#         rlist = file.Get(self.fName)
#         self.fResponseMatrix
#         self.fEfficiency
#         self.fReconstructedTruth
#         self.fTruth
#         self.fMeasured
#         if self.fJetInfo and len(self.fAxis) == 1 and "pt" in self.fAxis[0].fTruthAxis.fName:
#             self.fResolution
#             self.fEnergyScaleShift
#             self.fEnergyScaleShiftMedian
#         if len(self.fAxis) == 2 and "pt" in self.fAxis[0].fTruthAxis.fName:
#             self.fEfficiency1D.append()
#             self.fTruth1D.append()
#             self.fMeasured1D.append()
#             self.fReconstructedTruth1D.append()
#             self.fResponseMatrix1D.append()

    def GenerateRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        rlist.Add(self.fResponseMatrix)
        if self.fResponseMatrixUncertainty:
            rlist.Add(self.fResponseMatrixUncertainty)
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
        if self.fStatistics:
            slist = ROOT.TList()
            slist.SetName("DetectorResponse")
            for s in self.fStatistics.fStatisticSets:
                slist.Add(s.fHistogram)
            rlist.Add(slist)
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
        if not self.fTruth1D or not self.fReconstructedTruth1D or not self.fAxis[0].fCoarseResponseAxis:
            return
        for bin,(truth,recoTruth) in enumerate(zip(self.fTruth1D, self.fReconstructedTruth1D)):
            eff = self.GeneratePartialMultiEfficiencyForBin(self.fAxis[1], bin, truth, recoTruth)
            self.fEfficiency1D.append(eff)
            ratio = self.MakeComparisonRatio(eff, self.fEfficiency1D[0])
            self.fEfficiency1DRatios.append(ratio)

    def MakeComparisonRatio(self, num, den):
        hname = "{0}_Over_{1}".format(num.GetName(), den.GetName()) 
        ratio = num.Clone(hname)
        ratio.Divide(den)
        ratio.GetYaxis().SetTitle("ratio")
        return ratio

    def GeneratePartialMultiEfficiencyForBin(self, axis, bin, truth, recoTruth):
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
            self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins()+1], recoDmeson, truthDmeson, recoJet, truthJet, w)
            if self.fAxis[0].fTruthAxis.fName == "jet_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(truthJet.fPt)
            elif self.fAxis[0].fTruthAxis.fName == "d_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(truthDmeson.fPt)
            else:
                bin = 0
            if bin >= 1 and bin <= self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins():
                self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[bin], recoDmeson, truthDmeson, recoJet, truthJet, w)
                self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[0], recoDmeson, truthDmeson, recoJet, truthJet, w)

    def FillMeasured(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fDetectorAxis.fName for a in self.fAxis], self.fMeasured, dmeson, jet, w)

        if self.fMeasured1D:
            self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.GetNbins()+1], dmeson, jet, w)
            if self.fAxis[0].fDetectorAxis.fName == "jet_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.FindBin(jet.fPt)
            elif self.fAxis[0].fDetectorAxis.fName == "d_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.FindBin(dmeson.fPt)
            else:
                bin = 0
            if bin >= 1 and bin <= self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.GetNbins():
                self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[0], dmeson, jet, w)
                self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[bin], dmeson, jet, w)

    def FillTruth(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fTruthAxis.fName for a in self.fAxis], self.fTruth, dmeson, jet, w)

        if self.fTruth1D:
            self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins()+1], dmeson, jet, w)
            if self.fAxis[0].fTruthAxis.fName == "jet_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(jet.fPt)
            elif self.fAxis[0].fTruthAxis.fName == "d_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(dmeson.fPt)
            else:
                bin = 0
            if bin >= 1 and bin <= self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins():
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[bin], dmeson, jet, w)
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[0], dmeson, jet, w)

    def FillRecoTruth(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fTruthAxis.fName for a in self.fAxis], self.fReconstructedTruth, dmeson, jet, w)

        if self.fReconstructedTruth1D:
            self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fReconstructedTruth1D[self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins()+1], dmeson, jet, w)
            if self.fAxis[0].fTruthAxis.fName == "jet_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(jet.fPt)
            elif self.fAxis[0].fTruthAxis.fName == "d_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(dmeson.fPt)
            else:
                bin = 0
            if bin >= 1 and bin <= self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins():
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
        self.fCoarseResponseAxis = None

    def SetCoarseAxis(self, detector, truth):
        self.fCoarseResponseAxis = ResponseAxis(detector, truth)

class Axis:
    def __init__(self, name, bins, label = "", charged_jet = True):
        self.fName = name
        self.fBins = bins
        self.fLabel = label
        self.fChargedJet = charged_jet

    @classmethod
    def fromLimits(cls, name, start, stop, step, label = "", charged_jet = True):
        bins = []
        bins.extend(numpy.linspace(start, stop, (stop-start)/step, True))
        return cls(name, bins, label, charged_jet)

    def FindBin(self, x):
        for i,(min,max) in enumerate(zip(self.fBins[:-1], self.fBins[1:])):
            if x > min and x < max:
                return i
        return -1

    def GetNbins(self):
        return len(self.fBins)-1

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

        if self.fChargedJet:
            jetLabel = "ch jet"
        else:
            jetLabel = "jet"

        if self.fName == "jet_pt":
            if label:
                title = "#it{{p}}_{{T,{jetLabel}}}^{{{label}}}".format(jetLabel=jetLabel, label=label)
            else:
                title = "#it{{p}}_{{T,{jetLabel}}}".format(jetLabel=jetLabel)
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
                title = "#it{{z}}_{{||,D}}^{{{jetLabel},{label}}}".format(jetLabel=jetLabel, label=label)
            else:
                title = "#it{{z}}_{{||,D}}^{{{jetLabel}}}".format(jetLabel=jetLabel)

        return title

class AnalysisType(Enum):
    InvMassFit = 1
    SideBand = 2
    LikeSign = 3
    LikeSignFit = 4
    Truth = 5

class Spectrum:
    def __init__(self, config, dmeson, jtype, jradius, jtitle, binSet, effWeight):
        self.fDMeson = dmeson
        self.fJetType = jtype
        self.fJetRadius = jradius
        self.fJetTitle = jtitle
        if "suffix" in config:
            suffix = config["suffix"]
        else:
            suffix = None 
        self.fName = '_'.join(obj for obj in [self.fDMeson, self.fJetType, self.fJetRadius, config["name"], suffix] if obj)
        self.fBinSet = binSet
        self.fHistogram = None
        self.fNormHistogram = None
        self.fUncertainty = None
        self.fMass = None
        self.fMassWidth = None
        self.fBackground = None
        self.fAxis = []
        self.fSideBandWindowInvMassHistos = dict()
        self.fSignalWindowInvMassHistos = dict()
        self.fSkipBins = None
        self.fEfficiencyWeight = effWeight

        # S-B analysis
        self.fSideBandHistograms = None
        self.fSignalHistograms = None
        self.fSideBandWindowTotalHistogram = None
        self.fSignalWindowTotalHistogram = None

        # L-S analysis
        self.fLikeSignSubtractedBinSet = None
        self.fLikeSignHistograms = None
        self.fUnlikeSignHistograms = None
        self.fLikeSignTotalHistogram = None
        self.fUnlikeSignTotalHistogram = None

        self.fTitle = config["title"]
        if binSet.fTitle:
            if self.fTitle:
                self.fTitle += ", "
            self.fTitle += binSet.fTitle

        if "compare" in config:
            self.fCompare = config["compare"]
        else:
            self.fCompare = True

        if config["type"] == "side_band":
            self.fAnalysisType = AnalysisType.SideBand
            self.fSideBandMinSigmas = config["side_band"]["min_sigmas"]
            self.fSideBandMaxSigmas = config["side_band"]["max_sigmas"]
            self.fBinCountSignalSigmas = config["side_band"]["max_signal_sigmas"]
            self.fBackupSigma = config["side_band"]["backup_sigma"]
            self.fBackupMean = config["side_band"]["backup_mean"]
            if "skip_bins" in config["side_band"]:
                self.fSkipBins = config["side_band"]["skip_bins"]
        elif config["type"] == "like_sign":
            if config["like_sign"]["mode"] == "bin_count":
                self.fAnalysisType = AnalysisType.LikeSign
            elif config["like_sign"]["mode"] == "fit":
                self.fAnalysisType = AnalysisType.LikeSignFit
            else:
                print("Like Sign mode {0} not recognized!!".format(config["like_sign"]["mode"]))
                exit(1)
            self.fSideBandMinSigmas = config["like_sign"]["min_sigmas"]
            self.fSideBandMaxSigmas = config["like_sign"]["max_sigmas"]
            self.fBinCountSignalSigmas = config["like_sign"]["max_signal_sigmas"]
            self.fBinCountNormSignalSigmas = config["like_sign"]["max_signal_norm_sigmas"]
            self.fBackupSigma = config["like_sign"]["backup_sigma"]
            self.fBackupMean = config["like_sign"]["backup_mean"]
            self.fLikeSignTree = config["like_sign"]["name"]
            if "skip_bins" in config["like_sign"]:
                self.fSkipBins = config["like_sign"]["skip_bins"]
        elif config["type"] == "inv_mass_fit":
            self.fAnalysisType = AnalysisType.InvMassFit
        elif config["type"] == "truth":
            self.fAnalysisType = AnalysisType.Truth
        else:
            print("Analysis type {0} not recognized!".format(config["type"]))

        if (self.fAnalysisType == AnalysisType.SideBand or self.fAnalysisType == AnalysisType.LikeSign or self.fAnalysisType == AnalysisType.LikeSignFit):
            self.fAxis.append(binSet.fBinCountAnalysisAxis)
        else:
            for axis in binSet.fAxis:
                self.fAxis.append(axis)

        print("Spectrum {0} with analysis type {1} added".format(self.fName, self.fAnalysisType.name))
        self.BuildHistograms()

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
            for h in self.fSideBandWindowInvMassHistos.itervalues():
                SBlist.Add(h)
            for h in self.fSignalWindowInvMassHistos.itervalues():
                SBlist.Add(h)
            for h in self.fSideBandHistograms:
                SBlist.Add(h)
            for h in self.fSignalHistograms:
                SBlist.Add(h)
            SBlist.Add(self.fSideBandWindowTotalHistogram)
            SBlist.Add(self.fSignalWindowTotalHistogram)
            rlist.Add(SBlist)
        if self.fAnalysisType == AnalysisType.LikeSign and self.fLikeSignSubtractedBinSet:
            LSlist = self.fLikeSignSubtractedBinSet.GenerateInvMassRootList()
            LSlist.SetName("LikeSignAnalysis")
            for h in self.fLikeSignHistograms:
                LSlist.Add(h)
            for h in self.fUnlikeSignHistograms:
                LSlist.Add(h)
            LSlist.Add(self.fLikeSignTotalHistogram)
            LSlist.Add(self.fUnlikeSignTotalHistogram)
            rlist.Add(LSlist)
        return rlist

    def GenerateNormalizedSpectrum(self, events, weighted=False):
        if self.fHistogram:
            self.GenerateNormalizedSpectrumForHistogram(self.fHistogram, events, weighted)
        if self.fFDCorrHistogram:
            self.GenerateNormalizedSpectrumForHistogram(self.fFDCorrHistogram, events, weighted)

    def GenerateNormalizedSpectrumForHistogram(self, hist, events, weighted):
        hname = "{0}_Normalized".format(hist.GetName())
        self.fNormHistogram = hist.Clone(hname)
        self.fNormHistogram.SetTitle(hname)
        if len(self.fAxis) == 1:
            if weighted:
                axisTitle = "#frac{{d#sigma}}{{d{var}}}".format(var=self.fAxis[0].GetVariableName())
                if self.fAxis[0].GetVariableUnits():
                    axisTitle += " [mb ({0})^{{-1}}]".format(self.fAxis[0].GetVariableUnits())
                else:
                    axisTitle += " (mb)"
            else:
                axisTitle = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d#it{{N}}}}{{d{var}}}".format(var=self.fAxis[0].GetVariableName())
                if self.fAxis[0].GetVariableUnits():
                    axisTitle += " ({0})^{{-1}}".format(self.fAxis[0].GetVariableUnits())
            self.fNormHistogram.GetYaxis().SetTitle(axisTitle)
        elif len(self.fAxis) == 2:
            units1 = self.fAxis[0].GetVariableUnits()
            units2 = self.fAxis[1].GetVariableUnits()
            if weighted:
                axisTitle = "#frac{{d^{{2}}#sigma}}{{d{var1} d{var2}}}".format(var1=self.fAxis[0].GetVariableName(), var2=self.fAxis[1].GetVariableName())
                if units1 != units2:
                    axisTitle += " [mb"
                    if units1:
                        axisTitle += " ({0})^{{-1}}".format(units1)
                    if units2:
                        axisTitle += " ({0})^{{-1}}".format(units2)
                    axisTitle += "]"
                else:
                    if units1:
                        axisTitle += " [mb ({0})^{{-2}}]".format(units1)
            else:
                axisTitle = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d^{{2}}#it{{N}}}}{{d{var1} d{var2}}}".format(var1=self.fAxis[0].GetVariableName(), var2=self.fAxis[1].GetVariableName())
                if units1 != units2:
                    if units1:
                        axisTitle += " ({0})^{{-1}}".format(units1)
                    if units2:
                        axisTitle += " ({0})^{{-1}}".format(units2)
                else:
                    if units1:
                        axisTitle += " ({0})^{{-2}}".format(units1)
            self.fNormHistogram.GetZaxis().SetTitle(axisTitle)

        if not weighted and events > 0:
            self.fNormHistogram.Scale(1. / events, "width")
        else:
            self.fNormHistogram.Scale(1, "width")

    def BuildHistograms(self):
        if len(self.fAxis) == 1:
            self.fHistogram = self.BuildHistogram1D(self.fName, "counts")
            self.fUncertainty = self.BuildHistogram1D("{0}_Unc".format(self.fName), "relative statistical uncertainty")
            if self.fAnalysisType == AnalysisType.SideBand:
                self.fSideBandHistograms = []
                self.fSignalHistograms = []
                self.fMass = None
                self.fMassWidth = None
                self.fBackground = self.BuildHistogram1D("{0}_Bkg".format(self.fName), "background |#it{{m}} - <#it{{m}}>| < {0}#sigma".format(int(self.fBinCountSignalSigmas)))
                self.fSideBandWindowTotalHistogram = self.BuildHistogram1D("{0}_SideBandWindowTotal".format(self.fName), "counts")
                self.fSignalWindowTotalHistogram = self.BuildHistogram1D("{0}_SignalWindowTotal".format(self.fName), "counts")
            elif self.fAnalysisType == AnalysisType.LikeSign:
                self.fLikeSignHistograms = []
                self.fUnlikeSignHistograms = []
                self.fMass = None
                self.fMassWidth = None
                self.fBackground = self.BuildHistogram1D("{0}_Bkg".format(self.fName), "background |#it{{m}} - <#it{{m}}>| < {0}#sigma".format(int(self.fBinCountSignalSigmas)))
                self.fLikeSignTotalHistogram = self.BuildHistogram1D("{0}_LikeSignTotal".format(self.fName), "counts")
                self.fUnlikeSignTotalHistogram = self.BuildHistogram1D("{0}_UnlikeSignTotal".format(self.fName), "counts")
            elif self.fAnalysisType == AnalysisType.InvMassFit or self.fAnalysisType == AnalysisType.LikeSignFit:
                self.fMass = self.BuildHistogram1D("{0}_Mass".format(self.fName), "D^{0} mass (GeV/#it{c}^{2})")
                self.fMassWidth = self.BuildHistogram1D("{0}_MassWidth".format(self.fName), "D^{0} mass width (GeV/#it{c}^{2})")
                self.fBackground = self.BuildHistogram1D("{0}_Bkg".format(self.fName), "background |#it{m} - <#it{m}>| < 2#sigma")
        elif len(self.fAxis) == 2:
            self.fHistogram = self.BuildHistogram2D(self.fName, "counts")
            self.fUncertainty = self.BuildHistogram2D("{0}_Unc".format(self.fName), "relative statistical uncertainty")
            if self.fAnalysisType == AnalysisType.InvMassFit:
                self.fMass = self.BuildHistogram2D("{0}_Mass".format(self.fName), "mass")
                self.fBackground = self.BuildHistogram2D("{0}_Bkg".format(self.fName), "background |#it{m} - <#it{m}>| < 3#sigma")
        else:
            print("Cannot handle spectra with more than two axis!")
            return

    def BuildHistogram1D(self, name, yaxis):
        hist = ROOT.TH1D(name, name, len(self.fAxis[0].fBins)-1, array.array('d',self.fAxis[0].fBins))
        hist.GetXaxis().SetTitle(self.fAxis[0].GetTitle())
        hist.GetYaxis().SetTitle(yaxis)
        hist.Sumw2()
        return hist

    def BuildHistogram2D(self, name, zaxis):
        hist = ROOT.TH2D(name, name, len(self.fAxis[0].fBins)-1, array.array('d',self.fAxis[0].fBins), len(self.fAxis[1].fBins)-1, array.array('d',self.fAxis[1].fBins))
        hist.GetXaxis().SetTitle(self.fAxis[0].GetTitle())
        hist.GetYaxis().SetTitle(self.fAxis[1].GetTitle())
        hist.GetZaxis().SetTitle(zaxis)
        hist.Sumw2()
        return hist

    def Fill(self, dmeson, jet, w):
        values = []
        for axis in self.fAxis:
            if axis.fName == "jet_pt":
                values.append(jet.fPt)
            elif axis.fName == "jet_eta":
                values.append(jet.fEta)
            elif axis.fName == "jet_phi":
                values.append(jet.fPhi)
            elif axis.fName == "d_z":
                values.append(jet.fZ)
            elif axis.fName == "d_pt":
                values.append(dmeson.fPt)
            elif axis.fName == "d_eta":
                values.append(dmeson.fEta)
            elif axis.fName == "d_phi":
                values.append(dmeson.fPhi)
        if len(values) == 1:
            self.fHistogram.Fill(values[0], w)
        elif len(values) == 2:
            self.fHistogram.Fill(values[0], values[1], w)
        else:
            print("Cannot handle histograms with more than two axis!")

class DMesonJetCuts:
    def __init__(self, cutList):
        self.fCuts = cutList
        self.InitializeCuts()

    def Dummy(self, obj, min, max):
        return True

    def CutJet(self, dmeson, jet, varname, min, max):
        v = getattr(jet, varname)
        if v < min:
            return False
        if v > max:
            return False
        return True

    def CutJetPt(self, dmeson, jet, dummy, min, max):
        v = jet.fPt
        if v < min:
            return False
        if v > max:
            return False
        return True

    def CutJetEta(self, dmeson, jet, dummy, min, max):
        v = jet.fEta
        if v < min:
            return False
        if v > max:
            return False
        return True

    def CutJetPhi(self, dmeson, jet, dummy, min, max):
        v = jet.fPhi
        if v < min:
            return False
        if v > max:
            return False
        return True

    def CutJetZ(self, dmeson, jet, dummy, min, max):
        v = jet.fZ
        if v < min:
            return False
        if v > max:
            return False
        return True

    def CutD(self, dmeson, jet, varname, min, max): 
        v = getattr(dmeson, varname)
        if v < min:
            return False
        if v > max:
            return False
        return True

    def CutDPt(self, dmeson, jet, dummy, min, max):
        v = dmeson.fPt
        if v < min:
            return False
        if v > max:
            return False
        return True

    def CutDEta(self, dmeson, jet, dummy, min, max):
        v = dmeson.fEta
        if v < min:
            return False
        if v > max:
            return False
        return True

    def CutDPhi(self, dmeson, jet, dummy, min, max):
        v = dmeson.fPhi
        if v < min:
            return False
        if v > max:
            return False
        return True

    def InitializeCuts(self):
        self.fFastCuts = []
        for cut in self.fCuts:
            if "min" in cut: min = cut["min"]
            else: min = -sys.float_info.max
            if "max" in cut: max = cut["max"]
            else: max = +sys.float_info.max
            if cut["object"] == "d":
                if cut["variable"] == "fPt":
                    self.fFastCuts.append((self.CutDPt, None, min, max))
                elif cut["variable"] == "fEta":
                    self.fFastCuts.append((self.CutDEta, None, min, max))
                elif cut["variable"] == "fPhi":
                    self.fFastCuts.append((self.CutDPhi, None, min, max))
                else:
                    self.fFastCuts.append((self.CutD, cut["variable"], min, max))
            elif cut["object"] == "jet":
                if cut["variable"] == "fPt":
                    self.fFastCuts.append((self.CutJetPt, None, min, max))
                elif cut["variable"] == "fEta":
                    self.fFastCuts.append((self.CutJetEta, None, min, max))
                elif cut["variable"] == "fPhi":
                    self.fFastCuts.append((self.CutJetPhi, None, min, max))
                elif cut["variable"] == "fZ":
                    self.fFastCuts.append((self.CutJetZ, None, min, max))
                else:
                    self.fFastCuts.append((self.CutJet, cut["variable"], min, max))

    def ApplyCuts(self, dmeson, jet):
        for funct,var,min,max in self.fFastCuts:
            if not funct(dmeson, jet, var, min, max): return False
        return True

class BinMultiSet:
    def __init__(self):
        self.fBinSets = collections.OrderedDict()

    def GenerateInvMassRootLists(self):
        for binSet in self.fBinSets.itervalues():
            yield binSet.GenerateInvMassRootList()

    def AddBinSet(self, binSet):
        self.fBinSets[binSet.fBinSetName] = binSet
        
    def Initialize(self, dmeson, jtype, jradius, jtitle, inputPath):
        self.fDMeson = dmeson
        self.fJetType = jtype
        self.fJetRadius = jradius
        for k, binSet in self.fBinSets.items():
            r = binSet.Initialize(dmeson, jtype, jradius, jtitle, inputPath)
            if not r:
                del self.fBinSets[k]

    def SetWeightEfficiency(self, weight):
        for binSet in self.fBinSets.itervalues():
            binSet.fWeightEfficiency = weight

    def FindBin(self, dmeson, jet, dmesonDef):
        for binSet in self.fBinSets.itervalues():
            if not dmesonDef in binSet.fNeedInvMass:
                continue
            bins = binSet.FindBin(dmeson, jet)
            for bin, w in bins:
                yield bin, w

    def FindSpectra(self, dmeson, jet):
        for binSet in self.fBinSets.itervalues():
            if not binSet.fCuts.ApplyCuts(dmeson, jet):
                continue
            for s in binSet.fSpectra.itervalues():
                if s.fAnalysisType != AnalysisType.Truth:
                    continue
                yield s, binSet.fWeightEfficiency.GetEfficiencyWeight(dmeson, jet)

    def FindAllSpectra(self):
        for binSet in self.fBinSets.itervalues():
            for s in binSet.fSpectra.itervalues():
                yield s

class BinSet:
    def __init__(self, name, title, need_inv_mass, limitSetList, spectra, axis, cutList=[], bin_count_axis=None, weight=None, fitOptions=""):
        self.fBinSetName = name
        self.fTitle = title
        self.fBins = []
        self.fSpectraConfigs = spectra
        self.fSpectra = collections.OrderedDict()
        self.fAxis = axis
        self.fCuts = DMesonJetCuts(cutList)
        self.fFitOptions = fitOptions
        self.fWeightEfficiency = weight
        self.fNeedInvMass = need_inv_mass
        if bin_count_axis:
            self.fBinCountAnalysisAxis = Axis(bin_count_axis.keys()[0], bin_count_axis.values()[0], "", True)
        else:
            self.fBinCountAnalysisAxis = None
        limits = dict()
        self.AddBinsRecursive(limitSetList, limits)

    def AmIJetDependent(self):
        for axis in self.fAxis:
            if "jet" in axis.fName or axis.fName == "d_z":
                return True
        if self.fBinCountAnalysisAxis and ("jet" in self.fBinCountAnalysisAxis.fName or self.fBinCountAnalysisAxis.fName == "d_z"):
            return True
        for cut in self.fCuts.fCuts:
            if cut["object"] == "jet":
                return True
        return False

    def Initialize(self, dmeson, jtype, jradius, jtitle, inputPath):
        self.fName = "_".join(obj for obj in [dmeson, jtype, jradius, self.fBinSetName] if not obj == None)
        jetDep = self.AmIJetDependent()
        if jtype or jradius:
            if not jetDep:
                return False
        else:
            if jetDep:
                return False
        if jtype == "Full":
            if self.fBinCountAnalysisAxis:
                self.fBinCountAnalysisAxis.fChargedJet = False
            for a in self.fAxis:
                a.fChargedJet = False
        for s in self.fSpectraConfigs:
            if not dmeson in s["active_mesons"]:
                continue
            if "efficiency" in s and s["efficiency"] and not "MCTruth" in dmeson:
                effWeight = DMesonJetProjectors.EfficiencyWeightCalculator("{0}/{1}".format(inputPath, s["efficiency"]["file_name"]), s["efficiency"]["list_name"], s["efficiency"]["object_name"])
            else:
                effWeight = DMesonJetProjectors.SimpleWeight()
            spectrum = Spectrum(s, dmeson, jtype, jradius, jtitle, self, effWeight)
            self.fSpectra[spectrum.fName] = spectrum
        if "MCTruth" in dmeson and not isinstance(self.fWeightEfficiency , DMesonJetProjectors.SimpleWeight):
            self.fWeightEfficiency = DMesonJetProjectors.SimpleWeight()
        return True

    def GenerateInvMassRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        for bin in self.fBins:
            if bin.fInvMassHisto:
                rlist.Add(bin.fInvMassHisto)
            if bin.fMassFitter:
                rlist.Add(bin.fMassFitter)
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
        if not self.fCuts.ApplyCuts(dmeson, jet):
            return
        else:
            w = self.fWeightEfficiency.GetEfficiencyWeight(dmeson, jet)
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
        self.fCounts = 0
        self.fSumw2 = 0

    def FillInvariantMass(self, dmeson, jet, w):
        self.fCounts += w
        self.fSumw2 += w*w
        if self.fInvMassHisto:
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
        
        if not "MCTruth" in DMesonDef:
            self.fInvMassHisto = ROOT.TH1D(hname, htitle, nMassBins, minMass, maxMass)
            self.fInvMassHisto.Sumw2()
            self.fInvMassHisto.SetMarkerSize(0.9)
            self.fInvMassHisto.SetMarkerStyle(ROOT.kFullCircle)
            self.fInvMassHisto.SetMarkerColor(ROOT.kBlue+2)
            self.fInvMassHisto.SetLineColor(ROOT.kBlue+2)
        
        if self.fBinCountAnalysisAxis:
            self.fBinCountAnalysisHisto = ROOT.TH2D(hnameSB, htitleSB, nMassBins, minMass, maxMass, len(self.fBinCountAnalysisAxis.fBins)-1, array.array('d', self.fBinCountAnalysisAxis.fBins))
            self.fBinCountAnalysisHisto.Sumw2()
