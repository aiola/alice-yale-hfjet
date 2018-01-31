#!/usr/bin/env python
# python Class for detector response

import array
import math
import copy

import ROOT

import DMesonJetCuts
import StatisticSet
import BinSet


class DetectorResponse:

    def __init__(self, name, jetName, axis, cuts, weightEff, truthWeight):
        self.fWeightEfficiency = weightEff
        self.fTruthWeight = truthWeight
        self.fAxis = axis
        self.fCuts = DMesonJetCuts.DMesonJetCuts(cuts)
        self.fJetInfo = False
        print("Created response matrix {0} with following axis:".format(name))
        for a in self.fAxis:
            print("{0}: {1}".format(a.fTruthAxis.fName, a.fTruthAxis.fBins))
            if "jet" in a.fTruthAxis.fName or a.fTruthAxis.fName == "d_z": self.fJetInfo = True
        for c in cuts:
            if c["object"] == "jet": self.fJetInfo = True
        self.fName = name
        self.fJetName = jetName
        self.fJetTruthName = "{0}_truth".format(self.fJetName)
        self.fJetRecoName = "{0}_reco".format(self.fJetName)
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

    @classmethod
    def fromfile(cls, name, jetName, axis, cuts, weightEff, file):
        resp = cls(name, jetName, axis, cuts, weightEff)
        resp.LoadFromRootFile(file)

    def GenerateHistograms(self):
        self.fResponseMatrix = self.GenerateResponseMatrix(self.fAxis)
        if len(self.fAxis) == 1:
            self.fResponseMatrixUncertainty = self.GenerateResponseMatrix(self.fAxis)
            self.fResponseMatrixUncertainty.SetName(self.fResponseMatrixUncertainty.GetName().replace("DetectorResponse", "DetectorResponseUncertainty"))
            self.fResponseMatrixUncertainty.SetTitle(self.fResponseMatrixUncertainty.GetTitle().replace("DetectorResponse", "DetectorResponseUncertainty"))
            self.fResponseMatrixUncertainty.GetZaxis().SetTitle("relative statistical uncertainty")
        self.fTruth = self.GenerateTruth(self.fAxis)
        self.fMeasured = self.GenerateMeasured(self.fAxis)
        self.fReconstructedTruth = self.GenerateTruth(self.fAxis, "ReconstructedTruth")
        if self.fJetInfo and len(self.fAxis) == 1 and self.fAxis[0].fCoarseResponseAxis:
            self.SetupStatistics(self.fAxis[0].fCoarseResponseAxis)
        if len(self.fAxis) == 2 and self.fAxis[0].fCoarseResponseAxis and ("jet" in self.fAxis[0].fDetectorAxis.fName or self.fAxis[0].fDetectorAxis.fName == "d_z"):
            self.fResponseMatrix1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fCoarseResponseAxis.fDetectorAxis, [self.fAxis[1].fDetectorAxis, self.fAxis[1].fTruthAxis], bin, "DetectorResponse") for bin in xrange(0, len(self.fAxis[0].fCoarseResponseAxis.fTruthAxis.fBins) + 1)]
            self.fTruth1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fCoarseResponseAxis.fTruthAxis, [self.fAxis[1].fTruthAxis], bin, "Truth") for bin in xrange(0, len(self.fAxis[0].fCoarseResponseAxis.fTruthAxis.fBins) + 1)]
            self.fMeasured1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fCoarseResponseAxis.fDetectorAxis, [self.fAxis[1].fDetectorAxis], bin, "Measured") for bin in xrange(0, len(self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.fBins) + 1)]
            self.fReconstructedTruth1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fCoarseResponseAxis.fTruthAxis, [self.fAxis[1].fTruthAxis], bin, "ReconstructedTruth") for bin in xrange(0, len(self.fAxis[0].fCoarseResponseAxis.fTruthAxis.fBins) + 1)]
            self.fEfficiency1D = []
            self.fEfficiency1DRatios = []

    def SetupStatistics(self, axis):
        varName = "({0}-{1}) / {1}".format(axis.fDetectorAxis.GetVariableName(), axis.fTruthAxis.GetVariableName())
        sname = "{0}_DetectorResponse".format(self.fName)
        self.fStatistics = StatisticSet.StatisticMultiSet(sname, axis.fTruthAxis, varName)

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
        for x in xrange(0, self.fResponseMatrix.GetXaxis().GetNbins() + 2):
            for y in xrange(0, self.fResponseMatrix.GetYaxis().GetNbins() + 2):
                if self.fResponseMatrix.GetBinContent(x, y) == 0:
                    continue
                self.fResponseMatrixUncertainty.SetBinContent(x, y, self.fResponseMatrix.GetBinError(x, y) / self.fResponseMatrix.GetBinContent(x, y))

    def GenerateLowerDimensionHistogram(self, axisProj, axis, bin, name):
        if bin >= 0 and bin < len(axisProj.fBins):
            binLimits = BinSet.BinLimits()
            binLimits.AddFromAxis(axisProj, bin)
            binName = binLimits.GetName()
            binTitle = binLimits.GetTitle()
        else:
            binLimits = None
            binName = "NoJet"
            binTitle = "All, no jet requirement"

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
        if self.fResolution: rlist.Add(self.fResolution)
        if self.fEnergyScaleShift: rlist.Add(self.fEnergyScaleShift)
        if self.fEnergyScaleShiftMedian: rlist.Add(self.fEnergyScaleShiftMedian)
        if len(self.fAxis) == 2 and self.fAxis[0].fCoarseResponseAxis and ("jet" in self.fAxis[0].fDetectorAxis.fName or self.fAxis[0].fDetectorAxis.fName == "d_z"):
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
        # if len(self.fAxis) == 1:
        # self.fEfficiency = ROOT.TGraphAsymmErrors(self.fReconstructedTruth, self.fTruth, "b(1,1) mode")
        # self.fEfficiency.SetName("{0}_Efficiency".format(self.fName))
        # self.fEfficiency.SetTitle("{0}_Efficiency".format(self.fName))
        # else:
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
            maxAxis = len(self.fAxis) - 1
        if coord is None:
            coord = array.array('i', [-1] * (len(self.fAxis) * 2))
        if axis == -1:
            axis = minAxis
        for ibin in xrange(0, self.fResponseMatrix.GetAxis(axis).GetNbins() + 2):
            coord[axis] = ibin
            if axis == maxAxis:
                yield coord
            else:
                iterNextAxis = self.IterateResponseMatrix(coord, axis + 1, minAxis, maxAxis)
                for coord in iterNextAxis:
                    yield coord

    def IterateResponseMatrixMeasured(self, coord=None, axis=-1):
        return self.IterateResponseMatrix(coord, axis, 0, len(self.fAxis) - 1)

    def IterateResponseMatrixTruth(self, coord=None, axis=-1):
        return self.IterateResponseMatrix(coord, axis, len(self.fAxis), len(self.fAxis) * 2 - 1)

    def FoldResponse(self, truth):
        meas = self.GenerateMeasured(self.fAxis)

        if len(self.fAxis) == 1:
            for xbin in xrange(0, self.fResponseMatrix.GetNbinsX() + 2):
                binContent = 0
                binError = 0
                for ybin in xrange(0, self.fResponseMatrix.GetNbinsY() + 2):
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
                    normValue = self.GetHistogramBinContent(self.fReconstructedTruth, coord[len(self.fAxis):len(self.fAxis) * 2])
                    if normValue == 0:
                        continue
                    respMatrixValue = self.GetResponseMatrixBinContent(coord)
                    respMatrixError = self.GetResponseMatrixBinError(coord)
                    truthValue = self.GetHistogramBinContent(truth, coord[len(self.fAxis):len(self.fAxis) * 2])
                    truthError = self.GetHistogramBinError(truth, coord[len(self.fAxis):len(self.fAxis) * 2])
                    binContent += respMatrixValue * truthValue / normValue
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
        for bin, (truth, recoTruth) in enumerate(zip(self.fTruth1D, self.fReconstructedTruth1D)):
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
        # eff = ROOT.TGraphAsymmErrors(recoTruthProj, truthProj, "b(1,1) mode")
        # eff.SetName("{0}_{1}".format(binName, name))

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
            hist = ROOT.TH1D(hname, hname, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins))
            hist.GetXaxis().SetTitle(axis[0].GetTitle())
            hist.GetYaxis().SetTitle("counts")
        elif len(axis) == 2:
            hist = ROOT.TH2D(hname, hname, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins), len(axis[1].fBins) - 1, array.array('d', axis[1].fBins))
            hist.GetXaxis().SetTitle(axis[0].GetTitle())
            hist.GetYaxis().SetTitle(axis[1].GetTitle())
            hist.GetZaxis().SetTitle("counts")
        elif len(axis) == 3:
            hist = ROOT.TH3D(hname, hname, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins), len(axis[1].fBins) - 1, array.array('d', axis[1].fBins), len(axis[2].fBins) - 1, array.array('d', axis[2].fBins))
            hist.GetXaxis().SetTitle(axis[0].GetTitle())
            hist.GetYaxis().SetTitle(axis[1].GetTitle())
            hist.GetZaxis().SetTitle(axis[2].GetTitle())
        else:
            nbins = array.array('i', [len(a.fBins) for a in axis])
            hist = ROOT.THnSparseD(hname, hname, len(axis), nbins)
            for i, a in enumerate(axis):
                hist.GetAxis(i).Set(len(a.fBins) - 1, array.array('d', a.fBins))
                hist.GetAxis(i).SetTitle(axis[i].GetTitle())

        hist.Sumw2()

        return hist

    def FillResolution(self, recoDmeson, truthDmeson, recoJet, truthJet, w):
        if not self.fStatistics: return

        if self.fStatistics.fAxis.fName == "jet_pt":
            if (not (recoJet and truthJet)) or recoJet.fPt <= 0 or truthJet.fPt <= 0: return
            self.fStatistics.Fill(truthJet.fPt, (recoJet.fPt - truthJet.fPt) / truthJet.fPt, w)
        elif self.fStatistics.fAxis.fName == "d_z":
            if (not (recoJet and truthJet)) or recoJet.fZ <= 0 or truthJet.fZ <= 0: return
            self.fStatistics.Fill(truthJet.fZ, (recoJet.fZ - truthJet.fZ) / truthJet.fZ, w)
        elif self.fStatistics.fAxis.fName == "d_pt":
            if (not (recoDmeson and truthDmeson)) or recoDmeson.fPt <= 0 or truthDmeson.fPt <= 0: return
            self.fStatistics.Fill(truthDmeson.fPt, (recoDmeson.fPt - truthDmeson.fPt) / truthDmeson.fPt, w)

    def FillResponseMatrix(self, axis, resp, recoDmeson, truthDmeson, recoJet, truthJet, w):
        naxis = len(axis)
        values = array.array('d', [0] * (naxis * 2))
        for i, a in enumerate(axis):
            if a == "jet_pt":
                values[i] = recoJet.fPt
                values[i + naxis] = truthJet.fPt
            elif a == "jet_eta":
                values[i] = recoJet.fEta
                values[i + naxis] = truthJet.fEta
            elif a == "d_pt":
                values[i] = recoDmeson.fPt
                values[i + naxis] = truthDmeson.fPt
            elif a == "d_eta":
                values[i] = recoDmeson.fEta
                values[i + naxis] = truthDmeson.fEta
            elif a == "d_z":
                values[i] = recoJet.fZ
                values[i + naxis] = truthJet.fZ

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
        dimension = len(values)
        if dimension == 1:
            hist.Fill(values[0], w)
        elif dimension == 2:
            hist.Fill(values[0], values[1], w)
        elif dimension == 3:
            hist.Fill(values[0], values[1], values[2], w)
        else:
            hist.Fill(values, w)

    def FillDetectorResponse(self, recoDmeson, truthDmeson, recoJet, truthJet, w):
        if (not recoJet or recoJet.fPt > 0) and \
           (not truthJet or truthJet.fPt > 0):
            self.FillResponseMatrix([a.fDetectorAxis.fName for a in self.fAxis], self.fResponseMatrix, recoDmeson, truthDmeson, recoJet, truthJet, w)

        if self.fResponseMatrix1D:
            self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins() + 1], recoDmeson, truthDmeson, recoJet, truthJet, w)
            if self.fAxis[0].fDetectorAxis.fName == "jet_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.FindBin(recoJet.fPt) + 1
            elif self.fAxis[0].fDetectorAxis.fName == "d_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.FindBin(recoDmeson.fPt) + 1
            elif self.fAxis[0].fDetectorAxis.fName == "d_z":
                bin = self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.FindBin(recoJet.fZ) + 1
            else:
                bin = 0
            if bin >= 1 and bin <= self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins():
                self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[bin], recoDmeson, truthDmeson, recoJet, truthJet, w)
                self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[0], recoDmeson, truthDmeson, recoJet, truthJet, w)

    def FillMeasured(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fDetectorAxis.fName for a in self.fAxis], self.fMeasured, dmeson, jet, w)

        if self.fMeasured1D:
            self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.GetNbins() + 1], dmeson, jet, w)
            if self.fAxis[0].fDetectorAxis.fName == "jet_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.FindBin(jet.fPt) + 1
            elif self.fAxis[0].fDetectorAxis.fName == "d_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.FindBin(dmeson.fPt) + 1
            elif self.fAxis[0].fDetectorAxis.fName == "d_z":
                bin = self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.FindBin(jet.fZ) + 1
            else:
                bin = 0
            if bin >= 1 and bin <= self.fAxis[0].fCoarseResponseAxis.fDetectorAxis.GetNbins():
                self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[0], dmeson, jet, w)
                self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[bin], dmeson, jet, w)

    def FillTruth(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fTruthAxis.fName for a in self.fAxis], self.fTruth, dmeson, jet, w)

        if self.fTruth1D:
            self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins() + 1], dmeson, jet, w)
            if self.fAxis[0].fTruthAxis.fName == "jet_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(jet.fPt) + 1
            elif self.fAxis[0].fTruthAxis.fName == "d_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(dmeson.fPt) + 1
            elif self.fAxis[0].fTruthAxis.fName == "d_z":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(jet.fZ) + 1
            else:
                bin = 0
            if bin >= 1 and bin <= self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins():
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[bin], dmeson, jet, w)
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[0], dmeson, jet, w)

    def FillRecoTruth(self, recoDmeson, truthDmeson, recoJet, truthJet, w):
        if not truthJet or truthJet.fPt > 0:
            self.FillSpectrum([a.fTruthAxis.fName for a in self.fAxis], self.fReconstructedTruth, truthDmeson, truthJet, w)

        # The following code needs a bit of an explanation.
        # The "partial efficiency" that is calculated by this piece of code is the
        # efficiency of D-jet that are taken from a sub-sample in a certain range of a GENERATED variable.
        # In this way one can verify whether, e.g. a different TRUE jet momentum affects the reconstruction
        # efficiency of a D-jet as a function of D-meson momentum.
        # The resulting "partial efficiencies" DO NOT take into account feed-in and feed-out of the considered range:
        # it just refelcts the efficiency of a sub-sample, regardless of whether the object is reconstructed within that
        # same considered range in the reconstructed variable counterpart (this would be achieved by using reconstructed
        # variables in the following lines instead of generated variables). This is therefore different from setting
        # this cuts in the "cuts" section of the YAML file. Actually this allows to disentangle between differences in the
        # reconstruction efficiencies for different sub-samples (obtained from the code below) and differences in the reconstruction
        # efficiency coming from feed-in and feed-out of the considered range.
        if self.fReconstructedTruth1D:
            self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fReconstructedTruth1D[self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins() + 1], truthDmeson, truthJet, w)
            if self.fAxis[0].fDetectorAxis.fName == "jet_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(truthJet.fPt) + 1
            elif self.fAxis[0].fDetectorAxis.fName == "d_pt":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(truthDmeson.fPt) + 1
            elif self.fAxis[0].fDetectorAxis.fName == "d_z":
                bin = self.fAxis[0].fCoarseResponseAxis.fTruthAxis.FindBin(truthJet.fZ) + 1
            else:
                bin = 0
            if bin >= 1 and bin <= self.fAxis[0].fCoarseResponseAxis.fTruthAxis.GetNbins():
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fReconstructedTruth1D[bin], truthDmeson, truthJet, w)
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fReconstructedTruth1D[0], truthDmeson, truthJet, w)

    def Fill(self, dmeson, w):
        if self.fJetInfo:
            jetTruth = getattr(dmeson, self.fJetTruthName)
            jetMeasured = getattr(dmeson, self.fJetRecoName)
        else:
            jetTruth = None
            jetMeasured = None

        dMesonTruth = dmeson.DmesonJet.fGenerated
        dMesonMeasured = dmeson.DmesonJet.fReconstructed

        if self.fTruthWeight: w *= self.fTruthWeight.GetEfficiencyWeight(dMesonTruth, jetTruth)

        weff = self.fWeightEfficiency.GetEfficiencyWeight(dMesonMeasured, jetMeasured)

        if dMesonTruth.fPt <= 0: return
        if not self.fCuts.ApplyCutsGeneratorOnly(dMesonTruth, jetTruth): return

        if dMesonMeasured.fPt > 0 and self.fCuts.ApplyCuts(dMesonMeasured, jetMeasured):
            self.FillDetectorResponse(dMesonMeasured, dMesonTruth, jetMeasured, jetTruth, w * weff)
            self.FillRecoTruth(dMesonMeasured, dMesonTruth, jetMeasured, jetTruth, w * weff)
            self.FillResolution(dMesonMeasured, dMesonTruth, jetMeasured, jetTruth, w * weff)
            self.FillMeasured(dMesonMeasured, jetMeasured, w * weff)

        if self.fCuts.ApplyCuts(dMesonTruth, jetTruth):
            self.FillTruth(dMesonTruth, jetTruth, w)


class ResponseAxis:

    def __init__(self, detector, truth):
        self.fDetectorAxis = detector
        self.fTruthAxis = truth
        self.fCoarseResponseAxis = None

    def SetCoarseAxis(self, detector, truth):
        self.fCoarseResponseAxis = ResponseAxis(detector, truth)
