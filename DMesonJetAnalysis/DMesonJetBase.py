#!/usr/bin/env python
#python base classes and utilities for D Meson jet analysis

import ROOT
from array import array
import os
import math

def find_file(path, file_name):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file == file_name:
                yield os.path.join(root, file)

class DetectorResponse:
    def __init__(self, name, jetName, axis, weightEff):
        self.fWeightEfficiency = weightEff
        self.fAxis = axis
        self.fJetInfo = False
        for a in self.fAxis:
            if a.fTruthAxis.fName == "jet_pt" or a.fTruthAxis.fName == "d_z":
                self.fJetInfo = True
                break
        self.fName = name
        self.fJetName = jetName
        self.fResponseMatrix = self.GenerateResponseMatrix(axis)
        self.fTruth = self.GenerateTruth(axis)
        self.fMeasured = self.GenerateMeasured(axis)
        self.fReconstructedTruth = self.GenerateTruth(axis, "RecontructedTruth")
        self.fFoldedReconstructedTruth = None
        self.fFoldedTruth = None
        self.fEfficiency = None
        self.fFoldedEfficiency = None
        if len(self.fAxis) == 2:
            self.fResponseMatrix1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fDetectorAxis, [self.fAxis[1].fDetectorAxis, self.fAxis[1].fTruthAxis], bin, "DetectorResponse") for bin in range(0, len(self.fAxis[0].fTruthAxis.fBins)+1)]
            self.fTruth1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fTruthAxis, [self.fAxis[1].fTruthAxis], bin, "Truth") for bin in range(0, len(self.fAxis[0].fTruthAxis.fBins)+1)]
            self.fMeasured1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fDetectorAxis, [self.fAxis[1].fDetectorAxis], bin, "Measured") for bin in range(0, len(self.fAxis[0].fDetectorAxis.fBins)+1)]
            self.fReconstructedTruth1D = [self.GenerateLowerDimensionHistogram(self.fAxis[0].fDetectorAxis, [self.fAxis[1].fDetectorAxis], bin, "RecontructedTruth") for bin in range(0, len(self.fAxis[0].fDetectorAxis.fBins)+1)]
            self.fFoldedTruth1D = None
            self.fEfficiency1D = []
            self.fFoldedEfficiency1D = []
        else:
            self.fResponseMatrix1D = None
            self.fTruth1D = None
            self.fMeasured1D = None
            self.fReconstructedTruth1D = None
            self.fFoldedTruth1D = None
            self.fEfficiency1D = None
            self.fFoldedEfficiency1D = None
            
    def GenerateLowerDimensionHistogram(self, axisProj, axis, bin, name):
        if bin >= 0 and bin < len(axisProj.fBins):
            binLimits = BinLimits()
            binLimits.SetFromAxis(axisProj, bin)
            binName = binLimits.GetName()
            binTitle = binLimits.GetTitle()
        else:
            binLimits = None
            binName = "NoJet"
            binTitle = "All, no jet requirement"

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
        rlist.Add(self.fFoldedTruth)
        rlist.Add(self.fFoldedEfficiency)
        rlist.Add(self.fFoldedReconstructedTruth)
        if len(self.fAxis) == 2:
            for eff in self.fEfficiency1D:
                rlist.Add(eff)
            for eff in self.fFoldedEfficiency1D:
                rlist.Add(eff)
            for h in self.fTruth1D:
                rlist.Add(h)
            for h in self.fMeasured1D:
                rlist.Add(h)
            for h in self.fReconstructedTruth1D:
                rlist.Add(h)
            for h in self.fFoldedTruth1D:
                rlist.Add(h)
            for h in self.fResponseMatrix1D:
                rlist.Add(h)
        return rlist

    def GenerateFoldedTruth(self):
        self.fFoldedTruth = self.FoldResponse(self.fTruth)
        self.fFoldedTruth.SetName("{0}_FoldedTruth".format(self.fName))
        self.fFoldedTruth.SetTitle("{0}_FoldedTruth".format(self.fName))
        self.fFoldedReconstructedTruth = self.FoldResponse(self.fReconstructedTruth)
        self.fFoldedReconstructedTruth.SetName("{0}_FoldedReconstructedTruth".format(self.fName))
        self.fFoldedReconstructedTruth.SetTitle("{0}_FoldedReconstructedTruth".format(self.fName))
        self.fFoldedTruth1D = []
        if len(self.fAxis) == 2:
            binLimits = BinLimits()
            for bin, (binMin,binMax) in zip(range(0, self.fFoldedTruth.GetNbinsX()+1), zip([1]+range(1, self.fFoldedTruth.GetNbinsX()+1), [self.fFoldedTruth.GetNbinsX()]+range(1, self.fFoldedTruth.GetNbinsX()+1))):
                binLimits.SetFromAxis(self.fAxis[0].fDetectorAxis, bin)
                hist = self.fFoldedTruth.ProjectionY("{0}_FoldedTruth_{1}".format(self.fName, binLimits.GetName()), binMin, binMax)
                hist.SetTitle(binLimits.GetTitle())
                self.fFoldedTruth1D.append(hist)

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

        if self.fFoldedTruth:
            #if len(self.fAxis) == 1:
            #self.fFoldedEfficiency = ROOT.TGraphAsymmErrors(self.fMeasured, self.fFoldedTruth, "b(1,1) mode")
            #self.fFoldedEfficiency.SetName("{0}_FoldedEfficiency".format(self.fName))
            #self.fFoldedEfficiency.SetTitle("{0}_FoldedEfficiency".format(self.fName))
            #else:
            self.fFoldedEfficiency = self.GenerateMeasured(self.fAxis, "FoldedEfficiency")
            self.fFoldedEfficiency.Divide(self.fMeasured, self.fFoldedTruth)
            self.fFoldedEfficiency.GetXaxis().SetTitle(self.fAxis[0].fDetectorAxis.GetTitle())
            self.fFoldedEfficiency.GetYaxis().SetTitle("FoldedEfficiency")

        if len(self.fAxis) == 2:
            self.GeneratePartialMultiEfficiency()

    def IterateResponseMatrix(self, coord=None, axis=-1, minAxis=0, maxAxis=-1):
        if maxAxis == -1:
            maxAxis = len(self.fAxis)-1
        if coord is None:
            coord = array('i',[-1]*(len(self.fAxis)*2))
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
        for bin,(truth,recoTruth) in enumerate(zip(self.fTruth1D, self.fReconstructedTruth1D)):
            eff = self.GeneratePartialMultiEfficiencyForBin(self.fAxis[1], self.fAxis[0].fTruthAxis, bin, truth, recoTruth)
            self.fEfficiency1D.append(eff)

        for bin,(foldedTruth,measured) in enumerate(zip(self.fFoldedTruth1D, self.fMeasured1D)):
            eff = self.GeneratePartialMultiEfficiencyForBin(self.fAxis[1], self.fAxis[0].fDetectorAxis, bin, foldedTruth, measured)
            self.fFoldedEfficiency1D.append(eff)

    def GeneratePartialMultiEfficiencyForBin(self, axis, axisProj, bin, truth, recoTruth):
        hname = truth.GetName().replace("Truth", "Efficiency")

        #eff = ROOT.TGraphAsymmErrors(recoTruthProj, truthProj, "b(1,1) mode")
        #eff.SetName("{0}_{1}".format(binName, name))

        if "Folded" in hname:
            eff = self.GenerateMeasured([axis], hname)
        else:
            eff = self.GenerateTruth([axis], hname)

        eff.SetTitle(truth.GetTitle())
        eff.GetXaxis().SetTitle(axisProj.GetTitle())
        eff.GetYaxis().SetTitle("Efficiency")
        eff.Divide(recoTruth, truth)

        return eff

    def GenerateResponseMatrix(self, axis):
        hname = "{0}_DetectorResponse".format(self.fName)
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
            hist = ROOT.TH1D(hname, hname, len(axis[0].fBins)-1, array('d',axis[0].fBins))
            hist.GetXaxis().SetTitle(axis[0].GetTitle())
            hist.GetYaxis().SetTitle("counts")
        elif len(axis) == 2:
            hist = ROOT.TH2D(hname, hname, len(axis[0].fBins)-1, array('d',axis[0].fBins), len(axis[1].fBins)-1, array('d',axis[1].fBins))
            hist.GetXaxis().SetTitle(axis[0].GetTitle())
            hist.GetYaxis().SetTitle(axis[1].GetTitle())
            hist.GetZaxis().SetTitle("counts")
        elif len(axis) == 3:
            hist = ROOT.TH3D(hname, hname, len(axis[0].fBins)-1, array('d',axis[0].fBins), len(axis[1].fBins)-1, array('d',axis[1].fBins), len(axis[2].fBins)-1, array('d',axis[2].fBins))
            hist.GetXaxis().SetTitle(axis[0].GetTitle())
            hist.GetYaxis().SetTitle(axis[1].GetTitle())
            hist.GetZaxis().SetTitle(axis[2].GetTitle())
        else:
            nbins = array('i', [len(a.fBins) for a in axis])
            hist = ROOT.THnSparseD(hname, hname, len(axis), nbins)
            for i,a in enumerate(axis):
                hist.GetAxis(i).Set(len(a.fBins)-1, array('d',a.fBins))
                hist.GetAxis(i).SetTitle(axis[i].GetTitle())
                
        hist.Sumw2()
        
        return hist

    def FillResponseMatrix(self, axis, resp, recoDmeson, truthDmeson, recoJet, truthJet, w):
        naxis = len(axis)
        values = array('d', [0]*(naxis*2))
        for i,a in enumerate(axis):
            if a == "jet_pt":
                values[i] = recoJet.fPt
                values[i+naxis] = truthJet.fPt
            elif a == "d_pt":
                values[i] = recoDmeson.fPt
                values[i+naxis] = truthDmeson.fPt
            elif a == "d_z":
                values[i] = recoJet.fZ
                values[i+naxis] = truthJet.fZ

        self.FillHistogram(resp, values, w)

    def FillSpectrum(self, axis, hist, dmeson, jet, w):
        values = array('d')
        for a in axis:
            if a == "jet_pt":
                values.append(jet.fPt)
            elif a == "d_pt":
                values.append(dmeson.fPt)
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

        if self.fResponseMatrix1D and recoJet and truthJet:
            self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[self.fResponseMatrix.GetAxis(0).GetNbins()+1], recoDmeson, truthDmeson, recoJet, truthJet, w)
            bin = self.fResponseMatrix.GetAxis(0).FindBin(recoJet.fPt)
            if bin >= 1 and bin <= self.fResponseMatrix.GetAxis(0).GetNbins():
                self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[bin], recoDmeson, truthDmeson, recoJet, truthJet, w)
                self.FillResponseMatrix([self.fAxis[1].fDetectorAxis.fName], self.fResponseMatrix1D[0], recoDmeson, truthDmeson, recoJet, truthJet, w)

    def FillMeasured(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fDetectorAxis.fName for a in self.fAxis], self.fMeasured, dmeson, jet, w)

        if self.fMeasured1D and jet:
            self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[self.fMeasured.GetNbinsX()+1], dmeson, jet, w)
            bin = self.fMeasured.GetXaxis().FindBin(jet.fPt)
            if bin >= 1 and bin <= self.fMeasured.GetNbinsX():
                self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[0], dmeson, jet, w)
                self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fMeasured1D[bin], dmeson, jet, w)

    def FillTruth(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fTruthAxis.fName for a in self.fAxis], self.fTruth, dmeson, jet, w)

        if self.fTruth1D and jet:
            self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[self.fTruth.GetNbinsX()+1], dmeson, jet, w)
            bin = self.fTruth.GetXaxis().FindBin(jet.fPt)
            if bin >= 1 and bin <= self.fTruth.GetNbinsX():
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[bin], dmeson, jet, w)
                self.FillSpectrum([self.fAxis[1].fTruthAxis.fName], self.fTruth1D[0], dmeson, jet, w)

    def FillRecoTruth(self, dmeson, jet, w):
        if not jet or jet.fPt > 0:
            self.FillSpectrum([a.fDetectorAxis.fName for a in self.fAxis], self.fReconstructedTruth, dmeson, jet, w)

        if self.fReconstructedTruth1D and jet:
            self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fReconstructedTruth1D[self.fReconstructedTruth.GetNbinsX()+1], dmeson, jet, w)
            bin = self.fReconstructedTruth.GetXaxis().FindBin(jet.fPt)
            if bin >= 1 and bin <= self.fReconstructedTruth.GetNbinsX():
                self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fReconstructedTruth1D[bin], dmeson, jet, w)
                self.FillSpectrum([self.fAxis[1].fDetectorAxis.fName], self.fReconstructedTruth1D[0], dmeson, jet, w)

    def Fill(self, dmeson, w):
        if self.fJetInfo:
            jetTruth = getattr(dmeson, "{0}_truth".format(self.fJetName))
            jetMeasured = getattr(dmeson, "{0}_reco".format(self.fJetName))
        else:
            jetTruth = None
            jetMeasured = None

        dMesonTruth = dmeson.DmesonJet.fGenerated
        dMesonMeasured = dmeson.DmesonJet.fReconstructed

        weff = self.fWeightEfficiency.GetEfficiencyWeight(dmeson, self.fJetName)

        if dMesonTruth.fPt > 0 and dMesonMeasured.fPt > 0:
            self.FillDetectorResponse(dMesonMeasured, dMesonTruth, jetMeasured, jetTruth, w*weff)
            self.FillMeasured(dMesonMeasured, jetMeasured, w*weff)
            self.FillRecoTruth(dMesonTruth, jetTruth, w*weff)

        if dMesonTruth.fPt > 0:
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

    def GetTitle(self, label = ""):
        if not label:
            label = self.fLabel

        if label == "nolabel":
            label = ""

        if self.fName == "jet_pt":
            if label:
                title = "#it{{p}}_{{T,jet}}^{{{0}}} (GeV/#it{{c}})".format(label)
            else:
                title = "#it{p}_{T,jet} (GeV/#it{c})"
        if self.fName == "d_pt":
            if label:
                title = "#it{{p}}_{{T,D}}^{{{0}}} (GeV/#it{{c}})".format(label)
            else:
                title = "#it{p}_{T,D} (GeV/#it{c})"
        if self.fName == "d_z":
            if label:
                title = "#it{{z}}_{{||,D}}^{{{0}}}".format(label)
            else:
                title = "#it{z}_{||,D}"

        return title

class Spectrum:
    def __init__(self, config, name):
        self.fName = name
        self.fBins = config["bins"]
        self.fAxis = []
        if config["jet_pt"]:
            self.fAxis.append(Axis("jet_pt", config["jet_pt"]))
        if config["d_pt"]:
            self.fAxis.append(Axis("d_pt", config["d_pt"]))
        if config["d_z"]:
            self.fAxis.append(Axis("d_z", config["d_z"]))
        self.fHistogram = None

class BinSet:
    def __init__(self):
        self.fBins = dict()
        
    def AddBins(self, name, jetPtLimits = [0, -1], ZLimits = [0, -1], DPtLimits = [0, -1], _showJetPt = True, _showDPt = True, _showDZ = True):
        self.fBins[name] = []
        for _DPtMin,_DPtMax in zip(DPtLimits[:-1],DPtLimits[1:]):
            for _jetPtMin,_jetPtMax in zip(jetPtLimits[:-1],jetPtLimits[1:]):
                for _DZMin,_DZMax in zip(ZLimits[:-1],ZLimits[1:]):
                    self.fBins[name].append(BinLimits(jetPtMin=_jetPtMin, jetPtMax=_jetPtMax, DPtMin = _DPtMin, DPtMax = _DPtMax, DZMin = _DZMin, DZMax = _DZMax, 
                                                      showJetPt = _showJetPt, showDPt = _showDPt, showDZ = _showDZ))

    def FindBin(self, dmeson, jetDef):
        for bins in self.fBins.itervalues():
            for bin in bins:
                if bin.IsInBinLimits(dmeson, jetDef):
                    yield bin

class BinLimits:
    def __init__(self, jetPtMin = 0, jetPtMax = -1, DPtMin = 0, DPtMax = -1, DZMin = 0, DZMax = -1, showJetPt = True, showDPt = True, showDZ = True):
        self.fJetPtMin = jetPtMin
        self.fJetPtMax = jetPtMax
        self.fDPtMin = DPtMin
        self.fDPtMax = DPtMax
        self.fDZMin = DZMin
        self.fDZMax = DZMax
        self.fShowJetPt = showJetPt
        self.fShowDPt = showDPt
        self.fShowDZ = showDZ
        self.fInvMassHisto = None
        self.fMassFitter = None
        
    def SetFromAxis(self, axis, binIndex):
        self.fShowJetPt = False
        self.fShowDPt = False
        self.fShowDZ = False

        if axis.fName == "jet_pt":
            self.fShowJetPt = True
            if binIndex == 0:
                self.fJetPtMin = axis.fBins[0]
                self.fJetPtMax = axis.fBins[-1]
            else:
                self.fJetPtMin = axis.fBins[binIndex-1]
                self.fJetPtMax = axis.fBins[binIndex]                

        if axis.fName == "d_pt":
            self.fShowDPt = True
            if binIndex == 0:
                self.fDPtMin = axis.fBins[0]
                self.fDPtMax = axis.fBins[-1]
            else:
                self.fDPtMin = axis.fBins[binIndex-1]
                self.fDPtMax = axis.fBins[binIndex]

        if axis.fName == "d_z":
            self.fShowDZ = True
            if binIndex == 0:
                self.fDZMin = axis.fBins[0]
                self.fDZMax = axis.fBins[-1]
            else:
                self.fDZMin = axis.fBins[binIndex-1]
                self.fDZMax = axis.fBins[binIndex]
      
    def SetMassFitter(self, fitter):
        self.fMassFitter = fitter
    
    def SetJetPtLimits(self, min, max):
        self.fJetPtMin = min
        self.fJetPtMax = max
        
    def SetDPtLimits(self, min, max):
        self.fDPtMin = min
        self.fDPtMax = max
        
    def SetDZLimits(self, min, max):
        self.fDZMin = min
        self.fDZMax = max
        
    def IsInBinLimits(self, dmeson, jetDef):
        if self.fDPtMax > self.fDPtMin and (dmeson.DmesonJet.fPt < self.fDPtMin or dmeson.DmesonJet.fPt > self.fDPtMax):
            return False
        
        jetName = "Jet_AKT{0}{1}_pt_scheme".format(jetDef["type"], jetDef["radius"])
        jet = getattr(dmeson, jetName)
        
        jetPt = jet.fPt
        if self.fJetPtMax > self.fJetPtMin and (jetPt < self.fJetPtMin or jetPt > self.fJetPtMax):
            return False
        
        DZ = jet.fZ
        if self.fDZMax > self.fDZMin and (DZ < self.fDZMin or DZ > self.fDZMax):
            return False
        
        return True
    
    def GetBinCenter(self, axis):
        if axis == "jet_pt":
            if self.fJetPtMax > self.fJetPtMin:
                return (self.fJetPtMax + self.fJetPtMin) / 2
            else:
                return -1
            
        if axis == "d_pt":
            if self.fDPtMax > self.fDPtMin:
                return (self.fDPtMax + self.fDPtMin) / 2
            else:
                return -1
            
        if axis == "d_z":
            if self.fDZMax > self.fDZMin:
                return (self.fDZMax + self.fDZMin) / 2
            else:
                return -1
    
    def GetName(self):
        name = ""
        if self.fDPtMax > self.fDPtMin:
            name += "DPt_{0}_{1}_".format(int(self.fDPtMin*100), int(self.fDPtMax*100))
        
        if self.fJetPtMax > self.fJetPtMin:
            name += "JetPt_{0}_{1}_".format(int(self.fJetPtMin*100), int(self.fJetPtMax*100))
            
        if self.fDZMax > self.fDZMin:
            name += "DZ_{0}_{1}_".format(int(self.fDZMin*100), int(self.fDZMax*100))
        
        #remove last "_"
        if name:
            name = name[:-1]
        return name
        
    def GetTitle(self):
        title = ""
        if self.fDPtMax > self.fDPtMin and self.fShowDPt:
            title += "{0:.1f} < #it{{p}}_{{T,D}} < {1:.1f} GeV/#it{{c}}, ".format(self.fDPtMin, self.fDPtMax)
        
        if self.fJetPtMax > self.fJetPtMin and self.fShowJetPt:
            title += "{0:.0f} < #it{{p}}_{{T,jet}} < {1:.0f} GeV/#it{{c}}, ".format(self.fJetPtMin, self.fJetPtMax)
            
        if self.fDZMax > self.fDZMin and self.fShowDZ:
            title += "{0:.1f} < #it{{z}}_{{||, D}} < {1:.1f}, ".format(self.fDZMin, self.fDZMax)
        
        #remove last ", "
        if title:
            title = title[:-2]
        return title
    
    def Print(self):
        print(self.GetTitle())
    
    def CreateInvMassHisto(self, trigger, DMesonDef, xAxis, yAxis, nMassBins, minMass, maxMass):
        hname = "InvMass_{0}_{1}_{2}".format(trigger, DMesonDef, self.GetName())
        htitle = "{0} - {1} Invariant Mass: {2};{3};{4}".format(trigger, DMesonDef, self.GetTitle(), xAxis, yAxis)
        self.fInvMassHisto = ROOT.TH1D(hname, htitle, nMassBins, minMass-(maxMass-minMass)/2, maxMass+(maxMass-minMass)/2)
        self.fInvMassHisto.Sumw2()
        self.fInvMassHisto.SetMarkerSize(0.9)
        self.fInvMassHisto.SetMarkerStyle(ROOT.kFullCircle)
        self.fInvMassHisto.SetMarkerColor(ROOT.kBlue+2)
        self.fInvMassHisto.SetLineColor(ROOT.kBlue+2)