#!/usr/bin/env python
#python base classes and utilities for D Meson jet analysis

import ROOT
from array import array
import os

def find_file(path, file_name):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file == file_name:
                yield os.path.join(root, file)
                
class DetectorResponse:
    def __init__(self, name, axis):
        self.fAxis = axis
        self.fName = name
        self.fResponseMatrix = self.GenerateResponseMatrix(axis)
        self.fEfficiency = None
        self.fTruth = self.GenerateTruth(axis)
        self.fMeasured = self.GenerateMeasured(axis)
        self.fMissedTruth = self.GenerateTruth(axis, "MissedTruth")
        if len(self.fAxis) > 1:
            self.fEfficiency1D = dict()
            for axis in self.fAxis:
                self.fEfficiency1D[axis.fTruthAxis.fName] = None
            
    def Generate1DEfficiency(self):
        recoTruth = self.fResponseMatrix.ProjectionY()
        self.fEfficiency = ROOT.TGraphAsymmErrors(recoTruth, self.fTruth, "b(1,1) mode")
        self.fEfficiency.SetName("{0}_Efficiency".format(self.fName))
        self.fEfficiency.GetXaxis().SetTitle(self.fAxis[0].fTruthAxis.GetTitle())
        self.fEfficiency.GetYaxis().SetTitle("Efficiency")

    def GenerateMultiEfficiency(self):
        if  len(self.fAxis) == 2:
            self.fEfficiency = self.fResponseMatrix.Projection(3, 1)
            self.fEfficiency.GetXaxis().SetTitle(self.fAxis[0].fTruthAxis.GetTitle())
            self.fEfficiency.GetYaxis().SetTitle(self.fAxis[1].fTruthAxis.GetTitle())
            self.fEfficiency.GetZaxis().SetTitle("Efficiency")
        elif len(self.fAxis) == 3:
            self.fEfficiency = self.fResponseMatrix.Projection(5, 3, 1)
            self.fEfficiency.GetXaxis().SetTitle(self.fAxis[0].fTruthAxis.GetTitle())
            self.fEfficiency.GetYaxis().SetTitle(self.fAxis[1].fTruthAxis.GetTitle())
            self.fEfficiency.GetZaxis().SetTitle(self.fAxis[2].fTruthAxis.GetTitle())
        else:
            dim = array.array('i')
            for i in reversed(range(0, len(self.fAxis)*2, 2)):
                dim.append(i+1)
            self.fEfficiency = self.fResponseMatrix.Projection(len(dim), dim)
        self.fEfficiency.Divide(self.fTruth)
        self.fEfficiency.SetName("{0}_Efficiency".format(self.fName))
        
        for i,a in enumerate(self.fAxis):
            if len(self.fAxis) < 4:
                if i == 0:
                    truth = self.fTruth.ProjectionX()
                elif i == 1:
                    truth = self.fTruth.ProjectionY()
                elif i == 2:
                    truth = self.fTruth.ProjectionZ()
            else:
                truth = self.fTruth.Projection(i)
                
            recoTruth = self.fResponseMatrix.Projection(i*2+1)
            self.fEfficiency1D[a.fTruthAxis.fName] = ROOT.TGraphAsymmErrors(recoTruth, truth, "b(1,1) mode")
            self.fEfficiency1D[a.fTruthAxis.fName].SetName("{0}_{1}_Efficiency".format(self.fName, a.fTruthAxis.fName))
            self.fEfficiency1D[a.fTruthAxis.fName].GetXaxis().SetTitle(a.fTruthAxis.GetTitle())
            self.fEfficiency1D[a.fTruthAxis.fName].GetYaxis().SetTitle("Efficiency")

    def GenerateEfficiency(self):
        if len(self.fAxis) == 1:
            self.Generate1DEfficiency()
        else:
            self.GenerateMultiEfficiency()

    def GenerateResponseMatrix(self, axis):
        hname = "{0}_DetectorResponse".format(self.fName)
        haxis = []
        for a in axis:
            haxis.append(a.fDetectorAxis)
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

    def FillResponseMatrix(self, resp, recoDmeson, truthDmeson, recoJet, truthJet, w):
        values = array('d')
        for a in self.fAxis:
            if a.fDetectorAxis.fName == "jet_pt":
                values.append(recoJet.fPt)
                values.append(truthJet.fPt)
            elif a.fDetectorAxis.fName == "d_pt":
                values.append(recoDmeson.fPt)
                values.append(truthDmeson.fPt)
            elif a.fDetectorAxis.fName == "d_z":
                values.append(recoJet.fZ)
                values.append(truthJet.fZ)

        self.FillHistogram(resp, values, w)

    def FillSpectrum(self, hist, dmeson, jet, w):
        values = array('d')
        for a in self.fAxis:
            if a.fDetectorAxis.fName == "jet_pt":
                values.append(jet.fPt)
            elif a.fDetectorAxis.fName == "d_pt":
                values.append(dmeson.fPt)
            elif a.fDetectorAxis.fName == "d_z":
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

    def Fill(self, dmeson, jetName, w):
        jetTruth = getattr(dmeson, "{0}_truth".format(jetName))
        jetMeasured = getattr(dmeson, "{0}_reco".format(jetName))

        if jetTruth.fPt > 0 and jetMeasured.fPt > 0:
            self.FillResponseMatrix(self.fResponseMatrix, dmeson.DmesonJet.fReconstructed, dmeson.DmesonJet.fGenerated, jetMeasured, jetTruth, w)
        
        if jetTruth.fPt > 0 and not jetMeasured.fPt > 0:
            self.FillSpectrum(self.fMissedTruth, dmeson.DmesonJet.fGenerated, jetTruth, w)
        
        if jetTruth.fPt > 0:
            self.FillSpectrum(self.fTruth, dmeson.DmesonJet.fGenerated, jetTruth, w)

        if jetMeasured.fPt > 0:
            self.FillSpectrum(self.fMeasured, dmeson.DmesonJet.fReconstructed, jetMeasured, w)

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