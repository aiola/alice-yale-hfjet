#!/usr/bin/env python
#python program to project D meson jet trees into invariant mass histograms

import ROOT
import math
import os

def find_file(path, file_name):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file == file_name:
                yield os.path.join(root, file)

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
        
class DMesonJetDataProjector:
    def __init__(self, inputPath, train, fileName, taskName, maxEvents):
        self.fInputPath = inputPath
        self.fTrain = train
        self.fFileName = fileName
        self.fTaskName = taskName
        self.fMaxEvents = maxEvents
        self.fMassAxisTitle = "#it{m}(K#pi) GeV/#it{c}^{2}"
        self.fYieldAxisTitle = "counts"
        self.fChains = dict()
        
    def GenerateChain(self, treeName):
        path = "{0}/{1}".format(self.fInputPath, self.fTrain)
        
        files = find_file(path, self.fFileName)
        
        chain = ROOT.TChain(treeName)
        
        self.fChains[treeName] = chain

        for file in files:
            print("Adding file {0}...".format(file))
            chain.Add(file)
        
        #chain.Print()
        return chain
        
    def GetInvMassHisograms(self, trigger, DMesonDef, jetDefinitions, binSet, nMassBins, minMass, maxMass):
        if trigger:
            treeName = "{0}_{1}_{2}".format(self.fTaskName, trigger, DMesonDef)
        else:
            treeName = "{0}_{1}".format(self.fTaskName, DMesonDef)
        
        chain = self.GenerateChain(treeName)
        
        print("Running analysis on tree {0}. Total number of entries is {1}".format(treeName, chain.GetEntries()))
        if self.fMaxEvents > 0:
            print("The analysis will stop at the {0} entry.".format(self.fMaxEvents))
        
        for i,dmeson in enumerate(chain):
            if i % 10000 == 0:
                print("D meson candidate n. {0}".format(i))
                if self.fMaxEvents > 0 and i > self.fMaxEvents:
                    print("Stopping the analysis.")
                    break
            for jetDef in jetDefinitions:
                bins = binSet.FindBin(dmeson, jetDef)
                
                for bin in bins:
                    if not bin.fInvMassHisto:
                        bin.CreateInvMassHisto(trigger, DMesonDef, self.fMassAxisTitle, self.fYieldAxisTitle, nMassBins, minMass, maxMass)
                    bin.fInvMassHisto.Fill(dmeson.DmesonJet.fInvMass)
        