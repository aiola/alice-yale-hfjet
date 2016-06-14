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
        self.fBins = []
        
    def AddBins(self, jetPtLimits = [0, -1], ZLimits = [0, -1], DPtLimits = [0, -1]):
        for _DPtMin,_DPtMax in zip(DPtLimits[:-1],DPtLimits[1:]):
            for _jetPtMin,_jetPtMax in zip(jetPtLimits[:-1],jetPtLimits[1:]):
                for _DZMin,_DZMax in zip(ZLimits[:-1],ZLimits[1:]):
                    self.fBins.append(BinLimits(jetPtMin=_jetPtMin, jetPtMax=_jetPtMax, DPtMin = _DPtMin, DPtMax = _DPtMax, DZMin = _DZMin, DZMax = _DZMax))

    def FindBin(self, dmeson, jetDef):
        for bin in self.fBins:
            if bin.IsInBinLimits(dmeson, jetDef):
                yield bin

class BinLimits:
    def __init__(self, jetPtMin = 0, jetPtMax = -1, DPtMin = 0, DPtMax = -1, DZMin = 0, DZMax = -1):
        self.fJetPtMin = jetPtMin
        self.fJetPtMax = jetPtMax
        self.fDPtMin = DPtMin
        self.fDPtMax = DPtMax
        self.fDZMin = DZMin
        self.fDZMax = DZMax
        self.fInvMassHisto = None
    
    def SetJetPtLimits(self, min, max):
        self.fJetPtMin = min
        self.fJetPtMax = max
        
    def SetDPtLimits(self, min, max):
        self.fDPtMin = min
        self.fDPtMax = max
        
    def SetDZLimits(self, min, max):
        self.fDZMin = min
        self.fDZMax = max
        
    def IsInBinLimits(dmeson, jetDef):
        if self.fDPtMax > self.fDPtMin and (dmeson.fPt < self.fDPtMin or dmeson.fPt > self.fDPtMax):
            return false
        
        jetPt = getattr(dmeson, jetDef.GetPtLeafName())
        if self.fJetPtMax > self.fJetPtMin and (jetPt < self.JetDPtMin or jetPt > self.fJetPtMax):
            return false
        
        DZ = getattr(dmeson, jetDef.GetDZLeafName())
        if self.fDZMax > self.fDZMin and (DZ < self.fDZMin or DZ > self.fDZMax):
            return false
        
        return true
    
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
        if self.fDPtMax > self.fDPtMin:
            title += "{0:.2f} < #it{p}_{T,D} < {0:.2f} GeV/#it{c}, ".format(self.fDPtMin, self.fDPtMax)
        
        if self.fJetPtMax > self.fJetPtMin:
            title += "{0:.2f} < #it{p}_{T,jet} < {0:.2f} GeV/#it{c}, ".format(self.fJetPtMin, self.fJetPtMax)
            
        if self.fDZMax > self.fDZMin:
            title += "{0:.2f} < #it{z}_{||, D} < {0:.2f}, ".format(self.fDZMin, self.fDZMax)
        
        #remove last ", "
        if title:
            title = title[:-2]
        return title
    
    def CreateInvMassHisto(self, trigger, DMesonDef, nMassBins, minMass, maxMass):
        hname = "InvMass_{0}_{1}_{2}".format(trigger, DMesonDef, self.GetName())
        htitle = "{0} - {1} Invariant Mass: {2}".format(trigger, DMesonDef, self.GetTitle())
        self.fInvMassHisto = ROOT.TH1D(hname, htitle, nMassBins, minMass, maxMass)
        globalList.append(self.fInvMassHisto)
        
class DMesonJetDataProjector:
    def __init__(self, inputPath, train, fileName, taskName):
        self.fInputPath = inputPath
        self.fTrain = train
        self.fFileName = fileName
        self.fTaskName = taskName
        self.fChains = dict()
        
    def GenerateChain(self, treeName):
        path = "{0}/{1}".format(self.fInputPath, self.fTrain)
        
        files = find_file(path, self.fFileName)
        
        print(list(files))
        
        chain = ROOT.TChain(treeName)
        
        self.fChains[treeName] = chain

        for file in files:
            print("Adding file {0}...".format(file))
            chain.Add(file)
        
        chain.Print()
        return chain
        
        
    def GetInvMassHisograms(self, trigger, DMesonDef, jetDefinitions, binSet, nMassBins, minMass, maxMass):
        treeName = "{0}_{1}_{2}".format(self.fTaskName, trigger, DMesonDef)
        
        chain = self.GenerateChain(treeName)
        
        print("Running analysis on tree {0}. Total number of entries is {1}".format(treeName, chain.GetEntries()))
        
        for dmeson in chain:
            for jetDef in jetDefinitions:
                bins = binSet.FindBin(dmeson, jetDef)
                for bin in bins:
                    if not bin.fInvMassHisto:
                        bin.CreateInvMassHisto(trigger, DMesonDef, nMassBins, minMass, maxMass)
                    bin.fInvMassHisto.Fill(dmeson.fInvMass)
        