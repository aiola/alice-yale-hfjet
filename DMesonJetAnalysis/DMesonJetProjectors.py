#!/usr/bin/env python
#python program to project D meson jet trees into invariant mass histograms

import ROOT
import math
import os
from DMesonJetBase import *
import array
from bisect import bisect

class SimpleWeight:
    def GetEfficiencyWeight(self, dmeson):
        return 1.

class EfficiencyWeightCalculator:
    def __init__(self, filename="", listname="", objectname=""):
        self.fBreakpoints = []
        self.fEfficiencyValues = []
        if filename and listname and objectname:
            self.LoadEfficiency(filename,listname,objectname)

    def GetObjectFromRootFile(self, file, listname, objectname):
        rlist = file.Get(listname)
        if not rlist:
            print("Could not find list '{0}' in file '{1}'".format(listname, file.GetName()))
            return None
        obj = rlist.FindObject(objectname)
        if not obj:
            print("Could not find object '{0}' in list '{1}'".format(objectname, listname))
            return None
        return obj

    def LoadEfficiency(self, filename, listname, objectname):
        self.fBreakpoints = []
        self.fEfficiencyValues = []
        file = ROOT.TFile.Open(filename)
        if not file or file.IsZombie():
            print("Could not open file '{0}'! Effieciency could not be loaded!".format(filename))
            return
        obj = self.GetObjectFromRootFile(file, listname, objectname)
        if not obj:
            return
        self.fRootObject = obj.Clone()
        file.Close()
        print("Efficiency weights loaded from object '{0}/{1}' in file '{2}'".format(listname, objectname, filename))
#         if isinstance(self.fRootObject, ROOT.TH1):
#             for ibin in range(1, self.fRootObject.GetNbinsX()):
#                 self.fBreakpoints.append(self.fRootObject.GetXaxis().GetBinUpEdge(ibin))
#                 self.fEfficiencyValues.append(1. / self.fRootObject.GetBinContent(ibin))
# 
#             self.fEfficiencyValues.append(self.fRootObject.GetBinContent(self.fRootObject.GetNbinsX()))
# 
#         elif isinstance(self.fRootObject, ROOT.TGraph):
#             for ipoint in range(0, self.fRootObject.GetN()-1):
#                 self.fBreakpoints.append(self.fRootObject.GetX()[ipoint] + self.fRootObject.GetErrorX(ipoint))
#                 self.fEfficiencyValues.append(1. / self.fRootObject.GetY()[ipoint])
# 
#             self.fEfficiencyValues.append(1. / self.fRootObject.GetY()[self.fRootObject.GetN()-1])
# 
#         print(len(self.fBreakpoints), self.fBreakpoints)
#         print(len(self.fEfficiencyValues), self.fEfficiencyValues)

    def GetEfficiencyWeight(self, dmeson):
#         i = bisect(self.fBreakpoints, dmeson.fPt)
#         print("Efficiency weight for pt {0} is at position {1}".format(dmeson.fPt, i))
#         print("Efficiency is {0}".format(self.fEfficiencyValues[i]))
#         return self.fEfficiencyValues[i]
        if isinstance(self.fRootObject, ROOT.TH1):
            return 1. / self.fRootObject.Interpolate(dmeson.fPt)

        elif isinstance(self.fRootObject, ROOT.TGraph):
            return 1. / self.fRootObject.Eval(dmeson.fPt)

class DMesonJetDataProjector:
    def __init__(self, inputPath, train, fileName, taskName, maxEvents, effWeight=SimpleWeight()):
        self.fInputPath = inputPath
        self.fTrain = train
        self.fFileName = fileName
        self.fTaskName = taskName
        self.fMaxEvents = maxEvents
        self.fMassAxisTitle = "#it{m}(K#pi) GeV/#it{c}^{2}"
        self.fYieldAxisTitle = "counts"
        self.fChain = None
        self.fEfficiencyWeight = effWeight

    def GenerateChain(self, treeName):
        path = "{0}/{1}".format(self.fInputPath, self.fTrain)
        
        files = find_file(path, self.fFileName)
        
        self.fChain = ROOT.TChain(treeName)

        for file in files:
            print("Adding file {0}...".format(file))
            self.fChain.Add(file)

        #self.fChain.Print()

    def GetInvMassHisograms(self, trigger, DMesonDef, jetDefinitions, binSet, nMassBins, minMass, maxMass):
        if trigger:
            treeName = "{0}_{1}_{2}".format(self.fTaskName, trigger, DMesonDef)
        else:
            treeName = "{0}_{1}".format(self.fTaskName, DMesonDef)

        self.GenerateChain(treeName)

        print("Running analysis on tree {0}. Total number of entries is {1}".format(treeName, self.fChain.GetEntries()))
        if self.fMaxEvents > 0:
            print("The analysis will stop at the {0} entry.".format(self.fMaxEvents))
        
        for i,dmeson in enumerate(self.fChain):
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
                    bin.fInvMassHisto.Fill(dmeson.DmesonJet.fInvMass, self.fEfficiencyWeight.GetEfficiencyWeight(dmeson.DmesonJet))

class DMesonJetResponseProjector:
    def __init__(self, inputPath, train, fileName, taskName, maxEvents, effWeight=SimpleWeight()):
        self.fInputPath = inputPath
        self.fTrain = train
        self.fFileName = fileName
        self.fTaskName = taskName
        self.fMaxEvents = maxEvents
        self.fChain = None
        self.fPtHardBin = -1
        self.fPeriod = ""
        self.fWeight = 1
        self.fEfficiencyWeight = effWeight

    def GenerateChain(self, treeName):
        path = "{0}/{1}".format(self.fInputPath, self.fTrain)

        files = find_file(path, self.fFileName)

        self.fChain = ROOT.TChain(treeName)

        for file in files:
            print("Adding file {0}...".format(file))
            self.fChain.Add(file)

        #self.fChain.Print()

    def ExtractCurrentFileInfo(self):
        fname = self.fChain.GetCurrentFile().GetName()
        lastSlash = fname.rfind('/')
        secondLastSlash = fname.rfind('/',0,lastSlash-1)
        thirdLastSlash = fname.rfind('/',0,secondLastSlash-1)
        self.fPtHardBin = int(fname[secondLastSlash+1:lastSlash])
        self.fPeriod = fname[thirdLastSlash+1:secondLastSlash]
        
    def ExtractWeightFromHistogramList(self, hlist):
        xsection = hlist.FindObject("fHistXsection")
        trials  = hlist.FindObject("fHistTrials")
        
        if not trials or not xsection:
            print("Could not find trail and x-section information!")
            hlist.Print() 
            self.fWeight = 1
            return

        valNTRIALS = trials.GetBinContent(self.fPtHardBin+1);
        valXSEC = xsection.GetBinContent(self.fPtHardBin+1);
        scalingFactor = 0;
        if valNTRIALS > 0:
            self.fWeight = valXSEC/valNTRIALS;

    def RecalculateWeight(self):
        oldPtHardBin = self.fPtHardBin
        oldPeriod = self.fPeriod

        self.ExtractCurrentFileInfo()

        if oldPtHardBin == self.fPtHardBin and oldPeriod == self.fPeriod:
            return

        listName = "{0}_histos".format(self.fTaskName)
        hlist = self.fChain.GetCurrentFile().Get(listName)

        if not hlist:
            print("Could not get list '{0}' from file '{1}'".format(listName, self.fChain.GetCurrentFile().GetName()))
            return 1

        self.ExtractWeightFromHistogramList(hlist)
        
        print("Period: {0}\nPt hard bin: {1}\nWeight: {2}".format(self.fPeriod, self.fPtHardBin, self.fWeight))

    def GetDetectorResponse(self, respDefinitions, DMesonDef, jetDefinitions):
        response = dict()
        for jetDef in jetDefinitions:
            jetName = "Jet_AKT{0}{1}_pt_scheme".format(jetDef["type"], jetDef["radius"])
            for axisName,axisDef in respDefinitions.iteritems():
                respName = "{0}_{1}_{2}".format(DMesonDef, jetName, axisName)
                response[respName] = DetectorResponse(respName, jetName, axisDef)
                treeName = "{0}_{1}".format(self.fTaskName, DMesonDef)

        self.GenerateChain(treeName)

        print("Running analysis on tree {0}. Total number of entries is {1}".format(treeName, self.fChain.GetEntries()))
        if self.fMaxEvents > 0:
            print("The analysis will stop at the {0} entry.".format(self.fMaxEvents))

        for i,dmeson in enumerate(self.fChain):
            if i % 10000 == 0:
                print("D meson candidate n. {0}".format(i))
                if self.fMaxEvents > 0 and i > self.fMaxEvents:
                    print("Stopping the analysis.")
                    break

            self.RecalculateWeight()
            if dmeson.DmesonJet.fReconstructed.fPt > 0:
                weff = self.fEfficiencyWeight.GetEfficiencyWeight(dmeson.DmesonJet.fReconstructed)
            else:
                weff = 1.
            for r in response.itervalues():
                r.Fill(dmeson, self.fWeight, weff)

        return response
