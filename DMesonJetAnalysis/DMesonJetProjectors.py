#!/usr/bin/env python
#python program to project D meson jet trees into invariant mass histograms

import ROOT
import math
import os
from DMesonJetBase import *
import array
from bisect import bisect

class SimpleWeight:
    def GetEfficiencyWeight(self, dmeson, jetDef):
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
        if isinstance(self.fRootObject, ROOT.TGraph):
            print("The ROOT object is of type TGraph")
            self.GetEfficiencyWeight = self.GetEfficiencyWeightTGraph
        elif isinstance(self.fRootObject, ROOT.TH2):
            print("The ROOT object is of type TH2")
            self.GetEfficiencyWeight = self.GetEfficiencyWeightTH2
        elif isinstance(self.fRootObject, ROOT.TH1):
            print("The ROOT object is of type TH1")
            self.GetEfficiencyWeight = self.GetEfficiencyWeightTH1

    def GetEfficiencyWeightTGraph(self, dmeson, jetDef):
        eff = self.fRootObject.Eval(dmeson.DmesonJet.fReconstructed.fPt)
        #print("pT,D = {0}, eff = {1}".format(dmeson.DmesonJet.fReconstructed.fPt, eff))
        if eff == 0:
            return 0
        else:
            return 1. / eff

    def GetEfficiencyWeightTH2(self, dmeson, jetDef):
        jet = getattr(dmeson, "{0}_reco".format(jetDef))
#        if jet.fPt < self.fRootObject.GetXaxis().GetBinLowEdge(1) or \
#           jet.fPt >= self.fRootObject.GetXaxis().GetBinUpEdge(self.fRootObject.GetXaxis().GetNbins()) or \
#           dmeson.DmesonJet.fReconstructed.fPt < self.fRootObject.GetYaxis().GetBinLowEdge(1) or \
#           dmeson.DmesonJet.fReconstructed.fPt >= self.fRootObject.GetYaxis().GetBinUpEdge(self.fRootObject.GetYaxis().GetNbins()):
#            eff = 0
#        else:
#            eff = self.fRootObject.Interpolate(jet.fPt, dmeson.DmesonJet.fReconstructed.fPt)
        eff = self.fRootObject.GetBinContent(self.fRootObject.FindBin(jet.fPt, dmeson.DmesonJet.fReconstructed.fPt))
        #print("pT,D = {0}, pTjet = {1}, eff = {2}".format(dmeson.DmesonJet.fReconstructed.fPt, jet.fPt, eff))

        if eff == 0:
            return 0
        else:
            return 1. / eff

    def GetEfficiencyWeightTH1(self, dmeson, jetDef):
#        if dmeson.DmesonJet.fReconstructed.fPt < self.fRootObject.GetXaxis().GetBinLowEdge(1) or \
#           dmeson.DmesonJet.fReconstructed.fPt >= self.fRootObject.GetXaxis().GetBinUpEdge(self.fRootObject.GetXaxis().GetNbins()):
#            eff = 0
#        else:
#            eff = self.fRootObject.Interpolate(dmeson.DmesonJet.fReconstructed.fPt)
        eff = self.fRootObject.GetBinContent(self.fRootObject.FindBin(dmeson.DmesonJet.fReconstructed.fPt))
        #print("pT,D = {0}, eff = {1}".format(dmeson.DmesonJet.fReconstructed.fPt, eff))

        if eff == 0:
            return 0
        else:
            return 1. / eff

class DMesonJetDataProjector:
    def __init__(self, inputPath, train, fileName, taskName, maxEvents):
        self.fInputPath = inputPath
        self.fTrain = train
        self.fFileName = fileName
        self.fTaskName = taskName
        self.fMaxEvents = maxEvents
        self.fMassAxisTitle = "#it{m}(K#pi) GeV/#it{c}^{2}"
        self.fYieldAxisTitle = "counts"
        self.fChain = None

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
                    bin.fInvMassHisto.Fill(dmeson.DmesonJet.fInvMass)

class DMesonJetResponseProjector:
    def __init__(self, inputPath, train, fileName, taskName, maxEvents):
        self.fInputPath = inputPath
        self.fTrain = train
        self.fFileName = fileName
        self.fTaskName = taskName
        self.fMaxEvents = maxEvents
        self.fChain = None
        self.fPtHardBin = -1
        self.fPeriod = ""
        self.fWeight = 1

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
            for axisName,(axisDef, weightEff) in respDefinitions.iteritems():
                respName = "{0}_{1}_{2}".format(DMesonDef, jetName, axisName)
                response[respName] = DetectorResponse(respName, jetName, axisDef, weightEff)
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
            for r in response.itervalues():
                r.Fill(dmeson, self.fWeight)

        return response
