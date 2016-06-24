#!/usr/bin/env python
#python program to project D meson jet trees into invariant mass histograms

import ROOT
import math
import os
from DMesonJetBase import *
import array
        
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
        self.fPtHardBin = int(fname[lastSlash-1:lastSlash])
        self.fPeriod = fname[secondLastSlash+1:lastSlash-1]
        
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
            for r in response.itervalues():
                r.Fill(dmeson, self.fWeight)

        return response
