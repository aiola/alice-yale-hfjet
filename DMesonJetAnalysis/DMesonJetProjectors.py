#!/usr/bin/env python
# python program to project D meson jet trees into invariant mass histograms

import math
import os

import ROOT

import DMesonJetUtils
import DetectorResponse

class SimpleWeight:
    def GetEfficiencyWeight(self, dmeson, jet):
        return 1
    def GetEfficiencyWeightTH1ForPt(self, pt):
        return 1

class EfficiencyWeightCalculator:
    def __init__(self, filename="", listname="", objectname=""):
        self.fRootObject = None
        if filename and objectname:
            self.LoadEfficiency(filename, listname, objectname)

    def GetObjectFromRootFile(self, file, listname, objectname):
        if listname:
            rlist = file.Get(listname)
            if not rlist:
                print("Could not find list '{0}' in file '{1}'".format(listname, file.GetName()))
                file.ls()
                return None
            obj = rlist.FindObject(objectname)
            if not obj:
                print("Could not find object '{0}' in list '{1}'".format(objectname, listname))
                rlist.Print()
                return None
        else:
            obj = file.Get(objectname)
            if not obj:
                print("Could not find object '{0}' in file '{1}'".format(objectname, file.GetName()))
                rlist.Print()
                return None
        return obj

    def LoadEfficiency(self, filename, listname, objectname):
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
        elif isinstance(self.fRootObject, ROOT.TF1):
            print("The ROOT object is of type TF1")
            self.GetEfficiencyWeight = self.GetEfficiencyWeightTF1
        elif isinstance(self.fRootObject, ROOT.TH2):
            print("The ROOT object is of type TH2")
            self.GetEfficiencyWeight = self.GetEfficiencyWeightTH2
        elif isinstance(self.fRootObject, ROOT.TH1):
            print("The ROOT object is of type TH1")
            self.GetEfficiencyWeight = self.GetEfficiencyWeightTH1
        else:
            print("The ROOT object type is not recognized!!")
            self.fRootObject.Print()

    def GetEfficiencyWeightTGraph(self, dmeson, jet):
        eff = self.fRootObject.Eval(dmeson.fPt)
        if eff == 0:
            return 0
        else:
            return 1. / eff

    def GetEfficiencyWeightTF1(self, dmeson, jet):
        eff = self.fRootObject.Eval(dmeson.fPt)
        if eff == 0:
            return 0
        else:
            return 1. / eff

    def GetEfficiencyWeightTH2(self, dmeson, jet):
        eff = self.fRootObject.GetBinContent(self.fRootObject.FindBin(jet.fPt, dmeson.fPt))

        if eff == 0:
            return 0
        else:
            return 1. / eff

    def GetEfficiencyWeightTH1(self, dmeson, jet):
        eff = self.fRootObject.GetBinContent(self.fRootObject.FindBin(dmeson.fPt))

        if eff == 0:
            return 0
        else:
            return 1. / eff

    def GetEfficiencyWeightTH1ForPt(self, pt):
        eff = self.fRootObject.GetBinContent(self.fRootObject.FindBin(pt))
        if eff == 0:
            return 0
        else:
            return 1. / eff

class DMesonJetDataProjector:
    def __init__(self, inputPath, train, fileName, taskName, merging_type, maxEvents):
        self.fInputPath = inputPath
        self.fTrain = train
        self.fFileName = fileName
        self.fTaskName = taskName
        self.fMaxEvents = maxEvents
        self.fMassAxisTitle = "#it{m}(K#pi) GeV/#it{c}^{2}"
        self.fYieldAxisTitle = "counts"
        self.fChain = None
        self.fCurrentFileName = None
        self.fPeriod = None
        self.fWeight = 1
        self.fTotalEvents = 0
        self.fNFiles = 0
        self.fMergingType = merging_type

    def GenerateChain(self, treeName):
        path = "{0}/{1}".format(self.fInputPath, self.fTrain)

        print("Looking for file {0} in path {1}".format(self.fFileName, path))
        files = DMesonJetUtils.find_file(path, self.fFileName)

        self.fChain = ROOT.TChain(treeName)

        for file in files:
            print("Adding file {0}...".format(file))
            self.fChain.Add(file)
            self.fNFiles += 1

    def ExtractWeightFromHistogramList(self, hlist):
        xsection = hlist.FindObject("fHistXsectionAfterSel")
        trials = hlist.FindObject("fHistTrialsAfterSel")

        if not trials or not xsection:
            print("Could not find trial and x-section information (not necessarily a bad thing)!")
            hlist.Print()
            self.fWeight = 1
            return

        valNTRIALS = trials.Integral();
        valXSEC = xsection.GetMean(2);
        scalingFactor = 0;
        if valNTRIALS > 0:
            self.fWeight = valXSEC / valNTRIALS;

    def RecalculateWeight(self, trigger):
        if self.fMergingType == "simple_sum":
            self.fWeight = 1
        else:
            if trigger:
                listName = "{0}_{1}_histos".format(self.fTaskName, trigger)
            else:
                listName = "{0}_histos".format(self.fTaskName)
            hlist = self.fChain.GetCurrentFile().Get(listName)

            if not hlist:
                print("Could not get list '{0}' from file '{1}'".format(listName, self.fChain.GetCurrentFile().GetName()))
                self.fWeight = 1

            self.ExtractWeightFromHistogramList(hlist)

            if self.fMergingType == "average":
                averageFactor = 1. / self.fNFiles
                if self.fMaxEvents > 0:
                    averageFactor *= self.fChain.GetEntries() / self.fMaxEvents
                self.fWeight *= averageFactor

        print("File: {0}\nWeight: {1} (merging type = {2})".format(self.fCurrentFileName, self.fWeight, self.fMergingType))

    def ExtractCurrentFileInfo(self):
        fname = self.fChain.GetCurrentFile().GetName()
        lastSlash = fname.rfind('/')
        secondLastSlash = fname.rfind('/', 0, lastSlash - 1)
        thirdLastSlash = fname.rfind('/', 0, secondLastSlash - 1)
        self.fPeriod = fname[thirdLastSlash + 1:secondLastSlash]

    def ExtractEventsFromHistogramList(self, hlist):
        eventsHist = hlist.FindObject("fHistNEvents")

        if not eventsHist:
            print("Could not find event book-keeping histogram!")
            hlist.Print()
            return 0

        events = eventsHist.GetBinContent(eventsHist.GetXaxis().FindBin("Accepted"));
        return events

    def OnFileChange(self, DMesonDef, trigger):
        if self.fChain.GetCurrentFile().GetName() == self.fCurrentFileName:
            return

        self.fCurrentFileName = self.fChain.GetCurrentFile().GetName()

        self.ExtractCurrentFileInfo()
        self.RecalculateEvents(DMesonDef, trigger)
        self.RecalculateWeight(trigger)

    def RecalculateEvents(self, DMesonDef, trigger):
        if trigger:
            listName = "{0}_{1}_histos".format(self.fTaskName, trigger)
            subListName = "histos{0}_{1}".format(self.fTaskName, trigger)
        else:
            listName = "{0}_histos".format(self.fTaskName)
            subListName = "histos{0}".format(self.fTaskName)

        hlist = self.fChain.GetCurrentFile().Get(listName)

        if not hlist:
            print("Could not get list '{0}' from file '{1}'".format(listName, self.fChain.GetCurrentFile().GetName()))
            self.fChain.GetCurrentFile().ls()
            return

        hSubList = hlist.FindObject(subListName)

        if not hSubList:
            print("Could not get list '{0}' from list '{1}'".format(subListName, listName))
            hlist.Print()
            return

        hSubSubList = hSubList.FindObject(DMesonDef)

        if not hSubSubList:
            print("Could not get list '{0}' from list '{1}'".format(DMesonDef, subListName))
            hSubList.Print()
            return

        events = self.ExtractEventsFromHistogramList(hSubSubList)
        self.fTotalEvents += events

        print("Period: {0}\nEvents: {1}".format(self.fPeriod, events))

    def StartProjection(self, trigger, DMesonDef, binMultiSets, nMassBins, minMass, maxMass):
        self.fTotalEvents = 0
        if trigger:
            treeName = "{0}_{1}_{2}".format(self.fTaskName, trigger, DMesonDef)
        else:
            treeName = "{0}_{1}".format(self.fTaskName, DMesonDef)

        self.GenerateChain(treeName)

        print("Running analysis on tree {0}. Total number of entries is {1}".format(treeName, self.fChain.GetEntries()))
        if self.fMaxEvents > 0:
            print("The analysis will stop at the {0} entry.".format(self.fMaxEvents))

        for i, dmesonEvent in enumerate(self.fChain):
            if i % 10000 == 0:
                print("D meson candidate n. {0}".format(i))
                if self.fMaxEvents > 0 and i >= self.fMaxEvents:
                    print("Stopping the analysis.")
                    break
            dmeson = dmesonEvent.DmesonJet
            self.OnFileChange(DMesonDef, trigger)
            for (jtype, jradius), binMultiSet in binMultiSets.iteritems():
                if jtype or jradius:
                    jetName = "Jet_AKT{0}{1}_pt_scheme".format(jtype, jradius)
                    jet = getattr(dmesonEvent, jetName)
                    if jet.fPt == 0:
                        continue
                else:
                    jet = None

                bins = binMultiSet.FindBin(dmeson, jet, DMesonDef)
                for bin, weight in bins:
                    if bin.fCounts == 0: bin.CreateHistograms(trigger, DMesonDef, self.fMassAxisTitle, self.fYieldAxisTitle, nMassBins, minMass, maxMass)
                    bin.Fill(dmeson, jet, weight * self.fWeight)

                spectra = binMultiSet.FindSpectra(dmeson, jet)
                for spectrum, weight in spectra: spectrum.Fill(dmeson, jet, weight * self.fWeight)

        print("Total number of events: {0}".format(self.fTotalEvents))

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
        self.fCurrentFileName = ""

    def GenerateChain(self, treeName):
        path = "{0}/{1}".format(self.fInputPath, self.fTrain)

        files = DMesonJetUtils.find_file(path, self.fFileName)

        self.fChain = ROOT.TChain(treeName)

        for file in files:
            print("Adding file {0}...".format(file))
            self.fChain.Add(file)

    def ExtractCurrentFileInfo(self):
        if self.fCurrentFileName == self.fChain.GetCurrentFile().GetName(): return
        self.fCurrentFileName = self.fChain.GetCurrentFile().GetName()
        lastSlash = self.fCurrentFileName.rfind('/')
        secondLastSlash = self.fCurrentFileName.rfind('/', 0, lastSlash - 1)
        thirdLastSlash = self.fCurrentFileName.rfind('/', 0, secondLastSlash - 1)
        self.fPtHardBin = int(self.fCurrentFileName[secondLastSlash + 1:lastSlash])
        self.fPeriod = self.fCurrentFileName[thirdLastSlash + 1:secondLastSlash]

    def ExtractWeightFromHistogramList(self, hlist):
        xsection = hlist.FindObject("fHistXsection")
        trials = hlist.FindObject("fHistTrials")

        if not trials or not xsection:
            print("Could not find trial and x-section information!")
            hlist.Print()
            self.fWeight = 1
            return

        valNTRIALS = trials.GetBinContent(self.fPtHardBin + 1);
        valXSEC = xsection.GetBinContent(self.fPtHardBin + 1);
        scalingFactor = 0;
        if valNTRIALS > 0:
            self.fWeight = valXSEC / valNTRIALS;

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
            for axisName, (axisDef, weightEff, applyEffTruth, cuts) in respDefinitions.iteritems():
                respName = "{0}_{1}_{2}".format(DMesonDef, jetName, axisName)
                resp = DetectorResponse.DetectorResponse(respName, jetName, axisDef, cuts, weightEff, applyEffTruth)
                resp.GenerateHistograms()
                response[respName] = resp
                treeName = "{0}_{1}".format(self.fTaskName, DMesonDef)

        self.GenerateChain(treeName)

        print("Running analysis on tree {0}. Total number of entries is {1}".format(treeName, self.fChain.GetEntries()))
        if self.fMaxEvents > 0:
            print("The analysis will stop at the {0} entry.".format(self.fMaxEvents))

        for i, dmeson in enumerate(self.fChain):
            if i % 10000 == 0:
                print("D meson candidate n. {0}".format(i))
                if self.fMaxEvents > 0 and i > self.fMaxEvents:
                    print("Stopping the analysis.")
                    break

            self.RecalculateWeight()
            for r in response.itervalues():
                r.Fill(dmeson, self.fWeight)

        return response
