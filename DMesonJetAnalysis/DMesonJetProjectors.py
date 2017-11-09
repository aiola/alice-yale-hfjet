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
            print("Could not open file '{0}'! Efficiency could not be loaded!".format(filename))
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

class DMesonJetProjector:
    def __init__(self, inputPath, train, fileName, taskName, merging_type, norm_factor, maxEvents):
        self.fInputPath = inputPath
        self.fTrain = train
        self.fFileName = fileName
        self.fTaskName = taskName
        self.fMaxEvents = maxEvents
        self.fChain = None
        self.fCurrentFileName = None
        self.fCurrentRootFile = None
        self.fCurrentTreeNumber = -1
        self.fPeriod = None
        self.fPtHardBin = None
        self.fWeight = 1
        self.fNFiles = 0
        self.fMergingType = merging_type
        self.fHistEvents = None
        self.fNormFactor = norm_factor

    def GenerateChain(self, treeName):
        self.fChain = ROOT.TChain(treeName)
        path = "{0}/{1}".format(self.fInputPath, self.fTrain)

        print("Looking for file {0} in path {1}".format(self.fFileName, path))
        files = DMesonJetUtils.find_file(path, self.fFileName)

        for file in files:
            print("Adding file {0}...".format(file))
            self.fChain.Add(file)
            self.fNFiles += 1

    def ExtractWeightFromHistogramListOld(self, hlist):
        xsection = hlist.FindObject("fHistXsection")
        trials = hlist.FindObject("fHistTrials")

        if not trials or not xsection:
            print("Could not find trial and x-section information!")
            hlist.Print()
            self.fWeight = 1
            return

        valNTRIALS = trials.GetBinContent(self.fPtHardBin + 1);
        valXSEC = xsection.GetBinContent(self.fPtHardBin + 1);
        if valNTRIALS > 0:
            self.fWeight = valXSEC / valNTRIALS;

    def ExtractWeightFromHistogramList(self, hlist):
        xsection = hlist.FindObject("fHistXsectionVsPtHardNoSel")
        trials = hlist.FindObject("fHistTrialsVsPtHardNoSel")

        if not trials or not xsection:
            print("Could not find trial and x-section information (not necessarily a bad thing)!")
            hlist.Print()
            self.fWeight = 1
            return

        valNTRIALS = trials.Integral();
        valXSEC = xsection.GetMean(2);
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

            if "weighted_sum_old" in self.fMergingType:
                self.ExtractWeightFromHistogramListOld(hlist)
            else:
                self.ExtractWeightFromHistogramList(hlist)

            if self.fMergingType == "average":
                averageFactor = 1. / self.fNFiles
                totEvents = self.fChain.GetEntries()
                if self.fMaxEvents > 0 and self.fMaxEvents < totEvents:
                    averageFactor *= totEvents / self.fMaxEvents
                self.fWeight *= averageFactor
            else:
                self.fWeight *= self.fNormFactor

    def GetInfoFromFileName(self, fname):
        lastSlash = fname.rfind('/')
        secondLastSlash = fname.rfind('/', 0, lastSlash - 1)
        thirdLastSlash = fname.rfind('/', 0, secondLastSlash - 1)
        ptHardBinStr = fname[secondLastSlash + 1:lastSlash]
        if ptHardBinStr.isdigit():
            self.fPtHardBin = int(ptHardBinStr)
        else:
            self.fPtHardBin = -1
        self.fPeriod = fname[thirdLastSlash + 1:secondLastSlash]

    def ExtractCurrentFileInfo(self):
        self.fCurrentRootFile = self.fChain.GetCurrentFile()
        self.fCurrentFileName = self.fCurrentRootFile.GetName()
        self.GetInfoFromFileName(self.fCurrentFileName)

    def ExtractEventsFromHistogramList(self, hlist):
        accEventsHist = hlist.FindObject("fHistNEvents")
        rejEventsHist = hlist.FindObject("fHistEventRejectionReasons")

        if not (accEventsHist and rejEventsHist):
            print("Could not find event book-keeping histogram!")
            hlist.Print()
            return 0

        accepted = accEventsHist.GetBinContent(accEventsHist.GetXaxis().FindBin("Accepted"));
        nRecoVertVz_GT10 = rejEventsHist.GetBinContent(rejEventsHist.GetXaxis().FindBin("ZVtxOutFid"))
        nNoVert = rejEventsHist.GetBinContent(rejEventsHist.GetXaxis().FindBin("NoVertex")) + \
        rejEventsHist.GetBinContent(rejEventsHist.GetXaxis().FindBin("TooFewVtxContrib"))

        events = ROOT.TH1F("histEvents", "histEvents", 4, 0, 4)
        events.Fill("Accepted", accepted)
        events.Fill("ZVtxOutFid", nRecoVertVz_GT10)
        events.Fill("NoVert", nNoVert)

        print("Total number of accepted events with reconstructed vertex Vz < 10 cm: {0:e}".format(accepted))
        print("Total number of rejected events with reconstructed vertex Vz > 10 cm: {0:e}".format(nRecoVertVz_GT10))
        print("Total number of rejected events without reconstructed vertex: {0:e}".format(nNoVert))

        return events

    def CalculateNormalizedEvents(self):
        if not self.fHistEvents: return
        self.fAcceptedEvents = self.fHistEvents.GetBinContent(self.fHistEvents.GetXaxis().FindBin("Accepted"))
        self.fZVtxOutFidEvents = self.fHistEvents.GetBinContent(self.fHistEvents.GetXaxis().FindBin("ZVtxOutFid"))
        self.fNoVertEvents = self.fHistEvents.GetBinContent(self.fHistEvents.GetXaxis().FindBin("NoVert"))

        nRecoVert = self.fAcceptedEvents + self.fZVtxOutFidEvents
        self.fNormalizedEvents = self.fAcceptedEvents + self.fNoVertEvents * (1 - self.fZVtxOutFidEvents / nRecoVert)
        self.fHistEvents.Fill("Normalized", self.fNormalizedEvents)
        print("Total number of accepted events with reconstructed vertex Vz < 10 cm: {0:e}".format(self.fAcceptedEvents))
        print("Total number of rejected events with reconstructed vertex Vz > 10 cm: {0:e}".format(self.fZVtxOutFidEvents))
        print("Total number of rejected events without reconstructed vertex: {0:e}".format(self.fNoVertEvents))
        print("Correction factor: {}".format(self.fNormalizedEvents / self.fAcceptedEvents))
        print("Total number of normalized events: {0:e}".format(self.fNormalizedEvents))

    def OnFileChange(self, DMesonDef, DMesonDefSuffix, trigger):
        if self.fCurrentTreeNumber == self.fChain.GetTreeNumber(): return

        self.fCurrentTreeNumber = self.fChain.GetTreeNumber()
        oldPtHardBin = self.fPtHardBin
        oldPeriod = self.fPeriod

        self.ExtractCurrentFileInfo()

        if oldPtHardBin == self.fPtHardBin and oldPeriod == self.fPeriod:
            return

        if trigger:
            listName = "{0}_{1}_histos".format(self.fTaskName, trigger)
        else:
            listName = "{0}_histos".format(self.fTaskName)

        hlist = self.fChain.GetCurrentFile().Get(listName)

        if not hlist:
            print("Could not get list '{0}' from file '{1}'".format(listName, self.fChain.GetCurrentFile().GetName()))
            return 1

        self.RecalculateEvents(DMesonDef, DMesonDefSuffix, trigger)
        self.RecalculateWeight(trigger)

        print("File: {}\nPeriod: {}\nPt hard bin: {}\nWeight: {} (merging type = {})".format(self.fCurrentFileName, self.fPeriod, self.fPtHardBin, self.fWeight, self.fMergingType))

    def RecalculateEvents(self, DMesonDef, DMesonDefSuffix, trigger):
        if trigger:
            listName = "{0}_{1}_histos".format(self.fTaskName, trigger)
            subListName = "histos{0}_{1}".format(self.fTaskName, trigger)
        else:
            listName = "{0}_histos".format(self.fTaskName)
            subListName = "histos{0}".format(self.fTaskName)

        if DMesonDefSuffix:
            subSubListName = "{}_{}".format(DMesonDef, DMesonDefSuffix)
        else:
            subSubListName = DMesonDef

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

        hSubSubList = hSubList.FindObject(subSubListName)

        if not hSubSubList:
            print("Could not get list '{0}' from list '{1}'".format(subSubListName, subListName))
            hSubList.Print()
            return

        print("Calculating events for period: {}".format(self.fPeriod))
        events = self.ExtractEventsFromHistogramList(hSubSubList)

        if self.fHistEvents: self.fHistEvents.Add(events)
        else: self.fHistEvents = events

    def StartProjection(self, trigger, DMesonDef, DMesonDefSuffix, output, norm=1.0):
        if trigger:
            treeName = "{}_{}_{}_{}".format(self.fTaskName, trigger, DMesonDef, DMesonDefSuffix)
        else:
            treeName = "{}_{}_{}".format(self.fTaskName, DMesonDef, DMesonDefSuffix)

        self.fHistEvents = None
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
            self.OnFileChange(DMesonDef, DMesonDefSuffix, trigger)
            output.Fill(dmesonEvent, self.fWeight * norm)

        self.CalculateNormalizedEvents()
