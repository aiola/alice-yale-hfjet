#!/usr/bin/env python
# python script to project TTree containing inclusive jet spectra

import argparse
import IPython
import ROOT
import os
import DMesonJetCompare
import DMesonJetUtils
import numpy
import yaml

globalList = []


class Histogram:

    def __init__(self, var, h):
        self.fVariable = var
        self.fHistogram = h

    def Fill(self, obj, w):
        v = getattr(obj, self.fVariable)
        self.fHistogram.Fill(v, w)

    def Normalize(self):
        self.fHistogram.Scale(1.0, "width")


class ProjectInclusiveJetSpectra:

    def __init__(self, config, gen, proc, ts, stage, events):
        self.fConfig = config
        self.fGenerator = gen
        self.fProcess = proc
        self.fTS = ts
        self.fMergingStage = stage
        self.fCurrentTreeNumber = -1
        self.fTaskName = self.fConfig["task_name"]
        self.fTreeName = "{}_jets".format(self.fTaskName)
        self.fMergingType = config["merging_type"]
        self.fWeight = 1.0
        self.fMaxEvents = events

        self.fCuts = self.fConfig["cuts"]
        self.fJetBranch = self.fConfig["jet_branch"]

        if self.fConfig["train"] == "FastSim":
            suffix = "_".join([self.fGenerator, self.fProcess, self.fTS])

            self.fName = "_".join([self.fConfig["name"], suffix])

            self.fOutputPath = "{}/FastSim_{}/".format(config["input_path"], suffix)
            self.fInputPath = self.fOutputPath

            if self.fMergingStage >= 0:
                self.fInputPath += "stage_{}/output".format(self.fMergingStage)
                self.fFileName = "AnalysisResults_FastSim_{}.root".format(suffix)
            else:
                self.fInputPath += "output"
                self.fFileName = "AnalysisResults_FastSim_{}_{}.root".format(self.fGenerator, self.fProcess)
            self.fTrain = ""
        else:
            self.fName = config["name"]
            self.fTrain = config["train"]
            self.fInputPath = config["input_path"]
            self.fOutputPath = config["input_path"]
            self.fFileName = config["file_name"]

    def Terminate(self):
        for h in self.fHistograms.itervalues(): h.Normalize()

    def Start(self):
        self.GenerateChain()
        self.GenerateHistograms()
        self.ProjectTree()
        self.Terminate()

        fname = "{}/{}/Jets.root".format(self.fOutputPath, self.fTrain)
        file = ROOT.TFile(fname, "recreate")
        self.AddToTDirectory(file, self.fHistograms)
        file.Close()
        print("Results stored in {}".format(fname))

    def GenerateChain(self):
        self.fTree = ROOT.TChain(self.fTreeName)

        path = "{}/{}".format(self.fInputPath, self.fTrain)

        print("Looking for file {} in path {}".format(self.fFileName, path))
        files = DMesonJetUtils.find_file(path, self.fFileName)
        self.fNFiles = 0

        for file in files:
            print("Adding file {0}...".format(file))
            self.fTree.Add(file)
            self.fNFiles += 1

    def GenerateHistograms(self):
        self.fHistograms = dict()
        for hdef in config["histograms"]:
            hobj = ROOT.TH1D(hdef["name"], hdef["title"], len(hdef["bins"]) - 1, numpy.array(hdef["bins"], dtype=numpy.float32))
            h = Histogram(hdef["variable"], hobj)
            self.fHistograms[hdef["name"]] = h

    def ProjectTree(self):
        print("Total number of events: {}".format(self.fTree.GetEntries()))
        for i, entry in enumerate(self.fTree):
            if self.fMaxEvents > 0 and i > self.fMaxEvents: break
            self.OnFileChange()
            if i % 10000 == 0: print("Event #{}".format(i))
            jets = getattr(entry, self.fJetBranch)
            for jet in jets:
                bad = False
                for c in self.fCuts:
                    v = getattr(jet, c["variable"])
                    if v > c["min"] or v >= c["max"]:
                        bad = True
                        break
                if bad: break
                for h in self.fHistograms.itervalues():
                    h.Fill(jet, self.fWeight)

    def AddToTCollection(self, container, objects):
        for objName, obj in objects.iteritems():
            if isinstance(obj, Histogram):
                container.Add(obj.fHistogram)
            elif isinstance(obj, ROOT.TObject):
                container.Add(obj)
            elif isinstance(obj, dict):
                sub_container = ROOT.TList()
                sub_container.SetName(objName)
                self.AddToTCollection(sub_container, obj)
                container.Add(sub_container)
            else:
                print("Object type '{}' not recognized".format(obj))

    def AddToTDirectory(self, container, objects):
        container.cd()
        for objName, obj in objects.iteritems():
            if isinstance(obj, ROOT.TCollection):
                obj.Write(objName, ROOT.TObject.kSingleKey)
            elif isinstance(obj, Histogram):
                obj.fHistogram.Write()
            elif isinstance(obj, ROOT.TObject):
                obj.Write()
            elif isinstance(obj, dict):
                sub_container = ROOT.TList()
                sub_container.SetName(str(objName))
                self.AddToTCollection(sub_container, obj)
                sub_container.Write(str(objName), ROOT.TObject.kSingleKey)
            else:
                print("Object type '{}' not recognized".format(obj))

    def ExtractWeightFromHistogramList(self, hlist):
        xsection = hlist.FindObject("fHistXsectionVsPtHardNoSel")
        trials = hlist.FindObject("fHistTrialsVsPtHardNoSel")

        if not trials or not xsection:
            print("Falling back to secondary method for x-section and trials...")
            xsection = hlist.FindObject("fHistXsectionAfterSel")
            trials = hlist.FindObject("fHistTrialsAfterSel")

        if not trials or not xsection:
            print("Could not find trial and x-section information (not necessarily a bad thing)!")
            hlist.Print()
            self.fWeight = 1
            return

        valNTRIALS = trials.Integral();
        valXSEC = xsection.GetMean(2);
        if valNTRIALS > 0:
            self.fWeight = valXSEC / valNTRIALS;

    def RecalculateWeight(self):
        if self.fMergingType == "simple_sum":
            self.fWeight = 1
        else:
            listName = "{0}_histos".format(self.fTaskName)
            hlist = self.fTree.GetCurrentFile().Get(listName)

            if not hlist:
                print("Could not get list '{0}' from file '{1}'".format(listName, self.fTree.GetCurrentFile().GetName()))
                self.fWeight = 1

            self.ExtractWeightFromHistogramList(hlist)

            if self.fMergingType == "average":
                averageFactor = 1. / self.fNFiles
                self.fWeight *= averageFactor

    def OnFileChange(self):
        if self.fCurrentTreeNumber == self.fTree.GetTreeNumber(): return

        self.fCurrentTreeNumber = self.fTree.GetTreeNumber()

        listName = "{0}_histos".format(self.fTaskName)

        hlist = self.fTree.GetCurrentFile().Get(listName)

        if not hlist:
            print("Could not get list '{0}' from file '{1}'".format(listName, self.fTree.GetCurrentFile().GetName()))
            return

        self.RecalculateWeight()

        print("File: {}\nWeight: {} (merging type = {})".format(self.fTree.GetCurrentFile().GetName(), self.fWeight, self.fMergingType))


def main(config, gen, proc, ts, stage, events):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    eng = ProjectInclusiveJetSpectra(config, gen, proc, ts, stage, events)
    eng.Start()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Project TTree containing inclusive jet spectra.')
    parser.add_argument('config',
                        help='YAML file')
    parser.add_argument('--gen', metavar='GEN',
                        default=None)
    parser.add_argument('--proc', metavar='PROC',
                        default=None)
    parser.add_argument('--ts', metavar='TS',
                        default=None)
    parser.add_argument('--stage', metavar='N',
                        default=-1, type=int)
    parser.add_argument('-e', metavar='N',
                        default=-1, type=int)
    args = parser.parse_args()

    f = open(args.config, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.gen, args.proc, args.ts, args.stage, args.e)

    IPython.embed()
