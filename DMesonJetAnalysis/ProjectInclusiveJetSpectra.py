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
from enum import Enum

globalList = []

class VariableType(Enum):
    Jet = 0
    Event = 1

class Histogram:

    def __init__(self, eta_acceptance, hdef):
        self.fType = VariableType[hdef["type"]]
        self.fEtaAcceptance = eta_acceptance
        self.fVariable = hdef["variable"]
        self.fEtaDiff = hdef["eta_diff"]
        self.GenerateHistograms(hdef)

    def GenerateHistograms(self, hdef):
        if isinstance(hdef["bins"], dict):
            bins = range(hdef["bins"]["min"], hdef["bins"]["max"] + hdef["bins"]["step"], hdef["bins"]["step"])
        else:
            bins = hdef["bins"]
        self.fHistogram = ROOT.TH1D(hdef["name"], hdef["name"], len(bins) - 1, numpy.array(bins, dtype=numpy.float))
        xaxis_title = "{} ({})".format(hdef["variable_title"], hdef["units"])
        if hdef["eta_diff"]:
            yaxis_title = "#frac{{d^{{2}}#sigma}}{{d{}#it{{d#eta}}}} [mb ({})^{{-1}}]".format(hdef["variable_title"], hdef["units"])
        else:
            yaxis_title = "#frac{{d#sigma}}{{d{}}} [mb ({})^{{-1}}]".format(hdef["variable_title"], hdef["units"])
        self.fHistogram.GetXaxis().SetTitle(xaxis_title)
        self.fHistogram.GetYaxis().SetTitle(yaxis_title)
        self.fHistogramUnweighted = ROOT.TH1D("{}_Unweighted".format(hdef["name"]), "{}_Unweighted".format(hdef["name"]), len(bins) - 1, numpy.array(bins, dtype=numpy.float))
        if hdef["eta_diff"]:
            yaxis_title = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d^{{2}}N}}{{d{}d#it{{#eta}}}} ({})^{{-1}}".format(hdef["variable_title"], hdef["units"])
        else:
            yaxis_title = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{dN}}{{d{}}} ({})^{{-1}}".format(hdef["variable_title"], hdef["units"])
        self.fHistogramUnweighted.GetXaxis().SetTitle(xaxis_title)
        self.fHistogramUnweighted.GetYaxis().SetTitle(yaxis_title)

    def Fill(self, obj, w):
        v = getattr(obj, self.fVariable)
        self.fHistogram.Fill(v, w)
        self.fHistogramUnweighted.Fill(v)

    def Normalize(self, nevents):
        scale_factor = 1.0
        if self.fEtaDiff:
            scale_factor /= self.fEtaAcceptance
        self.fHistogram.Scale(scale_factor, "width")
        if nevents > 0:
            scale_factor /= nevents
        self.fHistogramUnweighted.Scale(scale_factor, "width")

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
        self.fMergingType = self.fConfig["merging_type"]
        self.fName = self.fConfig["name"]
        self.fWeight = 1.0
        self.fMaxEvents = events
        if "max_pt_hard" in self.fConfig:
            self.fMaxPtHard = self.fConfig["max_pt_hard"]
        else:
            self.fMaxPtHard = -1

        if "reject_outliers" in config and config["reject_outliers"]:
            self.fRejectOutliers = True
            self.fOutlierPtHardJetFactor = config["outlier_pt_hard_jet_factor"]
            self.fOutlierJetDef = config["outlier_jet_def"]
        else:
            self.fRejectOutliers = False

        self.fCuts = self.fConfig["cuts"]
        self.fJetBranches = self.fConfig["jets"]

        if self.fConfig["train"] == "FastSim" or self.fConfig["train"] == "FastSimOld":
            suffix = "_".join([self.fGenerator, self.fProcess, self.fTS])

            self.fOutputPath = "{}/FastSim_{}/".format(config["input_path"], suffix)
            self.fInputPath = self.fOutputPath

            if self.fConfig["train"] == "FastSimOld":
                self.fFileName = "AnalysisResults_FastSim_{}_{}_{}.root".format(self.fGenerator, self.fProcess, self.fTS)
            else:
                self.fFileName = "AnalysisResults_FastSim_{}_{}_jets.root".format(self.fGenerator, self.fProcess)
            if self.fMergingStage >= 0:
                self.fInputPath += "stage_{}/output".format(self.fMergingStage)
            else:
                self.fInputPath += "output"
            self.fTrain = ""
        elif self.fConfig["train"] == "debug":
            name = raw_input("Directory name: ")

            self.fOutputPath = "{}/{}/".format(config["input_path"], name)
            self.fInputPath = self.fOutputPath

            self.fFileName = "AnalysisResults_FastSim_{}_{}.root".format(self.fGenerator, self.fProcess)
            self.fTrain = ""
        else:
            self.fTrain = config["train"]
            self.fInputPath = config["input_path"]
            self.fOutputPath = config["input_path"]
            self.fFileName = config["file_name"]

    def Terminate(self):
        for histos in self.fHistograms.itervalues(): 
            for h in histos.itervalues():
                h.Normalize(self.fEvents)

    def Start(self):
        self.GenerateChain()
        self.GenerateHistograms()
        self.ProjectTree()
        self.Terminate()

        fname = "{}/{}/{}.root".format(self.fOutputPath, self.fTrain, self.fName)
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
        for jet_branch in self.fJetBranches:
            histograms = dict()
            for hdef in self.fConfig["histograms"]:
                h = Histogram(jet_branch["eta_acceptance"], hdef)
                histograms[hdef["name"]] = h

            self.fHistograms[jet_branch["name"]] = histograms

    def VerifyPtHardOutlier(self, event):
        pt_hard = event.Event.fPtHard
        jets = getattr(event, self.fOutlierJetDef)
        max_jet_pt = pt_hard / self.fOutlierPtHardJetFactor
        for jet in jets:
            jet_pt = jet.fPt
            if jet_pt > max_jet_pt:
                return False
        return True 

    def ProjectTree(self):
        print("Total number of events: {}".format(self.fTree.GetEntries()))
        if self.fMaxEvents > 0 and self.fMaxEvents < self.fTree.GetEntries():
            self.fEvents = self.fMaxEvents
        else:
            self.fEvents = self.fTree.GetEntries()
        for i, entry in enumerate(self.fTree):
            if self.fMaxEvents > 0 and i > self.fMaxEvents: 
                break
            if i % 10000 == 0: 
                print("Event #{}".format(i))
            self.OnFileChange()
            if self.fMergingType == "explicit_weight": 
                self.GetExplicitWeight(entry)
            if self.fMaxPtHard > 0 and entry.Event.fPtHard > self.fMaxPtHard:
                print("Skipping event {} with pT,hard = {} and weight = {}".format(i, entry.Event.fPtHard, self.fWeight))
                continue
            if self.fRejectOutliers and not self.VerifyPtHardOutlier(entry):
                print("Skipping event {} with pT,hard = {} and weight = {}".format(i, entry.Event.fPtHard, self.fWeight))
                continue
            for jet_branch in self.fJetBranches:
                histograms = self.fHistograms[jet_branch["name"]]
                for h in histograms.itervalues():
                    if h.fType != VariableType.Event: 
                        continue
                    h.Fill(entry.Event, self.fWeight)
                jets = getattr(entry, jet_branch["branch"])
                for jet in jets:
                    bad = False
                    for c in self.fCuts:
                        v = getattr(jet, c["variable"])
                        if v > c["min"] or v >= c["max"]:
                            bad = True
                            break
                    if bad:
                        break
                    for h in histograms.itervalues():
                        if h.fType != VariableType.Jet:
                            continue
                        h.Fill(jet, self.fWeight)

    def GetExplicitWeight(self, entry):
        event = entry.Event
        self.fWeight = event.fWeight / self.fEvents

    def AddToTCollection(self, container, objects):
        for objName, obj in objects.iteritems():
            if isinstance(obj, Histogram):
                container.Add(obj.fHistogram)
                container.Add(obj.fHistogramUnweighted)
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
                obj.fHistogramUnweighted.Write()
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
        xsection = hlist.FindObject("fHistXsectionNoSel")
        trials = hlist.FindObject("fHistTrialsVsPtHardNoSel")

        if not trials or not xsection:
            print("Falling back to secondary method for x-section and trials...")
            xsection = hlist.FindObject("fHistXsectionVsPtHardNoSel")
            trials = hlist.FindObject("fHistTrialsVsPtHardNoSel")

        if not trials or not xsection:
            print("Falling back to tertiary method for x-section and trials...")
            xsection = hlist.FindObject("fHistXsectionAfterSel")
            trials = hlist.FindObject("fHistTrialsAfterSel")

        if not trials or not xsection:
            print("Could not find trial and x-section information (not necessarily a bad thing)!")
            hlist.Print()
            self.fWeight = 1
            return

        valNTRIALS = trials.Integral()
        valXSEC = xsection.GetMean(2)
        if valNTRIALS > 0:
            self.fWeight = valXSEC / valNTRIALS

    def RecalculateWeight(self):
        if self.fMergingType == "simple_sum":
            self.fWeight = 1
        elif self.fMergingType == "explicit_weight":
            pass
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

        if self.fMergingType != "explicit_weight":
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
    parser.add_argument('--gen', metavar='GEN')
    parser.add_argument('--proc', metavar='PROC')
    parser.add_argument('--ts', metavar='TS')
    parser.add_argument('--stage', metavar='N',
                        default=-1, type=int)
    parser.add_argument('-e', metavar='N',
                        default=-1, type=int)
    args = parser.parse_args()

    f = open(args.config, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig, args.gen, args.proc, args.ts, args.stage, args.e)

    IPython.embed()
