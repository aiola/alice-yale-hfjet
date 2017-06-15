#!/usr/bin/env python
# python script to project TTree containing inclusive jet spectra

import argparse
import IPython
import ROOT
import os
import DMesonJetCompare
import DMesonJetUtils
import numpy

globalList = []

def main(path, train):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    tree = GenerateChain("AliAnalysisTaskEmcalJetTree_INT7_jets", path, train, "AnalysisResults_jets.root")

    centralityBins = [0, 10, 30, 50]
    jetPtBins = [0, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 150]
    jetCorrPtBins = [-80, -60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 70]
    leadingPtBins = [0, 3, 5]
    jetType = "Charged"
    jetAlgo = "AKT"
    jetRadius = "R040"

    histos = ProjectTree(tree, jetType, jetAlgo, jetRadius, centralityBins, jetPtBins, jetCorrPtBins, leadingPtBins)

    fname = "{}/{}/Jets.root".format(path, train)
    file = ROOT.TFile(fname, "recreate")
    AddToTDirectory(file, histos)
    file.Close()
    print("Results stored in {}".format(fname))

def AddToTCollection(container, objects):
    for objName, obj in objects.iteritems():
        if isinstance(obj, ROOT.TObject):
            container.Add(obj)
        elif isinstance(obj, dict):
            sub_container = ROOT.TList()
            sub_container.SetName(objName)
            AddToTCollection(sub_container, obj)
            container.Add(sub_container)
        else:
            print("Object type '{}' not recognized".format(obj))

def AddToTDirectory(container, objects):
    container.cd()
    for objName, obj in objects.iteritems():
        if isinstance(obj, ROOT.TCollection):
            obj.Write(objName, ROOT.TObject.kSingleKey)
        elif isinstance(obj, ROOT.TObject):
            obj.Write()
        elif isinstance(obj, dict):
            sub_container = ROOT.TList()
            sub_container.SetName(str(objName))
            AddToTCollection(sub_container, obj)
            sub_container.Write(str(objName), ROOT.TObject.kSingleKey)
        else:
            print("Object type '{}' not recognized".format(obj))

def ProjectTree(tree, jetType, jetAlgo, jetRadius, centralityBins, jetPtBins, jetCorrPtBins, leadingPtBins):
    if jetType == "Charged":
        jetBranch = "Jet_{}{}{}_tracks_pT0150_pt_scheme".format(jetAlgo, jetType, jetRadius)
    else:
        print("Not implemented for jet type different from 'Charged'!")
        exit(1)

    histos = dict()
    histos["CentRejected"] = ROOT.TH1F("CentRejected", "CentRejected", 100, 0, 100)

    for cmin, cmax in zip(centralityBins[:-1], centralityBins[1:]):
        for leadPt in leadingPtBins:
            chistos = dict()
            hname = "JetPt_Cent{}_{}_LeadHadPt{}".format(cmin, cmax, leadPt)
            chistos[hname] = ROOT.TH1F(hname, hname, len(jetPtBins) - 1, numpy.array(jetPtBins, dtype=numpy.float32))
            hname = "JetCorrPt_Cent{}_{}_LeadHadPt{}".format(cmin, cmax, leadPt)
            chistos[hname] = ROOT.TH1F(hname, hname, len(jetPtBins) - 1, numpy.array(jetCorrPtBins, dtype=numpy.float32))

            histos[(cmin, cmax), leadPt] = chistos

    print("Total number of events: {}".format(tree.GetEntries()))
    for i, entry in enumerate(tree):
        event = entry.Event
        if i % 1000 == 0: print("Event #{}".format(i))
        centOk = False
        for cmin, cmax in zip(centralityBins[:-1], centralityBins[1:]):
            if event.fCent >= cmin and event.fCent < cmax:
                centOk = True
                break
        if not centOk:
            histos["CentRejected"].Fill(event.fCent)
            continue
        jets = getattr(entry, jetBranch)
        for jet in jets:
            for leadPt in leadingPtBins:
                if jet.fLeadingPt >= leadPt:
                    chistos = histos[(cmin, cmax), leadPt]
                    hJetPt = chistos["JetPt_Cent{}_{}_LeadHadPt{}".format(cmin, cmax, leadPt)]
                    hJetCorrPt = chistos["JetCorrPt_Cent{}_{}_LeadHadPt{}".format(cmin, cmax, leadPt)]
                    hJetPt.Fill(jet.fPt)
                    hJetCorrPt.Fill(jet.fCorrPt)
    return histos

def GenerateChain(treeName, rootPath, train, fileName):
    chain = ROOT.TChain(treeName)

    path = "{}/{}".format(rootPath, train)

    print("Looking for file {0} in path {1}".format(fileName, path))
    files = DMesonJetUtils.find_file(path, fileName)

    for file in files:
        print("Adding file {0}...".format(file))
        chain.Add(file)

    return chain

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Project TTree containing inclusive jet spectra.')
    parser.add_argument('train',
                        help='Train name (e.g. Jets_EMC_PbPb_1872)')
    parser.add_argument('--path', metavar='PATH',
                        default="/Volumes/DATA/ALICE/JetResults")
    args = parser.parse_args()

    main(args.path, args.train)

    IPython.embed()
