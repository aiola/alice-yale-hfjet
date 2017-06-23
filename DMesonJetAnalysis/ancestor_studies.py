#!/usr/bin/env python

import ROOT
import math
import IPython
import argparse
import DMesonJetUtils

globalList = []

partonDict = { 4:"charm", 5:"beauty" }
ancestorDict = { 2212:"proton", 1:"uds", 2:"uds", 3:"uds", 4:"charm", 5:"beauty", 21:"gluon" }

def ProjectMyTree(tree, resp, hname):
    uniquePartons = set(partonDict.values())
    uniqueAncestor = set(ancestorDict.values())
    hQuarkVsAncestor = ROOT.TH2F("hQuarkVsAncestor", "hQuarkVsAncestor;ancestor;quark;counts", len(uniqueAncestor) + 1, 0, len(uniqueAncestor) + 1, len(uniquePartons) + 1, 0, len(uniquePartons) + 1)
    hQuarkVsAncestor.GetXaxis().SetBinLabel(1, "unknown");
    hQuarkVsAncestor.GetYaxis().SetBinLabel(1, "unknown");
    for i, b in enumerate(uniqueAncestor):
        hQuarkVsAncestor.GetXaxis().SetBinLabel(i + 2, b);
    for i, b in enumerate(uniquePartons):
        hQuarkVsAncestor.GetYaxis().SetBinLabel(i + 2, b);

    for i, obj in enumerate(tree):
        if i % 10000 == 0:
            print("Event {}".format(i))
        if resp:
            partonType = obj.DmesonJet.fGenerated.fPartonType
            ancestorPdg = obj.DmesonJet.fGenerated.fAncestorPDG
        else:
            partonType = obj.DmesonJet.fPartonType
            ancestorPdg = obj.DmesonJet.fAncestorPDG
        if partonType in partonDict:
            parton = partonDict[partonType]
        else:
            parton = "unknown"
            print("Unknown parton type: {}".format(partonType))
        if ancestorPdg in ancestorDict:
            ancestor = ancestorDict[ancestorPdg]
        else:
            ancestor = "unknown"
            print("Unknown ancestor PDG: {}".format(ancestorPdg))
        hQuarkVsAncestor.Fill(ancestor, parton, 1)
    return hQuarkVsAncestor

def GetQuarkVsAncestor(file, resp, dmeson_name):
    if resp:
        treename = "AliAnalysisTaskDmesonJetsDetectorResponse_{}".format(dmeson_name)
    else:
        treename = "AliAnalysisTaskDmesonJets_{}_MCTruth".format(dmeson_name)
    tree = DMesonJetUtils.GetObject(file, treename)
    h = ProjectMyTree(tree, resp, "hQuarkVsAncestor")
    return h

def single_file_analysis(path, filename, resp, dmeson_name):
    hQuarkVsAncestor = None
    file_names = DMesonJetUtils.find_file(path, filename)

    for fileName in file_names:
        print("Opening file {0}".format(fileName))
        file = ROOT.TFile(fileName)

        hQuarkVsAncestor_part = GetQuarkVsAncestor(file, resp, dmeson_name)
        hQuarkVsAncestor_part.Sumw2()
        if hQuarkVsAncestor:
            hQuarkVsAncestor.Add(hQuarkVsAncestor_part)
        else:
            hQuarkVsAncestor = hQuarkVsAncestor_part

    PlotQuarkVsAncestor(hQuarkVsAncestor)

def PlotQuarkVsAncestor(hQuarkVsAncestor):
    hQuarkVsAncestor.Scale(1. / hQuarkVsAncestor.Integral())
    hQuarkVsAncestor.GetZaxis().SetTitle("fraction")
    hQuarkVsAncestor.GetZaxis().SetRangeUser(1e-6, 1)

    c = ROOT.TCanvas("hQuarkVsAncestor", "hQuarkVsAncestor")
    c.SetLogz()
    hQuarkVsAncestor.Draw("colz text")

    globalList.append(c)
    globalList.append(hQuarkVsAncestor)

def GetWeight(file):
    tlist = file.Get("AliAnalysisTaskDmesonJets_histos")
    hTrials = tlist.FindObject("fHistTrialsVsPtHardNoSel")
    hXsec = tlist.FindObject("fHistXsectionVsPtHardNoSel")
    trials = hTrials.Integral()
    xsec = hXsec.GetMean(2)
    w = xsec / trials
    return w

def pthard_analysis(path, fname, resp, dmeson_name):
    hQuarkVsAncestor = None

    print("Looking for file '{}' in '{}'".format(fname, path))
    file_names = DMesonJetUtils.find_file(path, fname)
    for fileName in file_names:
        print("Opening file {0}".format(fileName))
        file = ROOT.TFile(fileName)

        hQuarkVsAncestor_pthard = GetQuarkVsAncestor(file, resp, dmeson_name)
        hQuarkVsAncestor_pthard.Sumw2()

        w = GetWeight(file)

        hQuarkVsAncestor_pthard.Scale(w)

        if hQuarkVsAncestor:
            hQuarkVsAncestor.Add(hQuarkVsAncestor_pthard)
        else:
            hQuarkVsAncestor = hQuarkVsAncestor_pthard

    PlotQuarkVsAncestor(hQuarkVsAncestor)

def ParseTrainNumbers(train):
    TrainNumberList = train.split(",")
    TrainNumbers = []
    for TrainNumberRange in TrainNumberList:
        Range = TrainNumberRange.split(":")
        if len(Range) > 1:
            TrainNumbers.extend(range(int(Range[0]), int(Range[len(Range) - 1]) + 1))
        else:
            TrainNumbers.extend(Range)
    TrainNumbersLabel = "_".join([str(i) for i in TrainNumbers])
    return TrainNumbersLabel

def main(train, test, pthard, train_name, fname, resp, dmeson_name):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if test:
        path = "../anaDev"
        single_file_analysis(path, fname, resp, dmeson_name)
    else:
        TrainNumbersLabel = ParseTrainNumbers(train)

        path = "/Volumes/DATA/ALICE/JetResults/{}_{}".format(train_name, TrainNumbersLabel)
        if pthard:
            pthard_analysis(path, fname, resp, dmeson_name)
        else:
            single_file_analysis(path, fname, resp, dmeson_name)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Ancestor studies.')
    parser.add_argument('--train',
                        help='Train number(s) (use : for a range and , for a list)')
    parser.add_argument('--test', action='store_const',
                        default=False, const=True,
                        help='Ignore train number to use ../anaDev/AnalysisResults.root')
    parser.add_argument('--pthard', action='store_const',
                        default=False, const=True,
                        help='Pt hard bin production')
    parser.add_argument('--resp', action='store_const',
                        default=False, const=True,
                        help='Pt hard bin production')
    parser.add_argument('--train-name',
                        default="Jets_EMC_pp_MC",
                        help='Train name')
    parser.add_argument('--d-meson',
                        default="D0",
                        help='D meson')
    parser.add_argument('--file-name',
                        default="AnalysisResults.root",
                        help='File name')
    args = parser.parse_args()

    main(args.train, args.test, args.pthard, args.train_name, args.file_name, args.resp, args.d_meson)

    IPython.embed()
