#!/usr/bin/env python

import ROOT
import math
import IPython
import argparse
# import DMesonJetCompare
import DMesonJetUtils

globalList = []


def GetD2HInvMass(file):
    listname = "PWG3_D2H_D0InvMass/coutputmassD0Mass0100"
    hlist = DMesonJetUtils.GetObject(file, listname)
    if not hlist:
        print("Could not get list {}".format(listname))
        exit(1)
    hSgn = None
    hRfl = None
    for i in range(0, 14):
        h = hlist.FindObject("histSgn_{}".format(i))
        if hSgn: hSgn.Add(h)
        else: hSgn = h.Clone("hSgn")
        h = hlist.FindObject("histRfl_{}".format(i))
        if hRfl: hRfl.Add(h)
        else: hRfl = h.Clone("hRfl")
    hSgn.Rebin(6)
    hRfl.Rebin(6)
    return hSgn, hRfl

def ProjectMyTree(tree):
    h = ROOT.TH1F("h", "h", 100, 1.6248, 2.2248)
    for i, obj in enumerate(tree):
        if i % 1000 == 0:
            print("Event {}".format(i))
        invMass = obj.DmesonJet.fInvMass
        h.Fill(invMass)
    h.Sumw2()
    return h

def GetMyInvMass(file):
    treename = "AliAnalysisTaskDmesonJets_D0_kSignalOnly"
    tree = DMesonJetUtils.GetObject(file, treename)
    hSgn = ProjectMyTree(tree)
    hSgn.SetName("hSgn")

    treename = "AliAnalysisTaskDmesonJets_D0_WrongPID"
    tree = DMesonJetUtils.GetObject(file, treename)
    hRfl = ProjectMyTree(tree)
    hRfl.SetName("hRfl")

    return hSgn, hRfl

def GetWeight(file):
    tlist = file.Get("AliAnalysisTaskDmesonJets_histos")
    hTrials = tlist.FindObject("fHistTrialsVsPtHard")
    hXsec = tlist.FindObject("fHistXsectionVsPtHard")
    trials = hTrials.Integral()
    xsec = hXsec.GetMean(2)
    w = xsec / trials
    return w

def single_file_analysis(filename):
    print("Opening file {0}".format(fileName))
    file = ROOT.TFile(fileName)

    hSgn_D2H, hRfl_D2H = GetD2HInvMass(file)
    hSgn, hRfl = GetMyInvMass(file)

    globalList.append(hSgn_D2H)
    globalList.append(hRfl_D2H)
    globalList.append(hSgn)
    globalList.append(hRfl)

    cSgn = ROOT.TCanvas("cSgn", "cSgn")
    globalList.append(cSgn)
    cSgn.cd()
    hSgn_D2H.Draw("hist")
    hSgn_D2H.SetMarkerSize(0.9)
    hSgn_D2H.SetMarkerColor(ROOT.kRed)
    hSgn_D2H.SetMarkerStyle(ROOT.kOpenCircle)
    hSgn_D2H.SetLineColor(ROOT.kRed)
    hSgn.Draw("same hist")
    hSgn.SetMarkerSize(0.9)
    hSgn.SetMarkerColor(ROOT.kBlue)
    hSgn.SetMarkerStyle(ROOT.kOpenSquare)
    hSgn.SetLineColor(ROOT.kBlue)

    cRfl = ROOT.TCanvas("cRfl", "cRfl")
    globalList.append(cRfl)
    cRfl.cd()
    hRfl_D2H.Draw("hist")
    hRfl_D2H.SetMarkerSize(0.9)
    hRfl_D2H.SetMarkerColor(ROOT.kRed)
    hRfl_D2H.SetMarkerStyle(ROOT.kOpenCircle)
    hRfl_D2H.SetLineColor(ROOT.kRed)
    hRfl.Draw("same hist")
    hRfl.SetMarkerSize(0.9)
    hRfl.SetMarkerColor(ROOT.kBlue)
    hRfl.SetMarkerStyle(ROOT.kOpenSquare)
    hRfl.SetLineColor(ROOT.kBlue)

def pthard_analysis(path, periods, pt_hard_bins):
    hSgn_D2H = ROOT.TH1F("h", "h;M (GeV/#it{c}^{2});arb. units", 100, 1.6248, 2.2248)
    hRfl_D2H = ROOT.TH1F("h", "h;M (GeV/#it{c}^{2});arb. units", 100, 1.6248, 2.2248)
    hSgn = ROOT.TH1F("h", "h;M (GeV/#it{c}^{2});arb. units", 100, 1.6248, 2.2248)
    hRfl = ROOT.TH1F("h", "h;M (GeV/#it{c}^{2});arb. units", 100, 1.6248, 2.2248)

    for pt_hard in pt_hard_bins:
        for period in periods:
            fileName = "{}/{}/{}/AnalysisResults.root".format(path, period, pt_hard)
            print("Opening file {0}".format(fileName))
            file = ROOT.TFile(fileName)

            hSgn_D2H_pthard, hRfl_D2H_pthard = GetD2HInvMass(file)
            hSgn_pthard, hRfl_pthard = GetMyInvMass(file)

            w = GetWeight(file)

            hSgn_D2H.Add(hSgn_D2H_pthard, w)
            hRfl_D2H.Add(hRfl_D2H_pthard, w)
            hSgn.Add(hSgn_pthard, w)
            hRfl.Add(hRfl_pthard, w)

    globalList.append(hSgn_D2H)
    globalList.append(hRfl_D2H)
    globalList.append(hSgn)
    globalList.append(hRfl)

    cSgn = ROOT.TCanvas("cSgn", "cSgn")
    globalList.append(cSgn)
    cSgn.cd()
    hSgn_D2H.Draw("hist")
    hSgn_D2H.SetMarkerSize(0.9)
    hSgn_D2H.SetMarkerColor(ROOT.kRed)
    hSgn_D2H.SetMarkerStyle(ROOT.kOpenCircle)
    hSgn_D2H.SetLineColor(ROOT.kRed)
    hSgn.Draw("same hist")
    hSgn.SetMarkerSize(0.9)
    hSgn.SetMarkerColor(ROOT.kBlue)
    hSgn.SetMarkerStyle(ROOT.kOpenSquare)
    hSgn.SetLineColor(ROOT.kBlue)

    cRfl = ROOT.TCanvas("cRfl", "cRfl")
    globalList.append(cRfl)
    cRfl.cd()
    hRfl_D2H.Draw("hist")
    hRfl_D2H.SetMarkerSize(0.9)
    hRfl_D2H.SetMarkerColor(ROOT.kRed)
    hRfl_D2H.SetMarkerStyle(ROOT.kOpenCircle)
    hRfl_D2H.SetLineColor(ROOT.kRed)
    hRfl.Draw("same hist")
    hRfl.SetMarkerSize(0.9)
    hRfl.SetMarkerColor(ROOT.kBlue)
    hRfl.SetMarkerStyle(ROOT.kOpenSquare)
    hRfl.SetLineColor(ROOT.kBlue)

def main(train):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if train == "test":
        fileName = "../anaDev/AnalysisResults.root"
        single_file_analysis(fileName)
    else:
        trainno = int(train)
        path = "/Volumes/DATA/ALICE/JetResults/Jets_EMC_pp_MC_{}_{}_{}_{}".format(trainno, trainno + 1, trainno + 2, trainno + 3)
        periods = ["LHC15i2b", "LHC15i2c", "LHC15i2d", "LHC15i2e"]
        pt_hard_bins = range(1, 9)
        pthard_analysis(path, periods, pt_hard_bins)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compare invariant mass with D2H.')
    parser.add_argument('train', metavar='train',
                        help='Train number (just the first number, e.g. 1033 -> Jets_EMC_pp_MC_1033_1034_1035_1036. Type test to use ../anaDev/AnalysisResults.root')
    args = parser.parse_args()

    main(args.train)

    IPython.embed()
