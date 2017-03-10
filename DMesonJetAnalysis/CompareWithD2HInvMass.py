#!/usr/bin/env python

import ROOT
import math
import IPython
# import argparse
# import DMesonJetCompare
import DMesonJetUtils

globalList = []


def GetD2HInvMass(file):
    listname = "PWG3_D2H_D0InvMassLoose/coutputmassD0MassLoose0100"
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
    for obj in tree:
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

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    fileName = "../anaDev/AnalysisResults.root"

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

if __name__ == '__main__':

    # parser = argparse.ArgumentParser(description='Compare statistical uncertainties.')
    # parser.add_argument('ts', metavar='ts',
    #                    help='Timestamps', nargs='*')
    # args = parser.parse_args()

    main()

    IPython.embed()
