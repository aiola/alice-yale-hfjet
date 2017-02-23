#!/usr/bin/env python

import ROOT
import math
import IPython

globalList = []

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    # fileNameTemp = "/Volumes/DATA/ALICE/JetResults/FastSim_pythia_charm_1487281181/{ptHard}/stage_0/output/001/AnalysisResults_FastSim_pythia_charm_1487281181.root"
    fileNameTemp = "/Volumes/DATA/ALICE/JetResults/FastSim_pythia_charm_1487374803/{ptHard}/stage_0/output/001/AnalysisResults_FastSim_pythia_charm_1487374803.root"

    ptHardBins = range(1, 10)
    histos = []
    hTot = ROOT.TH1D("htot", "htot", 100, 0, 100)
    hTot.Sumw2()
    globalList.append(hTot)

    colors = [ROOT.kRed + 2, ROOT.kBlue + 2, ROOT.kMagenta + 2, ROOT.kGreen + 2, ROOT.kOrange + 2, ROOT.kCyan + 2, ROOT.kPink + 1, ROOT.kTeal + 2, ROOT.kYellow + 2]
    canvas = ROOT.TCanvas("mycanv", "mycanv")
    canvas.cd()
    canvas.SetLogy()
    haxis = hTot.DrawCopy("axis")
    globalList.append(haxis)
    globalList.append(canvas)

    for ptHard in ptHardBins:
        fileName = fileNameTemp.format(ptHard=ptHard)
        print("Opening file {0}".format(fileName))
        file = ROOT.TFile(fileName)
        tlist = file.Get("AliAnalysisTaskDmesonJets_histos")
        hTrials = tlist.FindObject("fHistTrialsAfterSel")
        hXsec = tlist.FindObject("fHistXsectionAfterSel")
        trials = hTrials.Integral()
        xsec = hXsec.GetBinContent(ptHard + 1)
        w = xsec / trials
        print("PtHard = {0}, Trials = {1}, xsec = {2}, w = {3}".format(ptHard, trials, xsec, w))

        tree = file.Get("AliAnalysisTaskDmesonJets_D0_MCTruth")
        hname = "pthard{ptHard}".format(ptHard=ptHard)
        h = ROOT.TH1D(hname, hname, 100, 0, 100)
        for entry in tree:
            pt = entry.Jet_AKTChargedR040_pt_scheme.fPt
            h.Fill(pt)
        h.Sumw2()
        h.Scale(w)
        canvas.cd()
        h.SetMarkerStyle(ROOT.kOpenCircle)
        h.SetMarkerSize(0.6)
        h.SetMarkerColor(colors[ptHard])
        h.SetLineColor(colors[ptHard])
        h.Draw("same")
        hTot.Add(h)
        histos.append(h)
    print("max = {0}, min = {1}".format(hTot.GetMaximum(0), hTot.GetMinimum(0)))
    hTot.SetMarkerStyle(ROOT.kOpenSquare)
    hTot.SetMarkerSize(1)
    hTot.SetMarkerColor(ROOT.kBlack)
    hTot.SetLineColor(ROOT.kBlack)
    hTot.Draw("same")
    # span = hTot.GetMaximum(0) / hTot.GetMinimum(0)
    # haxis.SetMinimum(hTot.GetMinimum(0) / span ** 0.3)
    # haxis.SetMaximum(hTot.GetMaximum(0) * span ** 0.3)
    haxis.GetYaxis().SetRangeUser(1e-8, 10)
    globalList.extend(histos)

if __name__ == '__main__':

    main()

    IPython.embed()
