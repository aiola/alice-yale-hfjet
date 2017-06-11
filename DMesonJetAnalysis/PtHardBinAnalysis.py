#!/usr/local/bin/python

import ROOT
import math
import IPython
import argparse

globalList = []

def main(ts, nPtHard, stage):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    fileNameTemp = "/Volumes/DATA/ALICE/JetResults/FastSim_pythia_charm_{ts}/{ptHard}/stage_{stage}/output/001/AnalysisResults_FastSim_pythia_charm_{ts}.root"
    outfileNameTemp = "/Volumes/DATA/ALICE/JetResults/FastSim_pythia_charm_{ts}/{fname}.pdf"

    ptHardBins = range(0, nPtHard)
    histos = []

    hname = "ChargedR040JetPt"
    htitle = hname;
    htitle += ";#it{p}_{T,ch jet} (GeV/#it{c});#frac{d#sigma}{d#it{p}_{T}} (GeV/#it{c})^{-1}"
    hJetPtTot = ROOT.TH1D(hname, htitle, 100, 0, 100)
    hJetPtTot.Sumw2()
    globalList.append(hJetPtTot)

    hname = "DPt"
    htitle = hname;
    htitle += ";#it{p}_{T,D} (GeV/#it{c});#frac{d#sigma}{d#it{p}_{T}} (GeV/#it{c})^{-1}"
    hDPtTot = ROOT.TH1D(hname, htitle, 100, 0, 100)
    hDPtTot.Sumw2()
    globalList.append(hDPtTot)

    hname = "HardPt"
    htitle = hname;
    htitle += ";#it{p}_{T,hard} (GeV/#it{c});#frac{d#sigma}{d#it{p}_{T}} (GeV/#it{c})^{-1}"
    hHardPtTot = ROOT.TH1D(hname, htitle, 1000, 0, 1000)
    hHardPtTot.Sumw2()
    globalList.append(hHardPtTot)

    colors = [ROOT.kRed + 2, ROOT.kBlue + 2, ROOT.kMagenta + 2, ROOT.kGreen + 2, ROOT.kOrange + 2, ROOT.kCyan + 2, ROOT.kPink + 1, ROOT.kTeal + 2, ROOT.kYellow + 2]

    canvasJetPt = ROOT.TCanvas("JetPtSpectra", "JetPtSpectra")
    canvasJetPt.cd()
    canvasJetPt.SetLogy()
    haxis = hJetPtTot.DrawCopy("axis")
    haxis.GetYaxis().SetRangeUser(1e-8, 10)
    globalList.append(haxis)
    globalList.append(canvasJetPt)

    canvasDPt = ROOT.TCanvas("DPtSpectra", "DPtSpectra")
    canvasDPt.cd()
    canvasDPt.SetLogy()
    haxis = hDPtTot.DrawCopy("axis")
    haxis.GetYaxis().SetRangeUser(1e-8, 10)
    globalList.append(haxis)
    globalList.append(canvasDPt)

    canvasHardPt = ROOT.TCanvas("HardPtSpectra", "HardPtSpectra")
    canvasHardPt.cd()
    canvasHardPt.SetLogy()
    hHardPtTot.GetXaxis().SetRangeUser(0, 200)
    haxis = hHardPtTot.DrawCopy("axis")
    haxis.GetYaxis().SetRangeUser(1e-8, 10)
    globalList.append(haxis)
    globalList.append(canvasHardPt)

    for ptHard in ptHardBins:
        fileName = fileNameTemp.format(ptHard=ptHard, ts=ts, stage=stage)
        print("Opening file {0}".format(fileName))
        file = ROOT.TFile(fileName)
        tlist = file.Get("AliAnalysisTaskDmesonJets_histos")
        hTrials = tlist.FindObject("fHistTrialsVsPtHardNoSel")
        hXsec = tlist.FindObject("fHistXsectionVsPtHardNoSel")
        trials = hTrials.Integral()
        xsec = hXsec.GetMean(2)
        w = xsec / trials
        print("PtHard = {0}, Trials = {1}, xsec = {2}, w = {3}".format(ptHard, trials, xsec, w))

        tree = file.Get("AliAnalysisTaskDmesonJets_D0_MCTruth")
        hname = "ChargedR040JetPt_PtHard{ptHard}".format(ptHard=ptHard)
        htitle = hname;
        htitle += ";#it{p}_{T, ch jet} (GeV/#it{c});#frac{d#sigma}{d#it{p}_{T}} (GeV/#it{c})^{-1}"
        hJetPt = ROOT.TH1D(hname, htitle, 100, 0, 100)
        hname = "DPt_PtHard{ptHard}".format(ptHard=ptHard)
        htitle = hname;
        htitle += ";#it{p}_{T, D} (GeV/#it{c});#frac{d#sigma}{d#it{p}_{T}} (GeV/#it{c})^{-1}"
        hDPt = ROOT.TH1D(hname, htitle, 100, 0, 100)
        for entry in tree:
            jetpt = entry.Jet_AKTChargedR040_pt_scheme.fPt
            hJetPt.Fill(jetpt)
            dpt = entry.DmesonJet.fPt
            hDPt.Fill(dpt)

        hJetPt.Sumw2()
        hJetPt.Scale(w)
        canvasJetPt.cd()
        hJetPt.SetMarkerStyle(ROOT.kOpenCircle)
        hJetPt.SetMarkerSize(0.6)
        hJetPt.SetMarkerColor(colors[ptHard])
        hJetPt.SetLineColor(colors[ptHard])
        hJetPt.Draw("same")
        hJetPtTot.Add(hJetPt)
        histos.append(hJetPt)

        hDPt.Sumw2()
        hDPt.Scale(w)
        canvasDPt.cd()
        hDPt.SetMarkerStyle(ROOT.kOpenCircle)
        hDPt.SetMarkerSize(0.6)
        hDPt.SetMarkerColor(colors[ptHard])
        hDPt.SetLineColor(colors[ptHard])
        hDPt.Draw("same")
        hDPtTot.Add(hDPt)
        histos.append(hDPt)

        hHardPt = tlist.FindObject("fHistEventsVsPtHard")
        hHardPt.Sumw2()
        hHardPt.Scale(w)
        canvasHardPt.cd()
        hHardPt.SetMarkerStyle(ROOT.kOpenCircle)
        hHardPt.SetMarkerSize(0.6)
        hHardPt.SetMarkerColor(colors[ptHard])
        hHardPt.SetLineColor(colors[ptHard])
        hHardPt.Draw("same")
        hHardPtTot.Add(hHardPt)
        histos.append(hHardPt)

    canvasJetPt.cd()
    hJetPtTot.SetMarkerStyle(ROOT.kOpenSquare)
    hJetPtTot.SetMarkerSize(1)
    hJetPtTot.SetMarkerColor(ROOT.kBlack)
    hJetPtTot.SetLineColor(ROOT.kBlack)
    hJetPtTot.Draw("same")
    canvasJetPt.SaveAs(outfileNameTemp.format(ts=ts, fname=canvasJetPt.GetName()))

    canvasDPt.cd()
    hDPtTot.SetMarkerStyle(ROOT.kOpenSquare)
    hDPtTot.SetMarkerSize(1)
    hDPtTot.SetMarkerColor(ROOT.kBlack)
    hDPtTot.SetLineColor(ROOT.kBlack)
    hDPtTot.Draw("same")
    canvasDPt.SaveAs(outfileNameTemp.format(ts=ts, fname=canvasDPt.GetName()))

    canvasHardPt.cd()
    hHardPt.SetMarkerStyle(ROOT.kOpenSquare)
    hHardPt.SetMarkerSize(1)
    hHardPt.SetMarkerColor(ROOT.kBlack)
    hHardPt.SetLineColor(ROOT.kBlack)
    hHardPt.Draw("same")
    canvasHardPt.SaveAs(outfileNameTemp.format(ts=ts, fname=canvasHardPt.GetName()))

    globalList.extend(histos)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Pt hard bin analysis.')
    parser.add_argument('--ts', metavar='TS',
                        default=None)
    parser.add_argument('--pthard', metavar='N',
                        default=3, type=int)
    parser.add_argument('--stage', metavar='N',
                        default=0, type=int)
    args = parser.parse_args()

    main(args.ts, args.pthard, args.stage)

    IPython.embed()
