#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"

def GetMeasuredCrossSection():
    fname = "{0}/DataSystematics_LHC10.root".format(input_path)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "FinalSpectrum/CentralPointsStatisticalUncertainty")
    if not hStat:
        print("Cannot get measured cross section with statistical uncertainty!")
        exit(1)
    hSyst = DMesonJetUtils.GetObject(file, "FinalSpectrum/CentralPointsSystematicUncertainty")
    if not hSyst:
        print("Cannot get measured cross section with systematic uncertainty!")
        exit(1)
    return hStat, hSyst

def GetTheoryCrossSection():
    fname = "{0}/PromptDJetsPrediction.root".format(input_path)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "default/JetPtSpectrum_DPt_30/GeneratorLevel_JetPtSpectrum")
    if not hStat:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)
    hSystUp = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetPtSpectrum_DPt_30//GeneratorLevel_JetPtSpectrum/GeneratorLevel_JetPtSpectrum_UpperSyst")
    if not hSystUp:
        print("Cannot get theory cross section upper systematic uncertainty!")
        exit(1)
    hSystLow = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetPtSpectrum_DPt_30//GeneratorLevel_JetPtSpectrum/GeneratorLevel_JetPtSpectrum_LowerSyst")
    if not hSystUp:
        print("Cannot get theory cross section lower systematic uncertainty!")
        exit(1)

    # scale for the bin width and the antiparticle factor
    hStat.Scale(0.5, "width")
    hSystUp.Scale(0.5, "width")
    hSystLow.Scale(0.5, "width")
    return hStat, hSystUp, hSystLow

def PlotCrossSections(dataStat, dataSyst, theoryStat, theorySystUp, theorySystLow):
    cname = "FinalSpectrum_QM17"
    canvas = ROOT.TCanvas(cname, cname, 700, 700)
    globalList.append(canvas)
    canvas.Divide(1, 2)
    padMain = canvas.cd(1)
    padMain.SetPad(0, 0.35, 1, 1)
    padMain.SetBottomMargin(0)
    padMain.SetLeftMargin(0.15)
    padMain.SetTicks(1, 1)
    padRatio = canvas.cd(2)
    padRatio.SetPad(0, 0., 1, 0.35)
    padRatio.SetTopMargin(0)
    padRatio.SetBottomMargin(0.27)
    padRatio.SetLeftMargin(0.15)
    padRatio.SetTicks(1, 1)

    padMain.cd()
    padMain.SetLogy()
    dataSyst_copy = dataSyst.DrawCopy("e2")
    globalList.append(dataSyst_copy)
    dataSyst_copy.GetYaxis().SetRangeUser(2e-5, 5e-1)
    dataSyst_copy.GetYaxis().SetTitleFont(43)
    dataSyst_copy.GetYaxis().SetTitleSize(26)
    dataSyst_copy.GetYaxis().SetLabelFont(43)
    dataSyst_copy.GetYaxis().SetLabelSize(22)
    dataSyst_copy.GetYaxis().SetTitleOffset(1.6)
    dataStat_copy = dataStat.DrawCopy("same p e0 x0")
    globalList.append(dataStat_copy)

    theoryStat_copy = theoryStat.DrawCopy("same p e0 x0")
    globalList.append(theoryStat_copy)
    theoryStat_copy.SetLineColor(ROOT.kBlue + 2)
    theoryStat_copy.SetMarkerColor(ROOT.kBlue + 2)
    theoryStat_copy.SetMarkerStyle(ROOT.kOpenCircle)
    theoryStat_copy.SetMarkerSize(1.2)

    theorySystUp_copy = theoryStat.DrawCopy("same hist l")
    globalList.append(theorySystUp_copy)
    theorySystUp_copy.Add(theorySystUp, 1)
    theorySystUp_copy.SetLineColor(ROOT.kBlue + 2)
    theorySystUp_copy.SetLineWidth(1)

    theorySystLow_copy = theoryStat.DrawCopy("same hist l")
    globalList.append(theorySystLow_copy)
    theorySystLow_copy.Add(theorySystLow, -1)
    theorySystLow_copy.SetLineColor(ROOT.kBlue + 2)
    theorySystLow_copy.SetLineWidth(1)

    padRatio.cd()
    padRatio.SetGridy()

    ratioSyst = dataSyst_copy.DrawCopy("e2")
    globalList.append(ratioSyst)
    ratioSyst.GetYaxis().SetTitle("data / theory")
    ratioSyst.GetXaxis().SetTitleFont(43)
    ratioSyst.GetXaxis().SetTitleSize(26)
    ratioSyst.GetXaxis().SetLabelFont(43)
    ratioSyst.GetXaxis().SetLabelSize(22)
    ratioSyst.GetYaxis().SetTitleFont(43)
    ratioSyst.GetYaxis().SetTitleSize(26)
    ratioSyst.GetYaxis().SetLabelFont(43)
    ratioSyst.GetYaxis().SetLabelSize(22)
    ratioSyst.GetYaxis().SetTitleOffset(1.4)
    ratioSyst.GetXaxis().SetTitleOffset(2.9)
    ratioSyst.GetYaxis().SetRangeUser(0, 3.9)
    ratioStat = dataStat_copy.DrawCopy("same p e0 x0")
    globalList.append(ratioStat)
    ratioTheorySystUp = theorySystUp_copy.DrawCopy("same hist l")
    globalList.append(ratioTheorySystUp)
    ratioTheorySystLow = theorySystLow_copy.DrawCopy("same hist l")
    globalList.append(ratioTheorySystLow)

    for ibin in range(1, ratioSyst.GetNbinsX() + 1):
        ratioSyst.SetBinContent(ibin, ratioSyst.GetBinContent(ibin) / theoryStat_copy.GetBinContent(ibin))
        ratioSyst.SetBinError(ibin, ratioSyst.GetBinError(ibin) / theoryStat_copy.GetBinContent(ibin))
        ratioStat.SetBinContent(ibin, ratioStat.GetBinContent(ibin) / theoryStat_copy.GetBinContent(ibin))
        ratioStat.SetBinError(ibin, ratioStat.GetBinError(ibin) / theoryStat_copy.GetBinContent(ibin))
        ratioTheorySystUp.SetBinContent(ibin, ratioTheorySystUp.GetBinContent(ibin) / theoryStat_copy.GetBinContent(ibin))
        ratioTheorySystLow.SetBinContent(ibin, ratioTheorySystLow.GetBinContent(ibin) / theoryStat_copy.GetBinContent(ibin))

    padMain.cd()
    leg1 = ROOT.TLegend(0.34, 0.39, 0.80, 0.65, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(23)
    leg1.SetTextAlign(13)
    leg1.SetMargin(0.2)
    leg1.AddEntry(dataStat_copy, "ALICE", "pe")
    leg1.AddEntry(dataSyst_copy, "Systematic Uncertainty (data)", "f")
    leg1.AddEntry(theoryStat_copy, "POWHEG+PYTHIA6", "pe")
    leg1.AddEntry(theorySystUp_copy, "Systematic Uncertainty (theory)", "l")
    leg1.Draw()

    padMain.cd()
    paveALICE = ROOT.TPaveText(0.16, 0.71, 0.55, 0.90, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(20)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Preliminary, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4, |#eta_{jet}| < 0.5")
    paveALICE.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c., #it{p}_{T,D} > 3 GeV/#it{c}")
    paveALICE.Draw()

    return canvas

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    dataStat, dataSyst = GetMeasuredCrossSection()
    theoryStat, theorySystUp, theorySystDown = GetTheoryCrossSection()
    canvas = PlotCrossSections(dataStat, dataSyst, theoryStat, theorySystUp, theorySystDown)
    canvas.SaveAs("{0}/D0JetCrossSection_pp7TeV.pdf".format(input_path))

if __name__ == '__main__':
    main()

    IPython.embed()
