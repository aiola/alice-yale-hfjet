#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import RawYieldSpectrumLoader

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"


def GetBFeedDownSpectra(kincuts):
    loader = RawYieldSpectrumLoader.RawYieldSpectrumLoader(input_path, "Jets_EMC_pp_1116_1117_1118_1119", "LHC10_Train1116_efficiency")
    loader.fDMeson = "D0_D0toKpiCuts"
    loader.fJetType = "Charged"
    loader.fJetRadius = "R040"
    loader.fVariableName = "JetZ"
    loader.fKinematicCuts = kincuts
    loader.fRawYieldMethod = "SideBand"
    loader.fFDConfig = { "file_name": "BFeedDown_1505317519_1399.root",
                        "central_points": "default",
                        "spectrum": "DetectorLevel_JetZSpectrum_bEfficiencyMultiply_cEfficiencyDivide"}
    h = loader.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    hFDsub = loader.GetDefaultSpectrumFromMultiTrial(True, 0, 0)
    hFD = loader.GetFDCorrection()
    hFD_up = loader.GetFDCorrection(1)
    hFD_down = loader.GetFDCorrection(-1)
    hFDsyst = loader.GetFDCorrection("graph")
    return h, hFDsub, hFD, hFD_up, hFD_down, hFDsyst


def PrepareRatio(h, hFD, hFDsyst, marker, color, fillstyle, opt):
    ratio = hFD.DrawCopy(opt)
    ratio.Divide(h)
    ratio.SetMarkerStyle(marker)
    ratio.SetMarkerColor(color)
    ratio.SetLineColor(color)
    globalList.append(ratio)

    ratio.GetYaxis().SetTitle("B Feed-Down Fraction")
    ratio.GetXaxis().SetTitle("#it{z}_{||,D}^{ch jet}")
    ratio.GetXaxis().SetTitleFont(43)
    ratio.GetXaxis().SetTitleSize(26)
    ratio.GetXaxis().SetLabelFont(43)
    ratio.GetXaxis().SetLabelSize(22)
    ratio.GetYaxis().SetTitleFont(43)
    ratio.GetYaxis().SetTitleSize(26)
    ratio.GetYaxis().SetLabelFont(43)
    ratio.GetYaxis().SetLabelSize(22)
    ratio.GetYaxis().SetTitleOffset(0.9)
    ratio.GetYaxis().SetRangeUser(0, 1.3)

    ratioSyst = hFDsyst.Clone("ratioSyst")
    globalList.append(ratioSyst)
    ratioSyst.Draw("2")
    ratioSyst.SetLineColor(color)
    ratioSyst.SetFillColor(color)
    ratioSyst.SetLineWidth(0)
    ratioSyst.SetFillStyle(fillstyle)

    for ibin in range(1, h.GetNbinsX() + 1):
        # print(hFDsyst.GetY()[ibin - 1], hFD.GetBinContent(ibin))
        ratioSyst.SetPoint(ibin - 1, ratioSyst.GetX()[ibin - 1], ratioSyst.GetY()[ibin - 1] / h.GetBinContent(ibin))
        print(ratioSyst.GetErrorYlow(ibin - 1), ratioSyst.GetErrorYhigh(ibin - 1))
        ratioSyst.SetPointEYlow(ibin - 1, ratioSyst.GetErrorYlow(ibin - 1) / h.GetBinContent(ibin))
        ratioSyst.SetPointEYhigh(ibin - 1, ratioSyst.GetErrorYhigh(ibin - 1) / h.GetBinContent(ibin))
        print(ratioSyst.GetErrorYlow(ibin - 1), ratioSyst.GetErrorYhigh(ibin - 1))

    return ratio


def PlotBFeedDown():
    h_low, hFDsub_low, hFD_low, hFD_up_low, hFD_down_low, hFDsyst_low = GetBFeedDownSpectra("DPt_20_JetPt_5_15")
    h_high, hFDsub_high, hFD_high, hFD_up_high, hFD_down_high, hFDsyst_high = GetBFeedDownSpectra("DPt_60_JetPt_15_30")

    canvas = ROOT.TCanvas("BFeedDown_JetZSpectrum_Paper", "BFeedDown_JetZSpectrum_Paper")
    globalList.append(canvas)
    canvas.SetTicks(1, 1)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    canvas.SetBottomMargin(0.15)
    canvas.cd()

    ratio_low = PrepareRatio(h_low, hFD_low, hFDsyst_low, ROOT.kFullCircle, ROOT.kGreen + 2, 3245, "p")
    ratio_high = PrepareRatio(h_high, hFD_high, hFDsyst_high, ROOT.kFullSquare, ROOT.kOrange + 2, 3254, "psame")

    paveALICE = ROOT.TPaveText(0.17, 0.74, 0.56, 0.94, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(21)
    paveALICE.SetTextAlign(13)
    # paveALICE.AddText("ALICE Preliminary")
    paveALICE.AddText("pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4, |#eta_{jet}| < 0.5")
    paveALICE.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and charge conj.")
    # paveALICE.AddText("Raw B Feed-Down Fraction from POWHEG+PYTHIA6")
    # paveALICE.AddText("Not corrected for reconstruction efficiency and")
    # paveALICE.AddText("jet momentum resolution")
    paveALICE.Draw()

    leg1 = ROOT.TLegend(0.17, 0.54, 0.56, 0.72, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(21)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    leg1.AddEntry(ratio_low, "5 < #it{p}_{T,ch jet} < 15 GeV/#it{c}, #it{p}_{T,D} > 2 GeV/#it{c}", "p")
    leg1.AddEntry(ratio_high, "15 < #it{p}_{T,ch jet} < 30 GeV/#it{c}, #it{p}_{T,D} > 6 GeV/#it{c}", "p")
    entry = leg1.AddEntry(None, "POWHEG Systematic Uncertainty", "f")
    entry.SetLineColor(ROOT.kBlack)
    entry.SetFillColor(ROOT.kBlack)
    entry.SetFillStyle(3244)
    leg1.Draw()

    return canvas


def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    canvas = PlotBFeedDown()
    canvas.SaveAs("{0}/BFeedDown_JetZSpectrum_Paper.pdf".format(input_path))
    canvas.SaveAs("{0}/BFeedDown_JetZSpectrum_Paper.C".format(input_path))


if __name__ == '__main__':
    main()

    IPython.embed()
