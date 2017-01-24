#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import RawYieldSpectrumLoader

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"

def GetBFeedDownSpectra():
    loader = RawYieldSpectrumLoader.RawYieldSpectrumLoader(input_path, "Jets_EMC_pp_823_824_825_826", "LHC10_Train823_efficiency")
    loader.fDMeson = "D0"
    loader.fJetType = "Charged"
    loader.fJetRadius = "R040"
    loader.fSpectrumName = "JetPtSpectrum"
    loader.fKinematicCuts = "DPt_30"
    loader.fRawYieldMethod = "SideBand"
    h = loader.GetDefaultSpectrumFromMultiTrial(False)
    hFDsub = loader.GetDefaultSpectrumFromMultiTrial(True)
    hFD = loader.GetFDCorrection()
    hFD_up = loader.GetFDCorrection(1)
    hFD_down = loader.GetFDCorrection(-1)
    return h, hFDsub, hFD, hFD_up, hFD_down

def PlotBFeedDown():
    h, hFDsub, hFD, hFD_up, hFD_down = GetBFeedDownSpectra()

    cname = "BFeedDown_Internal"
    canvas1 = ROOT.TCanvas(cname, cname, 700, 700)
    globalList.append(canvas1)
    canvas1.Divide(1, 2)
    padMain = canvas1.cd(1)
    padMain.SetPad(0, 0.35, 1, 1)
    padMain.SetBottomMargin(0)
    padMain.SetLeftMargin(0.15)
    padMain.SetTicks(1, 1)
    padRatio = canvas1.cd(2)
    padRatio.SetPad(0, 0., 1, 0.35)
    padRatio.SetTopMargin(0)
    padRatio.SetBottomMargin(0.27)
    padRatio.SetLeftMargin(0.15)
    padRatio.SetTicks(1, 1)

    padMain.cd()
    padMain.SetLogy()
    h_axis = h.DrawCopy("axis")
    globalList.append(h_axis)
    h_axis.GetYaxis().SetTitle("yield (counts #times efficiency)")
    h_axis.GetYaxis().SetRangeUser(5e1, 3e4)
    h_axis.GetYaxis().SetTitleFont(43)
    h_axis.GetYaxis().SetTitleSize(26)
    h_axis.GetYaxis().SetLabelFont(43)
    h_axis.GetYaxis().SetLabelSize(22)
    h_axis.GetYaxis().SetTitleOffset(1.6)

    h_copy = h.DrawCopy("same p")
    h_copy.SetLineColor(ROOT.kBlue + 2)
    h_copy.SetMarkerColor(ROOT.kBlue + 2)
    h_copy.SetMarkerStyle(ROOT.kFullCircle)
    globalList.append(h_copy)

    hFDsub_copy = hFDsub.DrawCopy("same p")
    hFDsub_copy.SetLineColor(ROOT.kGreen + 2)
    hFDsub_copy.SetMarkerColor(ROOT.kGreen + 2)
    hFDsub_copy.SetMarkerStyle(ROOT.kFullSquare)
    globalList.append(hFDsub_copy)

    hFD_copy = hFD.DrawCopy("same p")
    globalList.append(hFD_copy)
    hFD_copy.SetLineColor(ROOT.kRed + 2)
    hFD_copy.SetMarkerColor(ROOT.kRed + 2)
    hFD_copy.SetMarkerStyle(ROOT.kOpenCircle)
    hFD_copy.SetMarkerSize(1.2)

    hFD_up_copy = hFD_up.DrawCopy("same hist l")
    globalList.append(hFD_up_copy)
    hFD_up_copy.SetLineColor(ROOT.kRed + 2)

    hFD_down_copy = hFD_down.DrawCopy("same hist l")
    globalList.append(hFD_down_copy)
    hFD_down_copy.SetLineColor(ROOT.kRed + 2)

    padRatio.cd()
    padRatio.SetGridy()

    ratio = hFD_copy.DrawCopy("p")
    ratio.Divide(h_copy)
    globalList.append(ratio)

    ratio_FDsub = hFDsub_copy.DrawCopy("same hist")
    ratio_FDsub.SetLineWidth(2)
    ratio_FDsub.Divide(h_copy)
    globalList.append(ratio_FDsub)

    ratio_up = hFD_up_copy.DrawCopy("same hist l")
    globalList.append(ratio_up)
    ratio_up.Divide(h_copy)

    ratio_down = hFD_down_copy.DrawCopy("same hist l")
    globalList.append(ratio_up)
    ratio_down.Divide(h_copy)

    ratio.GetYaxis().SetTitle("ratio")
    ratio.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
    ratio.GetXaxis().SetTitleFont(43)
    ratio.GetXaxis().SetTitleSize(26)
    ratio.GetXaxis().SetLabelFont(43)
    ratio.GetXaxis().SetLabelSize(22)
    ratio.GetYaxis().SetTitleFont(43)
    ratio.GetYaxis().SetTitleSize(26)
    ratio.GetYaxis().SetLabelFont(43)
    ratio.GetYaxis().SetLabelSize(22)
    ratio.GetYaxis().SetTitleOffset(1.4)
    ratio.GetXaxis().SetTitleOffset(2.9)
    ratio.GetYaxis().SetRangeUser(0, 0.99)

    padMain.cd()
    leg1 = ROOT.TLegend(0.42, 0.60, 0.88, 0.80, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(23)
    leg1.SetTextAlign(13)
    leg1.SetMargin(0.2)
    leg1.AddEntry(h_copy, "Uncorrected D^{0}-Jet Yield", "pe")
    leg1.AddEntry(hFDsub_copy, "FD Subtracted Yield", "pe")
    leg1.AddEntry(hFD_copy, "B Feed-Down (POWHEG)", "pe")
    leg1.AddEntry(ratio_down, "Systematic Uncertainty", "l")
    leg1.Draw()

    canvas2 = ROOT.TCanvas("BFeedDown_QM17", "BFeedDown_QM17")
    globalList.append(canvas2)
    canvas2.SetTicks(1, 1)
    canvas2.SetLeftMargin(0.13)
    canvas2.SetRightMargin(0.05)
    canvas2.SetTopMargin(0.05)
    canvas2.SetBottomMargin(0.15)
    canvas2.cd()

    ratio2 = hFD_copy.DrawCopy("p")
    ratio2.Divide(h_copy)
    globalList.append(ratio2)

    ratio2.GetYaxis().SetTitle("B Feed-Down Fraction")
    ratio2.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
    ratio2.GetXaxis().SetTitleFont(43)
    ratio2.GetXaxis().SetTitleSize(26)
    ratio2.GetXaxis().SetLabelFont(43)
    ratio2.GetXaxis().SetLabelSize(22)
    ratio2.GetYaxis().SetTitleFont(43)
    ratio2.GetYaxis().SetTitleSize(26)
    ratio2.GetYaxis().SetLabelFont(43)
    ratio2.GetYaxis().SetLabelSize(22)
    ratio2.GetYaxis().SetTitleOffset(0.9)
    ratio2.GetYaxis().SetRangeUser(0, 0.69)

    ratio2_up = hFD_up_copy.DrawCopy("same hist l")
    globalList.append(ratio2_up)
    ratio2_up.Divide(h_copy)

    ratio2_down = hFD_down_copy.DrawCopy("same hist l")
    globalList.append(ratio2_down)
    ratio2_down.Divide(h_copy)

    paveALICE = ROOT.TPaveText(0.17, 0.61, 0.56, 0.93, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(20)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Preliminary, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4, |#eta_{jet}| < 0.5")
    paveALICE.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c., #it{p}_{T,D} > 3 GeV/#it{c}")
    paveALICE.AddText("Raw B Feed-Down Fraction from POWHEG+PYTHIA6")
    paveALICE.AddText("Not corrected for reconstruction efficiency and")
    paveALICE.AddText("jet momentum resolution")
    paveALICE.Draw()

    return canvas1, canvas2

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    canvas1, canvas2 = PlotBFeedDown()
    canvas1.SaveAs("{0}/BFeedDown_InternalPlot.pdf".format(input_path))
    canvas2.SaveAs("{0}/BFeedDown_QM17.pdf".format(input_path))

if __name__ == '__main__':
    main()

    IPython.embed()