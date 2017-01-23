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
    hFD = loader.GetFDCorrection()
    hFD_up = loader.GetFDCorrection(1)
    hFD_down = loader.GetFDCorrection(-1)
    return h, hFD, hFD_up, hFD_down

def PlotBFeedDown():
    h, hFD, hFD_up, hFD_down = GetBFeedDownSpectra()

    cname = "BFeedDown_QM17"
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
    h_axis = h.DrawCopy("axis")
    globalList.append(h_axis)
    h_axis.GetYaxis().SetTitle("yield (arb. units)")
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

    ratio_up = hFD_up_copy.DrawCopy("same hist l")
    globalList.append(ratio_up)
    ratio_up.Divide(h_copy)

    ratio_down = hFD_down_copy.DrawCopy("same hist l")
    globalList.append(ratio_up)
    ratio_down.Divide(h_copy)

    ratio.GetYaxis().SetTitle("FD Fraction")
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
    ratio.GetYaxis().SetRangeUser(0, 0.49)

    padMain.cd()
    leg1 = ROOT.TLegend(0.42, 0.52, 0.88, 0.70, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(23)
    leg1.SetTextAlign(13)
    leg1.SetMargin(0.2)
    leg1.AddEntry(h_copy, "Uncorrected D^{0}-Jet Yield", "pe")
    leg1.AddEntry(hFD_copy, "B Feed-Down (POWHEG)", "pe")
    leg1.AddEntry(ratio_down, "Systematic Uncertainty", "l")
    leg1.Draw()

    return canvas

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    canvas = PlotBFeedDown()
    canvas.SaveAs("{0}/BFeedDown_QM17.pdf".format(input_path))

if __name__ == '__main__':
    main()

    IPython.embed()
