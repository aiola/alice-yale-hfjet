#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import DMesonJetCompare
import array

globalList = []


def ResolutionComparison(config):
    fname = "{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"])
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)

    spectrumName = "JetZSpectrum"
    jetName = "Jet_AKTChargedR040_pt_scheme"
    dmesonName = "D0_D0toKpiCuts"
    prefix = "Prompt_{}_{}_{}".format(dmesonName, jetName, spectrumName)

    pt_lim = [(0.2, 0.4, 5, 15, "DPt_20_JetPt_5_15"), (0.6, 0.8, 5, 15, "DPt_20_JetPt_5_15"), (0.8, 1.0, 5, 15, "DPt_20_JetPt_5_15"),
              (0.6, 0.8, 15, 30, "DPt_60_JetPt_15_30"), (0.8, 1.0, 15, 30, "DPt_60_JetPt_15_30")]
    histos = []
    for (minZ, maxZ, minJetPt, maxJetPt, kincuts) in pt_lim:
        resolutionName = "{0}_{1}/DetectorResponse/{0}_{1}_DetectorResponse_{2}_{3}".format(prefix, kincuts, int(minZ * 10), int(maxZ * 10))
        h = DMesonJetUtils.GetObject(file, resolutionName)
        h.SetTitle("{} < #it{{z}}_{{||}}^{{truth}} < {} and {} < #it{{p}}_{{T,ch jet}}^{{truth}} < {}".format(minZ, maxZ, minJetPt, maxJetPt))
        globalList.append(h)
        histos.append(h)

    cname = "ResolutionVsZ_Paper"
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = None
    comp.fDoSpectrumLegend = False
    comp.fLinUpperSpace = 0.50
    comp.fColors = [ROOT.kBlue + 2, ROOT.kRed + 2, ROOT.kGreen + 2,
                    ROOT.kRed + 2, ROOT.kGreen + 2]
    comp.fMarkers = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullDiamond,
                     ROOT.kOpenSquare, ROOT.kOpenDiamond]
    r = comp.CompareSpectra(histos[0], histos[1:])
    for obj in r:
        globalList.append(obj)

    canvas = comp.fCanvasSpectra
    canvas.SetTicks(1, 1)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    canvas.SetBottomMargin(0.15)
    canvas.cd()

    h = comp.fMainHistogram

    h.GetYaxis().SetTitle("probability density")
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleSize(26)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(26)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleOffset(0.9)
    h.GetXaxis().SetRangeUser(-0.4, 0.6)

    paveALICE = ROOT.TPaveText(0.14, 0.63, 0.53, 0.95, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(21)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Prompt D^{0} #rightarrow K^{-}#pi^{+} and charge conj.")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4")
    paveALICE.AddText("|#eta_{jet}| < 0.5")
    paveALICE.Draw()

    leg = ROOT.TLegend(0.13, 0.22, 0.43, 0.57)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(19)
    entry = leg.AddEntry(None, "5 < #it{p}_{T,ch jet}^{truth} < 15 GeV/#it{c}", "")
    entry = leg.AddEntry(None, "#it{p}_{T,D} > 2 GeV/#it{c}", "")
    entry = leg.AddEntry(None, "0.8 < #it{z}_{||}^{truth} < 1.0", "pe")
    entry.SetLineColor(ROOT.kGreen + 2)
    entry.SetMarkerColor(ROOT.kGreen + 2)
    entry.SetMarkerStyle(ROOT.kFullDiamond)
    entry = leg.AddEntry(None, "0.6 < #it{z}_{||}^{truth} < 0.8", "pe")
    entry.SetLineColor(ROOT.kRed + 2)
    entry.SetMarkerColor(ROOT.kRed + 2)
    entry.SetMarkerStyle(ROOT.kFullSquare)
    entry = leg.AddEntry(None, "0.2 < #it{z}_{||}^{truth} < 0.4", "pe")
    entry.SetLineColor(ROOT.kBlue + 2)
    entry.SetMarkerColor(ROOT.kBlue + 2)
    entry.SetMarkerStyle(ROOT.kFullCircle)
    leg.Draw()
    globalList.append(leg)

    leg = ROOT.TLegend(0.56, 0.31, 0.86, 0.57)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(19)
    entry = leg.AddEntry(None, "15 < #it{p}_{T,ch jet}^{truth} < 30 GeV/#it{c}", "")
    entry = leg.AddEntry(None, "#it{p}_{T,D} > 6 GeV/#it{c}", "")
    entry = leg.AddEntry(None, "0.8 < #it{z}_{||}^{truth} < 1.0", "pe")
    entry.SetLineColor(ROOT.kGreen + 2)
    entry.SetMarkerColor(ROOT.kGreen + 2)
    entry.SetMarkerStyle(ROOT.kOpenDiamond)
    entry = leg.AddEntry(None, "0.6 < #it{z}_{||}^{truth} < 0.8", "pe")
    entry.SetLineColor(ROOT.kRed + 2)
    entry.SetMarkerColor(ROOT.kRed + 2)
    entry.SetMarkerStyle(ROOT.kOpenSquare)
    leg.Draw()
    globalList.append(leg)

    return canvas


def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    f = open("LHC15i2response_Train1399_efficiency.yaml", 'r')
    config = yaml.load(f)
    f.close()

    canvas = ResolutionComparison(config)
    canvas.SaveAs("{0}/ResolutionVsZ_Paper.pdf".format(config["input_path"]))
    canvas.SaveAs("{0}/ResolutionVsZ_Paper.C".format(config["input_path"]))


if __name__ == '__main__':

    main()

    IPython.embed()
