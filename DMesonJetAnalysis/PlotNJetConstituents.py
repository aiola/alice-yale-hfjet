#!/usr/bin/env python
# python script to prepare B feed-down correction file

import argparse
import IPython
import ROOT
import array
import numpy
import math
import yaml
import DMesonJetUtils
import RawYieldSpectrumLoader

globalList = []

def main(data_config, mc_config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    data_loader = RawYieldSpectrumLoader.RawYieldSpectrumLoader(data_config["input_path"], data_config["train"], data_config["name"])
    data_loader.fUseReflections = False
    data_loader.fDMeson = "D0"
    data_loader.fJetType = "Charged"
    data_loader.fJetRadius = "R040"
    data_loader.fSpectrumName = "JetPtNConst"
    data_loader.fKinematicCuts = "DPt_30"
    data_loader.fRawYieldMethod = "SideBand"

    h = data_loader.GetDefaultSpectrumFromDMesonJetAnalysis()
    h2D_data = h.Clone("JetPtNConst_Data")

    mc_loader = RawYieldSpectrumLoader.RawYieldSpectrumLoader(mc_config["input_path"], mc_config["train"], mc_config["name"])
    mc_loader.fUseReflections = False
    mc_loader.fDMeson = "D0_kSignalOnly"
    mc_loader.fJetType = "Charged"
    mc_loader.fJetRadius = "R040"
    mc_loader.fSpectrumName = "JetPtNConst"
    mc_loader.fKinematicCuts = "DPt_30"
    mc_loader.fRawYieldMethod = "SignalOnly"

    h = mc_loader.GetDefaultSpectrumFromDMesonJetAnalysis()
    h2D_mc = h.Clone("JetPtNConst_MCdet")

    mc_loader_truth = RawYieldSpectrumLoader.RawYieldSpectrumLoader(mc_config["input_path"], mc_config["train"], mc_config["name"])
    mc_loader_truth.fUseReflections = False
    mc_loader_truth.fDMeson = "D0_MCTruth"
    mc_loader_truth.fJetType = "Charged"
    mc_loader_truth.fJetRadius = "R040"
    mc_loader_truth.fSpectrumName = "JetPtNConst"
    mc_loader_truth.fKinematicCuts = "DPt_30"
    mc_loader_truth.fRawYieldMethod = "Truth"

    h = mc_loader_truth.GetDefaultSpectrumFromDMesonJetAnalysis()
    h2D_mc_truth = h.Clone("JetPtNConst_MCgen")

    if h2D_data.GetNbinsX() <= 3:
        ncols = h2D_data.GetNbinsX()
        nrows = 1
    else:
        ncols = int((math.sqrt(h2D_data.GetNbinsX())).ceil())
        nrows = int((h2D_data.GetNbinsX() / ncols).floor())
    cname = "{0}/NJetConstituents_{1}_{2}".format(data_config["input_path"], data_config["name"], mc_config["name"])
    canvas = ROOT.TCanvas(cname, cname, ncols * 400, nrows * 400)
    canvas.Divide(ncols, nrows)
    globalList.append(canvas)

    histos_data = Plot(canvas, h2D_data, "", ROOT.kFullCircle, ROOT.kBlue + 2)
    histos_mc = Plot(canvas, h2D_mc, "same", ROOT.kOpenCircle, ROOT.kGreen + 2)
    histos_mctruth = Plot(canvas, h2D_mc_truth, "same", ROOT.kOpenSquare, ROOT.kRed + 2)

    canvas.cd(1)
    leg = ROOT.TLegend(0.50, 0.72, 0.72, 0.87, "", "NB NDC")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(16)
    leg.SetMargin(0.2)
    leg.AddEntry(histos_data[0], "Data", "pe")
    leg.AddEntry(histos_mc[0], "MC, det. level", "pe")
    leg.AddEntry(histos_mctruth[0], "MC, gen. level", "pe")
    leg.Draw()
    globalList.append(leg)

    canvas.SaveAs("{0}.pdf".format(canvas.GetName()))

def Plot(canvas, h2D, opt, marker, color):
    histos = []
    for ibin in range(1, h2D.GetNbinsX() + 1):
        ptmin = h2D.GetXaxis().GetBinLowEdge(ibin)
        ptmax = h2D.GetXaxis().GetBinUpEdge(ibin)
        h = h2D.ProjectionY("{0}_JetPt_{1:.0f}_{2:.0f}".format(h2D.GetName(), ptmin * 100, ptmax * 100), ibin, ibin)
        pad = canvas.cd(ibin)
        pad.SetLeftMargin(0.18)
        pad.SetRightMargin(0.02)
        pad.SetTopMargin(0.09)
        pad.SetBottomMargin(0.13)
        h_copy = h.DrawCopy(opt)
        globalList.append(h_copy)
        histos.append(h_copy)
        h_copy.GetYaxis().SetTitle("Probability Density")
        h_copy.GetXaxis().SetTitleFont(43)
        h_copy.GetXaxis().SetTitleOffset(1)
        h_copy.GetXaxis().SetTitleSize(19)
        h_copy.GetXaxis().SetLabelFont(43)
        h_copy.GetXaxis().SetLabelOffset(0.009)
        h_copy.GetXaxis().SetLabelSize(18)
        h_copy.GetYaxis().SetTitleFont(43)
        h_copy.GetYaxis().SetTitleOffset(1.5)
        h_copy.GetYaxis().SetTitleSize(19)
        h_copy.GetYaxis().SetLabelFont(43)
        h_copy.GetYaxis().SetLabelOffset(0.009)
        h_copy.GetYaxis().SetLabelSize(18)
        h_copy.Scale(1. / h_copy.Integral(), "width")
        h_copy.SetMarkerStyle(marker)
        h_copy.SetMarkerColor(color)
        h_copy.SetLineColor(color)
        h_copy.SetMarkerSize(1.2)
        htitle = ROOT.TPaveText(0.12, 0.99, 0.95, 0.93, "NB NDC")
        htitle.SetBorderSize(0)
        htitle.SetFillStyle(0)
        htitle.SetTextFont(43)
        htitle.SetTextSize(18)
        htitle.AddText("{0:.0f} < #it{{p}}_{{T,ch jet}} < {1:.0f} GeV/#it{{c}}".format(ptmin, ptmax))
        htitle.Draw()
        globalList.append(htitle)
    return histos

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare B feed-down correction file.')
    parser.add_argument('--data', metavar='conf.yaml')
    parser.add_argument('--mc', metavar='conf.yaml')

    args = parser.parse_args()

    f = open(args.data, 'r')
    data_config = yaml.load(f)
    f.close()

    f = open(args.mc, 'r')
    mc_config = yaml.load(f)
    f.close()

    main(data_config, mc_config)

    IPython.embed()
