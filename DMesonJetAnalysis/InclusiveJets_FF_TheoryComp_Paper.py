#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import argparse
import math
import numpy
import LoadInclusiveJetSpectrum

globalList = []


def GetD0JetCrossSection(input_path, file_name):
    fname = "{}/{}.root".format(input_path, file_name)
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


def GetInclJetCrossSection():
    return LoadInclusiveJetSpectrum.GetCrossSection()


def GetD0JetTheoryCrossSectionAll(config, axis):
    for t in config["theory"]:
        if not t["active"]: continue
        h = GetD0JetTheoryCrossSection(config["input_path"], t["gen"], t["proc"], t["ts"], config["theory_spectrum"], axis)
        t["histogram"] = h


def GetD0JetTheoryCrossSection(input_path, gen, proc, ts, spectrum, axis):
    fname = "{input_path}/FastSim_{gen}_{proc}_{ts}/FastSimAnalysis_ccbar_{gen}_{proc}_{ts}.root".format(input_path=input_path, gen=gen, proc=proc, ts=ts)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    h_orig = DMesonJetUtils.GetObject(file, "D0_MCTruth/Charged_R040/D0_MCTruth_Charged_R040_{spectrum}/D0_MCTruth_Charged_R040_{spectrum}".format(spectrum=spectrum))
    if not h_orig:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)

    h = DMesonJetUtils.Rebin1D(h_orig, axis)
    h.Scale(0.5, "width")  # particle/antiparticle

    return h


def GetInclusiveJetTheoryCrossSectionAll(config):
    for t in config["theory"]:
        if not t["inclusive"]: continue
        if not t["active"]: continue
        h = GetInclusiveJetTheoryCrossSection(config["input_path"], t["inclusive"]["gen"], t["inclusive"]["proc"], t["inclusive"]["ts"])
        t["inclusive_histogram"] = h


def GetInclusiveJetTheoryCrossSection(input_path, gen, proc, ts):
    fname = "{input_path}/FastSim_{gen}_{proc}_{ts}/Jets.root".format(input_path=input_path, gen=gen, proc=proc, ts=ts)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    h_orig = DMesonJetUtils.GetObject(file, "JetPt")
    if not h_orig:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)

    h = h_orig.Clone()
    return h


def PlotCrossSections(dataStat, dataSyst, theory, title, miny, maxy, minr, maxr, legx):
    cname = "D0JetVsInclusiveCrossSectionTheoryComp_Paper"
    canvas = ROOT.TCanvas(cname, cname, 700, 900)
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
    padRatio.SetGridy()
    padRatio.SetTicks(1, 1)

    padMain.cd()
    h = dataStat.DrawCopy("axis")
    h.GetYaxis().SetRangeUser(miny, maxy)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(26)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleOffset(1.6)
    # h.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]")
    # h.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")

    dataSyst_copy = dataSyst.Clone("{0}_copy".format(dataSyst.GetName()))
    dataSyst_copy.Draw("2")
    globalList.append(dataSyst_copy)

    dataStat_copy = dataStat.DrawCopy("same p e0 x0")
    globalList.append(dataStat_copy)

    for t in theory:
        if not t["active"]: continue
        h = t["histogram"].Clone(t["gen"])
        h.Draw("same e0")
        globalList.append(h)
        h.SetLineColor(t["color"])
        h.SetLineStyle(t["line"])
        h.SetLineWidth(2)
        h.SetMarkerStyle(1)
        h.SetMarkerColor(t["color"])
        t["histogram_plot"] = h

    padRatio.cd()

    hRatio = dataStat_copy.DrawCopy("axis")
    hRatio.GetYaxis().SetTitle("theory / data")
    hRatio.GetXaxis().SetTitleFont(43)
    hRatio.GetXaxis().SetTitleSize(26)
    hRatio.GetXaxis().SetLabelFont(43)
    hRatio.GetXaxis().SetLabelSize(22)
    hRatio.GetYaxis().SetTitleFont(43)
    hRatio.GetYaxis().SetTitleSize(26)
    hRatio.GetYaxis().SetLabelFont(43)
    hRatio.GetYaxis().SetLabelSize(22)
    hRatio.GetYaxis().SetTitleOffset(1.4)
    hRatio.GetXaxis().SetTitleOffset(2.9)
    hRatio.GetYaxis().SetRangeUser(minr, maxr)
    hRatio.GetYaxis().SetNdivisions(509)

    ratioSyst = dataSyst_copy.Clone("ratioSyst")
    globalList.append(ratioSyst)
    ratioSyst.Draw("2")

    ratioStat = dataStat_copy.Clone("ratioStat")
    globalList.append(ratioStat)
    ratioStat.Draw("same p e0 x0")
    # ratioStat.Draw("same e2")
    ratioStat.SetFillStyle(0)
    # ratioStat.SetMarkerSize(0)

    for ibin in range(0, ratioSyst.GetN()):
        ratioSyst.SetPointEYlow(ibin, ratioSyst.GetErrorYlow(ibin) / ratioSyst.GetY()[ibin])
        ratioSyst.SetPointEYhigh(ibin , ratioSyst.GetErrorYhigh(ibin) / ratioSyst.GetY()[ibin])
        ratioSyst.SetPoint(ibin, ratioSyst.GetX()[ibin], 1.0)
        ratioStat.SetBinError(ibin + 1, ratioStat.GetBinError(ibin + 1) / ratioStat.GetBinContent(ibin + 1))
        ratioStat.SetBinContent(ibin + 1, 1.0)

    for t in theory:
        if not t["active"]: continue
        r = t["histogram_plot"].Clone()
        for ibin in range(1, r.GetNbinsX() + 1):
            r.SetBinError(ibin, t["histogram_plot"].GetBinError(ibin) / t["histogram_plot"].GetBinContent(ibin))
            r.SetBinContent(ibin, t["histogram_plot"].GetBinContent(ibin) / dataStat_copy.GetBinContent(ibin))
        r.SetLineColor(t["color"])
        r.SetLineStyle(t["line"])
        r.SetLineWidth(2)
        r.SetMarkerStyle(20)
        r.SetMarkerSize(0)
        r.Draw("same e0")
        globalList.append(r)
        t["ratio"] = r

    padMain.cd()

    y2 = 0.66
    if len(title) > 3: y2 = 0.66 - 0.06 * (len(title) - 3)
    y1 = y2 - 0.036 * len(t) / 2
    if y1 < 0.1: y1 = 0.1
    leg1 = ROOT.TLegend(0.18, y1, 0.85, y2, "", "NB NDC")
    leg1.SetNColumns(2)
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(19)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    for t in theory:
        if not t["active"]: continue
        leg1.AddEntry(t["histogram_plot"], t["title"], "l")
    leg1.Draw()

    leg1 = ROOT.TLegend(legx, y1 - 0.12, legx + 0.30, y1, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(20)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    leg1.AddEntry(dataStat_copy, "Data", "p")
    leg1.AddEntry(dataSyst_copy, "Syst. Unc. (data)", "f")
    leg1.Draw()

    padMain.cd()
    paveALICE = ROOT.TPaveText(0.16, 0.90 - 0.07 * len(title), 0.55, 0.90, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(22)
    paveALICE.SetTextAlign(13)
    for line in title: paveALICE.AddText(line)
    paveALICE.Draw()

    padRatio.RedrawAxis("g")
    padRatio.RedrawAxis()
    padMain.RedrawAxis()

    return canvas


def NormalizeData(config, d0jet_stat, d0jet_syst, incl_stat, incl_syst):
    xsec_tot, stat_xsec_tot, syst_xsec_tot = GetTotalCrossSection(incl_stat, incl_syst, config["min_jet_pt"], config["max_jet_pt"])
    for ibin in range(1, d0jet_stat.GetNbinsX() + 1):
        y = d0jet_stat.GetBinContent(ibin) / xsec_tot
        stat_err = math.sqrt((d0jet_stat.GetBinError(ibin) / d0jet_stat.GetBinContent(ibin)) ** 2 + (stat_xsec_tot / xsec_tot) ** 2) * y
        syst_err = math.sqrt((d0jet_syst.GetErrorY(ibin - 1) / d0jet_stat.GetBinContent(ibin)) ** 2 + (syst_xsec_tot / xsec_tot) ** 2) * y
        d0jet_stat.SetBinContent(ibin, y)
        d0jet_stat.SetBinError(ibin, stat_err)
        d0jet_syst.SetPointError(ibin - 1, d0jet_syst.GetErrorX(ibin - 1), d0jet_syst.GetErrorX(ibin - 1), syst_err, syst_err)
        d0jet_syst.SetPoint(ibin - 1, d0jet_stat.GetXaxis().GetBinCenter(ibin), y)
    d0jet_stat.GetYaxis().SetTitle("R(#it{p}_{T, ch jet},#it{z}) / #Delta#it{z}")
    d0jet_syst.GetYaxis().SetTitle("R(#it{p}_{T, ch jet},#it{z}) / #Delta#it{z}")
    d0jet_stat.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")
    d0jet_syst.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")


def NormalizeTheory(config):
    for t in config["theory"]:
        if not t["active"]: continue
        if not "inclusive_histogram" in t or not "histogram" in t: continue
        xsec_tot, stat_xsec_tot, _ = GetTotalCrossSection(t["inclusive_histogram"], None, config["min_jet_pt"], config["max_jet_pt"])
        h = t["histogram"]
        for ibin in range(1, h.GetNbinsX() + 1):
            y = h.GetBinContent(ibin) / xsec_tot
            stat_err = math.sqrt((h.GetBinError(ibin) / h.GetBinContent(ibin)) ** 2 + (stat_xsec_tot / xsec_tot) ** 2) * y
            h.SetBinContent(ibin, y)
            h.SetBinError(ibin, stat_err)
        h.GetYaxis().SetTitle("R(#it{p}_{T, ch jet},#it{z}) / #Delta#it{z}")
        h.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")


def GetTotalCrossSection(stat, syst, minpt, maxpt):
    xsec_tot = 0
    stat_xsec_tot2 = 0
    syst_xsec_tot = 0
    for ibin in range(1, stat.GetNbinsX() + 1):
        if stat.GetXaxis().GetBinCenter(ibin) < minpt: continue
        if stat.GetXaxis().GetBinCenter(ibin) > maxpt: break
        binw = stat.GetXaxis().GetBinWidth(ibin)
        xsec_tot += stat.GetBinContent(ibin) * binw
        stat_xsec_tot2 += (stat.GetBinError(ibin) * binw) ** 2
        if syst: syst_xsec_tot += syst.GetErrorY(ibin - 1) * binw  # take the weighted average of the rel unc

    stat_xsec_tot = math.sqrt(stat_xsec_tot2)

    print("The total cross section for '{}' is {:.4f} +/- {:.4f} (stat) +/- {:.4f} (syst)".format(stat.GetName(), xsec_tot, stat_xsec_tot, syst_xsec_tot))

    return xsec_tot, stat_xsec_tot, syst_xsec_tot


def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    d0jet_stat, d0jet_syst = GetD0JetCrossSection(config["input_path"], config["data"])
    incl_stat, incl_syst = GetInclJetCrossSection()
    GetD0JetTheoryCrossSectionAll(config, d0jet_stat.GetXaxis())
    GetInclusiveJetTheoryCrossSectionAll(config)
    NormalizeData(config, d0jet_stat, d0jet_syst, incl_stat, incl_syst)
    NormalizeTheory(config)
    canvas = PlotCrossSections(d0jet_stat, d0jet_syst, config["theory"], config["title"],
                               config["miny"], config["maxy"], config["minr"], config["maxr"], config["legx"])
    canvas.SaveAs("{}/{}.pdf".format(config["input_path"], config["name"]))
    canvas.SaveAs("{}/{}.C".format(config["input_path"], config["name"]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Jet pt spectrum theory comparison.')
    parser.add_argument('yaml', metavar='conf.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config)

    IPython.embed()
