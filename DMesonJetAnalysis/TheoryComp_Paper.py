#!/usr/bin/env python
# python script to compare D0 jet spectrum with theory
# use it with
# JetPtSpectrum_ComparePYTHIA6vs8.yaml
# JetPtSpectrum_CompareHerwig.yaml
# JetPtSpectrum_TheoryComp_Paper.yaml
# JetPtSpectrum_PowhegComp_Paper.yaml
# JetPtSpectrum_FullChargedComp_Paper.yaml
# JetZSpectrum_JetPt_5_15_TheoryComp_Paper.yaml
# JetZSpectrum_JetPt_5_15_PowhegComp_Paper.yaml
# JetZSpectrum_JetPt_5_15_FullChargedComp_Paper.yaml
# JetZSpectrum_JetPt_15_30_TheoryComp_Paper.yaml
# JetZSpectrum_JetPt_15_30_PowhegComp_Paper.yaml
# JetZSpectrum_JetPt_15_30_FullChargedComp_Paper.yaml

import argparse
import yaml
import IPython
import ROOT
import DMesonJetUtils
import LoadTheoryCrossSections

globalList = []

def GetMeasuredCrossSection(input_path, file_name):
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

def GetTheoryCrossSectionAll(config, axis):
    return LoadTheoryCrossSections.GetD0JetTheoryCrossSectionAll(config, axis)

def PlotCrossSections(dataStat, dataSyst, config):
    cname = config["name"]
    if "canvas_h" in config:
        canvas_h = config["canvas_h"]
    else:
        canvas_h = 900
    if "canvas_w" in config:
        canvas_w = config["canvas_w"]
    else:
        canvas_w = 700
    canvas = ROOT.TCanvas(cname, cname, canvas_w, canvas_h)
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
    if config["logy"]: padMain.SetLogy()
    h = dataStat.DrawCopy("axis")
    h.GetYaxis().SetRangeUser(config["miny"], config["maxy"])
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

    for t in config["theory"]:
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
    hRatio.GetYaxis().SetRangeUser(config["minr"], config["maxr"])
    hRatio.GetYaxis().SetNdivisions(509)

    ratioSyst = dataSyst_copy.Clone("ratioSyst")
    globalList.append(ratioSyst)
    ratioSyst.Draw("2")

    ratioStat = dataStat_copy.Clone("ratioStat")
    globalList.append(ratioStat)
    ratioStat.Draw("same p e0 x0")
    ratioStat.SetFillStyle(0)

    for ibin in range(0, ratioSyst.GetN()):
        ratioSyst.SetPointEYlow(ibin, ratioSyst.GetErrorYlow(ibin) / ratioSyst.GetY()[ibin])
        ratioSyst.SetPointEYhigh(ibin , ratioSyst.GetErrorYhigh(ibin) / ratioSyst.GetY()[ibin])
        ratioSyst.SetPoint(ibin, ratioSyst.GetX()[ibin], 1.0)
        ratioStat.SetBinError(ibin + 1, ratioStat.GetBinError(ibin + 1) / ratioStat.GetBinContent(ibin + 1))
        ratioStat.SetBinContent(ibin + 1, 1.0)

    for t in config["theory"]:
        if not t["active"]: continue
        r = t["histogram_plot"].Clone()
        for ibin in range(1, r.GetNbinsX() + 1):
            if t["histogram_plot"].GetBinContent(ibin) <= 0: 
                r.SetBinError(ibin, 0)
                r.SetBinContent(ibin, 0)
                continue
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
    if len(config["title"]) > 3: y2 = 0.66 - 0.06 * (len(config["title"]) - 3)
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
    for t in config["theory"]: 
        if not t["active"]: continue
        leg1.AddEntry(t["histogram_plot"], t["title"], "l")
    leg1.Draw()

    leg1 = ROOT.TLegend(config["legx"], y1 - 0.12, config["legx"] + 0.30, y1, "", "NB NDC")
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
    paveALICE = ROOT.TPaveText(0.16, 0.90 - 0.07 * len(config["title"]), 0.55, 0.90, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(22)
    paveALICE.SetTextAlign(13)
    for line in config["title"]: paveALICE.AddText(line)
    paveALICE.Draw()

    padRatio.RedrawAxis("g")
    padRatio.RedrawAxis()
    padMain.RedrawAxis()

    return canvas

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    dataStat, dataSyst = GetMeasuredCrossSection(config["input_path"], config["data"])
    GetTheoryCrossSectionAll(config, dataStat.GetXaxis())
    canvas = PlotCrossSections(dataStat, dataSyst, config)
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
