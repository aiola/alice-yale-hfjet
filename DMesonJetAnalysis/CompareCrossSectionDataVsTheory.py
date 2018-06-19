#!/usr/bin/env python
# python script to compare D0 jet spectrum with theory
# use it with
# JetPtSpectrum_ComparePYTHIA6vs8.yaml
# JetPtSpectrum_CompareHerwig.yaml
# JetPtSpectrum_ComparePowheg.yaml
# JetPtSpectrum_CompareTheory.yaml
# JetPtSpectrum_CompareFullVsCharged.yaml
# JetPtSpectrum_ComparePtDCut.yaml

# JetZSpectrum_JetPt_5_15_ComparePYTHIA6vs8.yaml
# JetZSpectrum_JetPt_5_15_CompareHerwig.yaml
# JetZSpectrum_JetPt_5_15_ComparePowheg.yaml
# JetZSpectrum_JetPt_5_15_CompareTheory.yaml
# JetZSpectrum_JetPt_5_15_CompareFullVsCharged.yaml
# JetZSpectrum_JetPt_5_15_ComparePtDCut.yaml

# JetZSpectrum_JetPt_15_30_ComparePYTHIA6vs8.yaml
# JetZSpectrum_JetPt_15_30_CompareHerwig.yaml
# JetZSpectrum_JetPt_15_30_ComparePowheg.yaml
# JetZSpectrum_JetPt_15_30_CompareTheory.yaml
# JetZSpectrum_JetPt_15_30_CompareFullVsCharged.yaml
# JetZSpectrum_JetPt_15_30_ComparePtDCut.yaml

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
    cname = "{}_{}".format(config["name_prefix"], config["name"])
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
    if config["logy"]:
        padMain.SetLogy()
    
    if "minx" in config:
        minx = config["minx"]
    else:
        minx = dataStat.GetXaxis().GetXmin()

    if "maxx" in config:
        maxx = config["maxx"]
    else:
        maxx = dataStat.GetXaxis().GetXmax()

    h = ROOT.TH1I("h", "h", 1000, minx, maxx)
    h.Draw("axis")
    globalList.append(h)
    h.GetYaxis().SetRangeUser(config["miny"], config["maxy"])
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(26)
    h.GetYaxis().SetTitleOffset(1.8)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    if "y_axis_title" in config:
        h.GetYaxis().SetTitle(config["y_axis_title"])
    else:
        h.GetYaxis().SetTitle(dataStat.GetYaxis().GetTitle())

    dataSyst_copy = dataSyst.Clone("{0}_copy".format(dataSyst.GetName()))
    dataSyst_copy.Draw("2")
    globalList.append(dataSyst_copy)

    dataStat_copy = dataStat.DrawCopy("same p e0 x0")
    globalList.append(dataStat_copy)

    if "data_minx" in config and "data_maxx" in config:
        points_to_be_removed = set()
        for ipoint in range(0, dataSyst_copy.GetN()):
            if dataSyst_copy.GetX()[ipoint] < config["data_minx"] or dataSyst_copy.GetX()[ipoint] > config["data_maxx"]:
                points_to_be_removed.add(ipoint)
        for ipoint in points_to_be_removed:
            dataSyst_copy.RemovePoint(ipoint)
        dataStat_copy.GetXaxis().SetRangeUser(config["data_minx"], config["data_maxx"])

    for t in config["theory"]:
        if not t["active"]: continue
        h = t["histogram"].Clone(t["gen"])
        h.Draw("same e0")
        if "data_minx" in config and "data_maxx" in config:
            h.GetXaxis().SetRangeUser(config["data_minx"], config["data_maxx"])
        globalList.append(h)
        h.SetLineColor(t["color"])
        h.SetLineStyle(t["line"])
        h.SetLineWidth(2)
        h.SetMarkerStyle(1)
        h.SetMarkerColor(t["color"])
        t["histogram_plot"] = h

    padRatio.cd()

    hRatio = ROOT.TH1I("h", "h", 1000, minx, maxx)
    hRatio.Draw("axis")
    globalList.append(hRatio)
    hRatio.GetYaxis().SetTitle("MC / Data")
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
    if "x_axis_title" in config:
        hRatio.GetXaxis().SetTitle(config["x_axis_title"])
    else:
        hRatio.GetXaxis().SetTitle(dataStat.GetXaxis().GetTitle())

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
    
    for ibin in range(1, ratioStat.GetNbinsX()+1):
        ratioStat.SetBinError(ibin, ratioStat.GetBinError(ibin) / ratioStat.GetBinContent(ibin))
        ratioStat.SetBinContent(ibin, 1.0)

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

    y1 = 0.90
    y2 = y1 - 0.07 * len(config["title"])
    padMain.cd()
    paveALICE = ROOT.TPaveText(0.16, y1, 0.55, y2, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(22)
    paveALICE.SetTextAlign(13)
    for line in config["title"]: paveALICE.AddText(line)
    paveALICE.Draw()

    if "theory_legend" in config and "n_columns" in config["theory_legend"]:
        n_leg_columns = config["theory_legend"]["n_columns"]
    else:
        max_length = 0
        for t in config["theory"]: 
            if not t["active"]:
                continue
            if len(t["title"]) > max_length:
                max_length = len(t["title"])

        if max_length > 40:
            n_leg_columns = 1
        else:
            n_leg_columns = 2

    if "theory_legend" in config and "y" in config["theory_legend"]:
        y1 = config["theory_legend"]["y"]
    else:
        y1 = y2 - 0.03
    active_t = len([t for t in config["theory"] if "histogram_plot" in t])
    y2 = y1 - 0.06 * active_t / n_leg_columns
    if "theory_legend" in config and "x" in config["theory_legend"]:
        x1 = config["theory_legend"]["x"]
    else:
        x1 = 0.16
    x2 = 0.85

    leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetNColumns(n_leg_columns)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(19)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.15)
    for t in config["theory"]: 
        if not t["active"]:
            continue
        leg1.AddEntry(t["histogram_plot"], t["title"], "l")

    leg1.Draw()

    if "data_legend" in config and "y" in config["data_legend"]:
        y1 = config["data_legend"]["y"]
    else:
        y1 = y2 - 0.02
    y2 = y1 - 0.12
    if "data_legend" in config and "x" in config["data_legend"]:
        x1 = config["data_legend"]["x"]
    else:
        x1 = 0.16
    x2 = x1 + 0.30
    leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
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
    canvas.SaveAs("{}/{}.pdf".format(config["input_path"], canvas.GetName()))
    canvas.SaveAs("{}/{}.C".format(config["input_path"], canvas.GetName()))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Jet pt spectrum theory comparison.')
    parser.add_argument('yaml', metavar='conf.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig)

    IPython.embed()
