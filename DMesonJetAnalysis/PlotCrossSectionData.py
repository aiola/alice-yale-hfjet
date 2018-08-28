#!/usr/bin/env python
# python script to fit D0 jet spectrum
# use it with
# JetPtCrossSection_Fit_Thesis.yaml

import argparse
import yaml
import IPython
import ROOT
import DMesonJetUtils

globalList = []

if ROOT.gROOT.GetVersionInt() >= 60000: ROOT.ROOT.Math.IntegratorOneDimOptions.SetDefaultIntegrator("Gauss")

def GetMeasuredCrossSection(input_path, file_name, scale):
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
    hStat.Scale(scale)
    for ipoint in range(0, hSyst.GetN()):
        hSyst.SetPointEYlow(ipoint, hSyst.GetErrorYlow(ipoint) * scale)
        hSyst.SetPointEYhigh(ipoint, hSyst.GetErrorYhigh(ipoint) * scale)
        hSyst.SetPoint(ipoint, hSyst.GetX()[ipoint], hSyst.GetY()[ipoint] * scale)
    return hStat, hSyst

def PlotCrossSections(dataStat, dataSyst, config):
    no_box_sys = ["[]", "||", "0[]", "0||"]
    if not "data_systematics_style" in config:
        config["data_systematics_style"] = "2"
    config["data_box_systematics"] = not (config["data_systematics_style"] in no_box_sys)

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

    if len(config["functions"]) > 0:
        canvas.Divide(1, 2)
        padMain = canvas.cd(1)
        padMain.SetPad(0, 0.35, 1, 1)
        padMain.SetBottomMargin(0)
        padRatio = canvas.cd(2)
        padRatio.SetPad(0, 0., 1, 0.35)
        padRatio.SetTopMargin(0)
        padRatio.SetBottomMargin(config["bottom_margin"])
        padRatio.SetLeftMargin(config["left_margin"])
        padRatio.SetRightMargin(0.05)
        padRatio.SetGridy()
        padRatio.SetTicks(1, 1)
    else:
        padMain = canvas
        padRatio = None
        padMain.SetBottomMargin(config["bottom_margin"])

    padMain.SetLeftMargin(config["left_margin"])
    padMain.SetRightMargin(0.05)
    padMain.SetTicks(1, 1)

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
    h.GetYaxis().SetTitleSize(config["axis"]["font_size"])
    h.GetYaxis().SetTitleOffset(config["axis"]["y_offset"])
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(config["axis"]["font_size"] - 3)
    h.GetYaxis().SetTitle(config["axis"]["y_title"])

    if not padRatio:
        h.GetXaxis().SetTitle(config["axis"]["x_title"])
        h.GetXaxis().SetTitleFont(43)
        h.GetXaxis().SetTitleSize(config["axis"]["font_size"])
        h.GetXaxis().SetLabelFont(43)
        h.GetXaxis().SetLabelSize(config["axis"]["font_size"] - 3)
        h.GetXaxis().SetTitleOffset(config["axis"]["x_offset"])

    dataSyst_copy = dataSyst.Clone("{0}_copy".format(dataSyst.GetName()))
    globalList.append(dataSyst_copy)
    dataStat_copy = dataStat.Clone("{0}_copy".format(dataStat.GetName()))
    globalList.append(dataStat_copy)

    if "data_minx" in config and "data_maxx" in config:
        points_to_be_removed = set()
        for ipoint in range(0, dataSyst_copy.GetN()):
            if dataSyst_copy.GetX()[ipoint] < config["data_minx"] or dataSyst_copy.GetX()[ipoint] > config["data_maxx"]:
                points_to_be_removed.add(ipoint)
        for ipoint in points_to_be_removed:
            dataSyst_copy.RemovePoint(ipoint)
        dataStat_copy.GetXaxis().SetRangeUser(config["data_minx"], config["data_maxx"])

    # Plot data systematic uncertainties, if they are plotted as a filled box
    if config["data_box_systematics"]:
        dataSyst_copy.Draw(config["data_systematics_style"])

    # Plot data statistical uncertainties
    dataStat_copy.Draw("same p e0 x0")

    for f in config["functions"]:
        if not f["active"]:
            continue
        function = ROOT.TF1(f["name"], f["formula"], dataStat_copy.GetXaxis().GetBinLowEdge(1), dataStat_copy.GetXaxis().GetBinUpEdge(dataStat_copy.GetXaxis().GetNbins()))
        for ipar, (pname, pvalue, pfixvalue) in enumerate(zip(f["parameter_names"], f["parameter_init_values"], f["parameter_fixed_values"])):
            function.SetParName(ipar, pname)
            if not pfixvalue is None:
                function.FixParameter(ipar, pfixvalue)
            elif not pvalue is None:
                function.SetParameter(ipar, pvalue)
        fit_results = dataStat_copy.Fit(function, "N S R")
        fit_status = int(fit_results)
        if fit_status != 0:
            print("WARNING: fit status != 0")
        f["function"] = function
        f["fit_results"] = fit_results
        print("Chi2 = {:.4f}".format(fit_results.Chi2()))
        print("Ndf = {:.4f}".format(fit_results.Ndf()))
        print("Chi2 / Ndf = {:.4f}".format(fit_results.Chi2() / fit_results.Ndf()))
        function.SetLineColor(f["color"])
        function.SetLineStyle(f["line"])
        function.SetLineWidth(2)
        function.Draw("same")
        
    # Plot data systematic uncertainties, if they are not plotted as filled boxes
    if not config["data_box_systematics"]:
        for ipoint in range(0, dataSyst_copy.GetN()):
            dataSyst_copy.SetPointEXhigh(ipoint, 0)
            dataSyst_copy.SetPointEXlow(ipoint, 0)
        dataSyst_copy.SetLineColor(dataStat_copy.GetLineColor())
        dataSyst_copy.Draw(config["data_systematics_style"])

    # RATIO

    ratioSyst = dataSyst_copy.Clone("ratioSyst")
    globalList.append(ratioSyst)

    ratioStat = dataStat_copy.Clone("ratioStat")
    ratioStat.SetFillStyle(0)
    globalList.append(ratioStat)

    for ibin in range(0, ratioSyst.GetN()):
        ratioSyst.SetPointEYlow(ibin, ratioSyst.GetErrorYlow(ibin) / ratioSyst.GetY()[ibin])
        ratioSyst.SetPointEYhigh(ibin , ratioSyst.GetErrorYhigh(ibin) / ratioSyst.GetY()[ibin])
        ratioSyst.SetPoint(ibin, ratioSyst.GetX()[ibin], 1.0)
    
    for ibin in range(1, ratioStat.GetNbinsX()+1):
        if ratioStat.GetBinContent(ibin) == 0:
            continue
        ratioStat.SetBinError(ibin, ratioStat.GetBinError(ibin) / ratioStat.GetBinContent(ibin))
        ratioStat.SetBinContent(ibin, 1.0)

    for f in config["functions"]:
        if not f["active"]:
            continue
        h = ROOT.TH1D("{}_hist".format(f["name"]), f["title"], ratioStat.GetNbinsX(), ratioStat.GetXaxis().GetXbins().GetArray())
        h.SetLineColor(f["color"])
        h.SetLineStyle(f["line"])
        h.SetLineWidth(2)
        h.Add(f["function"])
        h.Divide(dataStat_copy)
        f["ratio"] = h

    if padRatio:
        padRatio.cd()

        hRatio = ROOT.TH1I("h", "h", 1000, minx, maxx)
        hRatio.Draw("axis")
        globalList.append(hRatio)
        
        hRatio.GetXaxis().SetTitle(config["axis"]["x_title"])
        hRatio.GetXaxis().SetTitleFont(43)
        hRatio.GetXaxis().SetTitleSize(config["axis"]["font_size"])
        hRatio.GetXaxis().SetLabelFont(43)
        hRatio.GetXaxis().SetLabelSize(config["axis"]["font_size"] - 3)
        hRatio.GetXaxis().SetTitleOffset(config["axis"]["x_offset"])

        hRatio.GetYaxis().SetTitle("Fit / Data")
        hRatio.GetYaxis().SetTitleFont(43)
        hRatio.GetYaxis().SetTitleSize(config["axis"]["font_size"])
        hRatio.GetYaxis().SetLabelFont(43)
        hRatio.GetYaxis().SetLabelSize(config["axis"]["font_size"] - 3)
        hRatio.GetYaxis().SetTitleOffset(config["axis"]["y_offset"])
        hRatio.GetYaxis().SetRangeUser(config["minr"], config["maxr"])
        hRatio.GetYaxis().SetNdivisions(509)
        
        # Plot data systematic uncertainties, if they are plotted as a filled box
        if config["data_box_systematics"]:
            ratioSyst.Draw(config["data_systematics_style"])

        # Plot data statistical uncertainties
        ratioStat.Draw("same p e0 x0")

        for f in config["functions"]:
            if not f["active"]:
                continue
            f["ratio"].Draw("same hist")

        # Plot data systematic uncertainties, if they are not plotted as filled boxes
        if not config["data_box_systematics"]:
            for ipoint in range(0, ratioSyst.GetN()):
                ratioSyst.SetPointEXhigh(ipoint, 0)
                ratioSyst.SetPointEXlow(ipoint, 0)
            ratioSyst.SetLineColor(ratioStat.GetLineColor())
            ratioSyst.Draw(config["data_systematics_style"])

    padMain.cd()

    if "y" in config["title"]:
        y1 = config["title"]["y"]
    else:
        y1 = 0.90
    if "x" in config["title"]:
        x1 = config["title"]["x"]
    else:
        x1 = 0.19
    y2 = y1 - 0.07 * len(config["title"]["text"])
    x2 = x1 + 0.36
    if x2 > 0.99:
        x2 = 0.99
    padMain.cd()
    paveALICE = ROOT.TPaveText(x1, y1, x2, y2, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(config["title"]["font_size"])
    paveALICE.SetTextAlign(13)
    for line in config["title"]["text"]: 
        paveALICE.AddText(line)
    paveALICE.Draw()

    if len(config["functions"]) > 0:
        if "theory_legend" in config and "n_columns" in config["theory_legend"]:
            n_leg_columns = config["theory_legend"]["n_columns"]
        else:
            max_length = 0
            for f in config["functions"]: 
                if not "function" in f:
                    continue
                if len(f["title"]) > max_length:
                    max_length = len(f["title"])

            if max_length > 40:
                n_leg_columns = 1
            else:
                n_leg_columns = 2

        if "theory_legend" in config and "y" in config["theory_legend"]:
            y1 = config["theory_legend"]["y"]
        else:
            y1 = y2 - 0.03
        active_t = len([f for f in config["functions"] if "function" in f])
        y2 = y1 - 0.07 * active_t / n_leg_columns
        if "theory_legend" in config and "x" in config["theory_legend"]:
            x1 = config["theory_legend"]["x"]
        else:
            x1 = 0.16
        x2 = x1 + 0.50
        if x2 > 0.95:
            x2 = 0.95

        leg1 = ROOT.TLegend(x1, y1, x2, y2, "", "NB NDC")
        globalList.append(leg1)
        leg1.SetNColumns(n_leg_columns)
        leg1.SetBorderSize(0)
        leg1.SetFillStyle(0)
        leg1.SetTextFont(43)
        leg1.SetTextSize(config["theory_legend"]["font_size"])
        leg1.SetTextAlign(12)
        leg1.SetMargin(0.1)
        for f in config["functions"]: 
            if not "function" in f:
                continue
            leg1.AddEntry(f["function"], f["title"], "l")
        leg1.Draw()

    if "data_legend" in config and "y" in config["data_legend"]:
        y1 = config["data_legend"]["y"]
    else:
        y1 = y2 - 0.02
    y2 = y1 - 0.06
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
    leg1.SetTextSize(config["data_legend"]["font_size"])
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)

    if config["data_box_systematics"]:
        entry = leg1.AddEntry(None, "Data", "pf")
        entry.SetFillStyle(dataSyst_copy.GetFillStyle())
        entry.SetFillColor(dataSyst_copy.GetFillColor())
        entry.SetLineColor(dataSyst_copy.GetFillColor())
        entry.SetMarkerColor(dataStat_copy.GetMarkerColor())
        entry.SetMarkerStyle(dataStat_copy.GetMarkerStyle())
    else:
        entry = leg1.AddEntry(None, "Data", "p")
        entry.SetMarkerColor(dataStat_copy.GetMarkerColor())
        entry.SetMarkerStyle(dataStat_copy.GetMarkerStyle())

    leg1.Draw()

    if padRatio:
        padRatio.RedrawAxis("g")
        padRatio.RedrawAxis()
    padMain.RedrawAxis()

    return canvas

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if "scale" in config:
        scale = config["scale"]
    else:
        scale = 1.0
    dataStat, dataSyst = GetMeasuredCrossSection(config["input_path"], config["data"], scale)
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
