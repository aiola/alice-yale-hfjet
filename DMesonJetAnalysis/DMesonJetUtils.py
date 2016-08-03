#!/usr/bin/env python
#python utilities for D Meson jet analysis

import os
import ROOT

def find_file(path, file_name):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file == file_name:
                yield os.path.join(root, file)

def CompareSpectra(baseline, spectra, comparisonName, opt="", optRatio=""):
    results = []
    baselineRatio = None
    maxRatio = 0
    minRatio = 999

    colors = [ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kAzure+2, ROOT.kMagenta+2, ROOT.kTeal+2]
    markers = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond, ROOT.kFullStar]
    lines = range(2,11)
    c = ROOT.TCanvas(comparisonName, comparisonName)
    results.append(c)
    c.cd()
    c.SetLogy()
    if "hist" in opt:
        baseline.SetLineColor(ROOT.kBlue+2)
        baseline.SetLineWidth(2)
        baseline.SetLineStyle(1)
    else: 
        baseline.SetMarkerColor(ROOT.kBlue+2)
        baseline.SetLineColor(ROOT.kBlue+2)
        baseline.SetMarkerStyle(ROOT.kOpenCircle)
        baseline.SetMarkerSize(1.2)

    cname = "{0}_Ratio".format(comparisonName)
    cRatio = ROOT.TCanvas(cname, cname)
    results.append(cRatio)
    cRatio.cd()
    leg = ROOT.TLegend(0.55, 0.87-0.04*(len(spectra)+1), 0.85, 0.87)
    results.append(leg)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(16)
    leg.AddEntry(baseline, baseline.GetTitle(), "pe")

    legRatio = ROOT.TLegend(0.55, 0.87-0.04*len(spectra), 0.85, 0.87)
    results.append(legRatio)
    legRatio.SetFillStyle(0)
    legRatio.SetBorderSize(0)
    legRatio.SetTextFont(43)
    legRatio.SetTextSize(16)

    if not "same" in opt:
        opt += "same"
    if not "same" in optRatio:
        optRatio += "same"
    for color, marker, line, h in zip(colors, markers, lines, spectra):
        c.cd()
        h.Draw(opt)
        if "hist" in opt:
            h.SetLineColor(color)
            h.SetLineWidth(2)
            h.SetLineStyle(line)
            leg.AddEntry(h, h.GetTitle(), "l")
        else:
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            h.SetMarkerStyle(marker)
            h.SetMarkerSize(1.2)
            leg.AddEntry(h, h.GetTitle(), "pe")

        cRatio.cd()
        hRatio = h.Clone("{0}_Ratio".format(h.GetName()))
        hRatio.GetYaxis().SetTitle("ratio")
        if not baselineRatio:
            baselineRatio = hRatio
        if "hist" in optRatio:
            hRatio.SetLineColor(color)
            hRatio.SetLineWidth(2)
            hRatio.SetLineStyle(line)
            legRatio.AddEntry(hRatio, h.GetTitle(), "l")
        else:
            hRatio.SetMarkerColor(color)
            hRatio.SetLineColor(color)
            hRatio.SetMarkerStyle(marker)
            hRatio.SetMarkerSize(1.2)
            legRatio.AddEntry(hRatio, h.GetTitle(), "pe")
        results.append(hRatio)
        hRatio.SetTitle("{0} Ratio".format(h.GetTitle()))
        hRatio.Divide(baseline)
        if minRatio > hRatio.GetMinimum():
            minRatio = hRatio.GetMinimum()
        if maxRatio < hRatio.GetMaximum():
            maxRatio = hRatio.GetMaximum()
        hRatio.Draw(optRatio)
        opt.replace("same","")
        optRatio.replace("same","")
    c.cd()
    baseline.Draw("same")
    maxRatio *= 1.5
    if minRatio < 0.2:
        minRatio = 0
    else:
        minRatio *= 0.6
    baselineRatio.SetMinimum(minRatio)
    baselineRatio.SetMaximum(maxRatio)
    leg.Draw()
    
    cRatio.cd()
    legRatio.Draw()
    return results
