#!/usr/bin/env python
#python utilities for D Meson jet analysis

import os
import ROOT
import math
import array

def ConvertDMesonName(dmeson):
    if "D0" in dmeson:
        return "D^{0} #rightarrow K^{-}#pi^{+} and c.c."
    else:
        return dmeson

def binom(n, k):
    return math.factorial(n) / math.factorial(n-k) / math.factorial(k)

def GenerateMultiCanvas(cname, n):
    rows = int(math.floor(math.sqrt(n)))
    cols = int(math.ceil(float(n) / rows))
    c = ROOT.TCanvas(cname, cname, cols*400, rows*400)
    c.Divide(cols, rows)
    return c

def find_file(path, file_name):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file == file_name:
                yield os.path.join(root, file)

def CompareSpectra(baseline, spectra, comparisonName, opt="", optRatio="", yaxisRatio="ratio", doSpectra="logy", doRatio="lineary", c=None, cRatio=None, leg=None, legRatio=None, styles=None):
    results = []
    baselineRatio = None
    mainRatioHist = None
    maxRatio = 0
    minRatio = 999
    mainHist = None
    
    print("CompareSpectra: {0}".format(comparisonName))
    print(baseline.GetName())
    for s in spectra:
        print(s.GetName())

    if styles:
        colors = styles["colors"]
        markers = styles["markers"]
        lines = styles["lines"]        
    else:
        colors = [ROOT.kBlack, ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kAzure+2, ROOT.kMagenta+2, ROOT.kCyan+2, ROOT.kPink+1, ROOT.kTeal+2]
        markers = [ROOT.kOpenCircle, ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond, ROOT.kFullStar, ROOT.kStar, ROOT.kOpenCircle]
        lines = [1, 2, 9, 5, 7, 10, 4, 3, 6, 8, 9]

    if doSpectra:
        if not c:
            c = ROOT.TCanvas(comparisonName, comparisonName)

        results.append(c)
        c.cd()
        if doSpectra == "logy":
            c.SetLogy()

        if leg:
            leg.SetY1(leg.GetY1()-0.04*(len(spectra)+1))
        else:
            leg = ROOT.TLegend(0.55, 0.87-0.04*(len(spectra)+1), 0.85, 0.87)
            leg.SetName("{0}_legend".format(c.GetName()))
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextFont(43)
            leg.SetTextSize(16)

        results.append(leg)

        if "hist" in opt:
            baseline.SetLineColor(colors[0])
            baseline.SetLineWidth(2)
            baseline.SetLineStyle(lines[0])
            leg.AddEntry(baseline, baseline.GetTitle(), "l")
        else: 
            baseline.SetMarkerColor(colors[0])
            baseline.SetLineColor(colors[0])
            baseline.SetMarkerStyle(markers[0])
            baseline.SetMarkerSize(1.2)
            leg.AddEntry(baseline, baseline.GetTitle(), "pe")

        baseline.Draw(opt)
        if "frac" in baseline.GetYaxis().GetTitle():
            c.SetLeftMargin(0.12)
            baseline.GetYaxis().SetTitleOffset(1.4)
        if not "same" in opt:
            opt += "same"

        for obj in c.GetListOfPrimitives():
            if isinstance(obj, ROOT.TH1):
                mainHist = obj
                mainHist.SetMinimum(-1111)
                mainHist.SetMaximum(-1111)
                print("Main histogram is: {0}".format(mainHist.GetName()))
                break
        minY = min(mainHist.GetMinimum(0), baseline.GetMinimum(0))
        maxY = max(mainHist.GetMaximum(), baseline.GetMaximum())

    if doRatio:
        cname = "{0}_Ratio".format(comparisonName)
        if not cRatio:
            cRatio = ROOT.TCanvas(cname, cname)
        results.append(cRatio)
        cRatio.cd()
        if doRatio == "logy":
            cRatio.SetLogy()

        if legRatio:
            legRatio.SetY1(legRatio.GetY1()-0.04*(len(spectra)))
        else:
            legRatio = ROOT.TLegend(0.55, 0.87-0.04*(len(spectra)+1), 0.85, 0.87)
            legRatio.SetName("{0}_legend".format(cRatio.GetName()))
            legRatio.SetFillStyle(0)
            legRatio.SetBorderSize(0)
            legRatio.SetTextFont(43)
            legRatio.SetTextSize(16)

        results.append(legRatio)

    for color, marker, line, h in zip(colors[1:], markers[1:], lines[1:], spectra):
        if doSpectra:
            c.cd()
            if minY > h.GetMinimum(0):
                minY = h.GetMinimum(0)
            if maxY < h.GetMaximum():
                maxY = h.GetMaximum()
            h.Draw(opt)
            if "hist" in opt:
                h.SetLineColor(color)
                h.SetLineWidth(3)
                h.SetLineStyle(line)
                leg.AddEntry(h, h.GetTitle(), "l")
            else:
                h.SetMarkerColor(color)
                h.SetLineColor(color)
                h.SetMarkerStyle(marker)
                h.SetMarkerSize(1.2)
                leg.AddEntry(h, h.GetTitle(), "pe")

        if doRatio:
            cRatio.cd()
            hRatio = h.Clone("{0}_Ratio".format(h.GetName()))
            hRatio.GetYaxis().SetTitle(yaxisRatio)
            if not baselineRatio:
                baselineRatio = hRatio
            if "hist" in optRatio:
                hRatio.SetLineColor(color)
                hRatio.SetLineWidth(3)
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
            hRatio.Draw(optRatio)
            if not mainRatioHist:
                for obj in cRatio.GetListOfPrimitives():
                    if isinstance(obj, ROOT.TH1):
                        mainRatioHist = obj
                        print("Main ratio histogram is: {0}".format(mainRatioHist.GetName()))
                        mainRatioHist.SetMinimum(-1111)
                        mainRatioHist.SetMaximum(-1111)
                        minRatio = mainRatioHist.GetMinimum(0)
                        maxRatio = mainRatioHist.GetMaximum()
                        break
            if minRatio > hRatio.GetMinimum(0):
                minRatio = hRatio.GetMinimum(0)
            if maxRatio < hRatio.GetMaximum():
                maxRatio = hRatio.GetMaximum()

            if not "same" in optRatio:
                optRatio += "same"

    if doRatio:
        if doRatio == "logy":
            maxRatio *= 10
            minRatio /= 5
        else:
            maxRatio *= 1.5
            if minRatio < 0.2:
                minRatio = 0
            else:
                minRatio *= 0.6
        mainRatioHist.SetMinimum(minRatio)
        mainRatioHist.SetMaximum(maxRatio)
        cRatio.cd()
        legRatio.Draw()

    if doSpectra:
        if doSpectra == "logy":
            maxY *= 10
            minY /= 5
        else:
            maxY *= 1.5
            if minY < 0.2:
                minY = 0
            else:
                minY *= 0.6
        mainHist.SetMinimum(minY)
        mainHist.SetMaximum(maxY)
        c.cd()
        leg.Draw()

    return results

def DivideNoErrors(ratio, den):
    if not ratio.GetNbinsX() == den.GetNbinsX():
        print("DMesonJetUtils.DivideNoErrors: histograms have different number of bins!")
        return False

    for ibin in range(0, ratio.GetNbinsX()+2):
        if den.GetBinContent(ibin) == 0:
            ratio.SetBinContent(ibin, 0)
        else:
            ratio.SetBinContent(ibin, ratio.GetBinContent(ibin)/den.GetBinContent(ibin))

    return True

def V2TH1(vect):
    result = ROOT.TH1D("vect", "vect", len(vect)-2, 1, len(vect)-2)
    for ibin in range(0, result.GetNbinsX()+2):
        result.SetBinContent(ibin, vect[ibin])
    return result

def BuildHistogram(axis, name, yaxis):
    if len(axis) == 1:
        hist = ROOT.TH1D(name, name, len(axis[0].fBins)-1, array.array('d',axis[0].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(yaxis)
        hist.Sumw2()
    elif len(axis) == 2:
        hist = ROOT.TH2D(name, name, len(axis[0].fBins)-1, array.array('d',axis[0].fBins), len(axis[1].fBins)-1, array.array('d',axis[1].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(axis[1].GetTitle())
        hist.GetZaxis().SetTitle(yaxis)
        hist.Sumw2()
    else:
        hist = ROOT.TH2D(name, name, len(axis[0].fBins)-1, array.array('d',axis[0].fBins), len(axis[1].fBins)-1, array.array('d',axis[1].fBins), len(axis[2].fBins)-1, array.array('d',axis[2].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(axis[1].GetTitle())
        hist.GetZaxis().SetTitle(axis[2].GetTitle())
        hist.Sumw2()
    return hist

def CompareAxis(axis1, axis2):
    if axis1.GetNbins() != axis2.GetNbins(): return False
    for ibin in range(0,axis1.GetNbins()+2):
        if axis1.GetLowEdge(ibin) != axis2.GetLowEdge(ibin): return False
    return True 
