#!/usr/bin/env python
# python utilities for D Meson jet analysis

import os
import ROOT
import math
import array
import DMesonJetCompare

def ConvertDMesonName(dmeson):
    if "D0" in dmeson:
        return "D^{0} #rightarrow K^{-}#pi^{+} and c.c."
    else:
        return dmeson

def binom(n, k):
    return math.factorial(n) / math.factorial(n - k) / math.factorial(k)

def GenerateMultiCanvas(cname, n):
    rows = int(math.floor(math.sqrt(n)))
    cols = int(math.ceil(float(n) / rows))
    c = ROOT.TCanvas(cname, cname, cols * 400, rows * 400)
    c.Divide(cols, rows)
    return c

def find_file(path, file_name):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file == file_name:
                yield os.path.join(root, file)

def FindMinimum(histogram, limit=0.):
    m = None
    for ibin in range(1, histogram.GetNbinsX() + 1):
        if histogram.GetBinContent(ibin) <= limit:
            continue
        if m is None or histogram.GetBinContent(ibin) < m: m = histogram.GetBinContent(ibin)
    return m

def FindMaximum(histogram, limit=0.):
    m = None
    for ibin in range(1, histogram.GetNbinsX() + 1):
        if histogram.GetBinContent(ibin) <= limit: continue
        if m is None or histogram.GetBinContent(ibin) > m: m = histogram.GetBinContent(ibin)
    return m

def CompareSpectra(baseline, spectra, comparisonName, opt="", optRatio="", yaxisRatio="ratio", doSpectra="logy", doRatio="lineary", c=None, cRatio=None, leg=None, legRatio=None, styles=None):
    comp = DMesonJetCompare.DMesonJetCompare(comparisonName)
    comp.fOptSpectrum = opt
    comp.fOptRatio = optRatio
    comp.fYaxisRatio = yaxisRatio
    comp.fDoSpectraPlot = doSpectra
    comp.fDoRatioPlot = doRatio
    comp.fCanvasSpectra = c
    comp.fCanvasRatio = cRatio
    comp.fLegendSpectra = leg
    comp.fLegendRatio = legRatio
    comp.fColors = [ROOT.kBlack, ROOT.kBlue + 2, ROOT.kRed + 2, ROOT.kGreen + 2, ROOT.kOrange + 2, ROOT.kAzure + 2, ROOT.kMagenta + 2, ROOT.kCyan + 2, ROOT.kPink + 1, ROOT.kTeal + 2]
    comp.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond, ROOT.kFullStar, ROOT.kStar, ROOT.kOpenCircle]
    comp.fLines = [1, 2, 9, 5, 7, 10, 4, 3, 6, 8, 9]
    if styles:
        if "colors" in styles: comp.fColors = styles["colors"]
        if "markers" in styles: comp.fMarkers = styles["markers"]
        if "lines" in styles: comp.fLines = styles["lines"]

    if c:
        for obj in c.GetListOfPrimitives():
            if isinstance(obj, ROOT.TH1):
                comp.fMainHistogram = obj
                comp.fMainHistogram.SetMinimum(-1111)
                comp.fMainHistogram.SetMaximum(-1111)
                print("Main histogram is: {0}".format(comp.fMainHistogram.GetName()))
                break

    if cRatio:
        for obj in cRatio.GetListOfPrimitives():
            if isinstance(obj, ROOT.TH1):
                comp.fMainRatioHistogram = obj
                print("Main ratio histogram is: {0}".format(comp.fMainRatioHistogram.GetName()))
                comp.fMainRatioHistogram.SetMinimum(-1111)
                comp.fMainRatioHistogram.SetMaximum(-1111)
                break

    return comp.CompareSpectra(baseline, spectra)

def DivideNoErrors(ratio, den):
    if not ratio.GetNbinsX() == den.GetNbinsX():
        print("DMesonJetUtils.DivideNoErrors: histograms have different number of bins!")
        return False

    for ibin in range(0, ratio.GetNbinsX() + 2):
        if den.GetBinContent(ibin) == 0:
            ratio.SetBinContent(ibin, 0)
        else:
            ratio.SetBinContent(ibin, ratio.GetBinContent(ibin) / den.GetBinContent(ibin))

    return True

def V2TH1(vect):
    result = ROOT.TH1D("vect", "vect", len(vect) - 2, 1, len(vect) - 2)
    for ibin in range(0, result.GetNbinsX() + 2):
        result.SetBinContent(ibin, vect[ibin])
    return result

def BuildHistogram(axis, name, yaxis):
    if len(axis) == 1:
        hist = ROOT.TH1D(name, name, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(yaxis)
        hist.Sumw2()
    elif len(axis) == 2:
        hist = ROOT.TH2D(name, name, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins), len(axis[1].fBins) - 1, array.array('d', axis[1].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(axis[1].GetTitle())
        hist.GetZaxis().SetTitle(yaxis)
        hist.Sumw2()
    else:
        hist = ROOT.TH2D(name, name, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins), len(axis[1].fBins) - 1, array.array('d', axis[1].fBins), len(axis[2].fBins) - 1, array.array('d', axis[2].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(axis[1].GetTitle())
        hist.GetZaxis().SetTitle(axis[2].GetTitle())
        hist.Sumw2()
    return hist

def CompareAxis(axis1, axis2):
    if axis1.GetNbins() < axis2.GetNbins(): return -1
    elif axis1.GetNbins() > axis2.GetNbins(): return 1
    for ibin in range(0, axis1.GetNbins() + 2):
        if axis1.GetBinLowEdge(ibin) > axis2.GetBinLowEdge(ibin): return -1
        elif axis1.GetBinLowEdge(ibin) < axis2.GetBinLowEdge(ibin): return 1
    return 0

def Rebin1D(hist, xaxis, warnings=False):
    return Rebin1D_fromBins(hist, hist.GetName(), xaxis.GetNbins(), xaxis.GetXbins().GetArray(), warnings)

def Rebin1D_fromBins(hist, name, nbinsX, binsX, warnings=False):
    r = ROOT.TH1D(name, name, nbinsX, binsX)
    r.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    r.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    for xbin in range(0, hist.GetXaxis().GetNbins() + 2):
        xbinCenter = hist.GetXaxis().GetBinCenter(xbin)
        rxbin = r.GetXaxis().FindBin(xbinCenter)
        binValue = hist.GetBinContent(xbin) + r.GetBinContent(rxbin)
        binError = math.sqrt(hist.GetBinError(xbin) ** 2 + r.GetBinError(rxbin) ** 2)
        r.SetBinContent(rxbin, binValue)
        r.SetBinError(rxbin, binError)
        if binValue > 0:
            relErr = binError / binValue
            if relErr > 0.9 and warnings:
                print("Bin ({0}) has rel stat err = {1}. This is VERY dangerous!".format(xbin, relErr))
    return r

def Rebin2D(hist, xaxis, yaxis, warnings=False):
    return Rebin2D_fromBins(hist, hist.GetName(), xaxis.GetNbins(), xaxis.GetXbins().GetArray(), yaxis.GetNbins(), yaxis.GetXbins().GetArray(), warnings)

def Rebin2D_fromBins(hist, name, nbinsX, binsX, nbinsY, binsY, warnings=False):
    r = ROOT.TH2D(name, name, nbinsX, binsX, nbinsY, binsY)
    r.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    r.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    for xbin in range(0, hist.GetXaxis().GetNbins() + 2):
        for ybin in range(0, hist.GetYaxis().GetNbins() + 2):
            xbinCenter = hist.GetXaxis().GetBinCenter(xbin)
            ybinCenter = hist.GetYaxis().GetBinCenter(ybin)
            rxbin = r.GetXaxis().FindBin(xbinCenter)
            rybin = r.GetYaxis().FindBin(ybinCenter)
            binValue = hist.GetBinContent(xbin, ybin) + r.GetBinContent(rxbin, rybin)
            binError = math.sqrt(hist.GetBinError(xbin, ybin) ** 2 + r.GetBinError(rxbin, rybin) ** 2)
            r.SetBinContent(rxbin, rybin, binValue)
            r.SetBinError(rxbin, rybin, binError)

    for xbin in range(0, r.GetXaxis().GetNbins() + 2):
        for ybin in range(0, r.GetYaxis().GetNbins() + 2):
            binValue = r.GetBinContent(xbin, ybin)
            binError = r.GetBinError(xbin, ybin)
            if binValue > 0:
                relErr = binError / binValue
                if relErr > 0.9 and warnings:
                    print("Bin ({0},{1}) has rel stat err = {2}. This is VERY dangerous!".format(xbin, ybin, relErr))
    return r

