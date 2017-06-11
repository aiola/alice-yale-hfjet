#!/usr/local/bin/python
# python utilities for D Meson jet analysis

import os
import ROOT
import math
import array
from collections import OrderedDict

def soft_clone(origin, name, title=None, yaxisTitle=None):
    if not title: title = name
    if not yaxisTitle: yaxisTitle = origin.GetYaxis().GetTitle()
    h = ROOT.TH1D(name, title, origin.GetNbinsX(), origin.GetXaxis().GetXbins().GetArray())
    h.GetXaxis().SetTitle(origin.GetXaxis().GetTitle())
    h.GetYaxis().SetTitle(yaxisTitle)
    return h

def GetRelativeUncertaintyHistogram(h):
    h_unc = soft_clone(h, "{0}_unc".format(h.GetName()), "{0} Rel. Unc.".format(h.GetTitle()), "rel. unc.")
    for ibin in xrange(0, h.GetNbinsX() + 2):
        if h.GetBinContent(ibin) == 0: continue
        h_unc.SetBinContent(ibin, h.GetBinError(ibin) / h.GetBinContent(ibin))
    return h_unc

def ConvertDMesonName(dmeson):
    if "D0" in dmeson:
        return "D^{0} #rightarrow K^{-}#pi^{+} and c.c."
    else:
        return dmeson

def binom(n, k):
    return math.factorial(n) / math.factorial(n - k) / math.factorial(k)

def GetObject(obj, name):
    slash = 0
    while(slash >= 0):
        slash = name.find("/", slash)
        if name[slash + 1] == '/':
            slash += 2
        else:
            break

    if slash < 0:
        name_lookup = name.replace("//", "/")
        name = None
    else:
        name_lookup = name[:slash].replace("//", "/")
        name = name[slash + 1:]
    if isinstance(obj, ROOT.TCollection):
        res = obj.FindObject(name_lookup)
        name_obj = obj.GetName()
    elif isinstance(obj, dict) or isinstance(obj, OrderedDict):
        res = obj[name_lookup]
        name_obj = obj
    elif isinstance(obj, ROOT.TDirectory):
        res = obj.Get(name_lookup)
        name_obj = obj.GetName()
    if not res:
        print("Could not find object {0} in collection '{1}'".format(name_lookup, name_obj))
        if isinstance(obj, ROOT.TObject): obj.ls()
        return None
    if name:
        return GetObject(res, name)
    else:
        return res

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

def FindMinimum(histogram, limit=0., errors=True):
    m = None
    if histogram.GetDimension() == 1:
        for ibin in xrange(1, histogram.GetNbinsX() + 1):
            if errors:
                cont = histogram.GetBinContent(ibin) - histogram.GetBinError(ibin)
            else:
                cont = histogram.GetBinContent(ibin)
            if cont <= limit: continue
            if m is None or cont < m: m = cont
    elif histogram.GetDimension() == 2:
        for xbin in xrange(1, histogram.GetNbinsX() + 1):
            for ybin in xrange(1, histogram.GetNbinsY() + 1):
                if errors:
                    cont = histogram.GetBinContent(xbin, ybin) - histogram.GetBinError(xbin, ybin)
                else:
                    cont = histogram.GetBinContent(xbin, ybin)
                if cont <= limit: continue
                if m is None or cont < m: m = cont
    return m

def FindMaximum(histogram, limit=0., errors=True):
    m = None
    if histogram.GetDimension() == 1:
        for ibin in xrange(1, histogram.GetNbinsX() + 1):
            if errors:
                cont = histogram.GetBinContent(ibin) + histogram.GetBinError(ibin)
            else:
                cont = histogram.GetBinContent(ibin)
            if cont <= limit: continue
            if m is None or cont > m: m = cont
    elif histogram.GetDimension() == 2:
        for xbin in xrange(1, histogram.GetNbinsX() + 1):
            for ybin in xrange(1, histogram.GetNbinsY() + 1):
                if errors:
                    cont = histogram.GetBinContent(xbin, ybin) + histogram.GetBinError(xbin, ybin)
                else:
                    cont = histogram.GetBinContent(xbin, ybin)
                if cont <= limit: continue
                if m is None or cont > m: m = cont
    return m

def DivideNoErrors(ratio, den):
    if not ratio.GetNbinsX() == den.GetNbinsX():
        print("DMesonJetUtils.DivideNoErrors: histograms have different number of bins!")
        return False

    for ibin in xrange(0, ratio.GetNbinsX() + 2):
        if den.GetBinContent(ibin) == 0:
            ratio.SetBinContent(ibin, 0)
        else:
            ratio.SetBinContent(ibin, ratio.GetBinContent(ibin) / den.GetBinContent(ibin))

    return True

def V2TH1(vect):
    result = ROOT.TH1D("vect", "vect", len(vect) - 2, 1, len(vect) - 2)
    for ibin in xrange(0, result.GetNbinsX() + 2):
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
        hist = ROOT.TH3D(name, name, len(axis[0].fBins) - 1, array.array('d', axis[0].fBins), len(axis[1].fBins) - 1, array.array('d', axis[1].fBins), len(axis[2].fBins) - 1, array.array('d', axis[2].fBins))
        hist.GetXaxis().SetTitle(axis[0].GetTitle())
        hist.GetYaxis().SetTitle(axis[1].GetTitle())
        hist.GetZaxis().SetTitle(axis[2].GetTitle())
        hist.Sumw2()
    return hist

def CompareAxis(axis1, axis2):
    if axis1.GetNbins() < axis2.GetNbins(): return -1
    elif axis1.GetNbins() > axis2.GetNbins(): return 1
    for ibin in xrange(0, axis1.GetNbins() + 2):
        if axis1.GetBinLowEdge(ibin) > axis2.GetBinLowEdge(ibin): return -1
        elif axis1.GetBinLowEdge(ibin) < axis2.GetBinLowEdge(ibin): return 1
    return 0

def Rebin1D(hist, xaxis, warnings=False):
    return Rebin1D_fromBins(hist, hist.GetName(), xaxis.GetNbins(), xaxis.GetXbins().GetArray(), warnings)

def Rebin1D_fromBins(hist, name, nbinsX, binsX, warnings=False):
    r = ROOT.TH1D(name, name, nbinsX, binsX)
    r.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
    r.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
    for xbin in xrange(0, hist.GetXaxis().GetNbins() + 2):
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
    r.GetZaxis().SetTitle(hist.GetZaxis().GetTitle())
    for xbin in xrange(0, hist.GetXaxis().GetNbins() + 2):
        for ybin in xrange(0, hist.GetYaxis().GetNbins() + 2):
            xbinCenter = hist.GetXaxis().GetBinCenter(xbin)
            ybinCenter = hist.GetYaxis().GetBinCenter(ybin)
            rxbin = r.GetXaxis().FindBin(xbinCenter)
            rybin = r.GetYaxis().FindBin(ybinCenter)
            binValue = hist.GetBinContent(xbin, ybin) + r.GetBinContent(rxbin, rybin)
            binError = math.sqrt(hist.GetBinError(xbin, ybin) ** 2 + r.GetBinError(rxbin, rybin) ** 2)
            r.SetBinContent(rxbin, rybin, binValue)
            r.SetBinError(rxbin, rybin, binError)

    for xbin in xrange(0, r.GetXaxis().GetNbins() + 2):
        for ybin in xrange(0, r.GetYaxis().GetNbins() + 2):
            binValue = r.GetBinContent(xbin, ybin)
            binError = r.GetBinError(xbin, ybin)
            if binValue > 0:
                relErr = binError / binValue
                if relErr > 0.9 and warnings:
                    print("Bin ({0},{1}) has rel stat err = {2}. This is VERY dangerous!".format(xbin, ybin, relErr))
    return r

