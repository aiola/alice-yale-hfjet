#!/usr/bin/env python
# python script to prepare B feed-down correction file

import argparse
import IPython
import ROOT
import array
import numpy
import math
import yaml
import DMesonJetUnfolding
import DMesonJetUtils
import DMesonJetCompare
import copy
from collections import OrderedDict

globalList = []

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)

    robjects = CompareVariations(config)

    if robjects:
        output_file_name = "{0}/{1}.root".format(config["input_path"], config["name"])
        outputFile = ROOT.TFile(output_file_name, "recreate")
        if not outputFile or outputFile.IsZombie():
            print("Could not open output file {0}".format(output_file_name))
            exit(1)
        for obj in robjects:
            obj.Write(obj.GetName(), ROOT.TObject.kSingleKey)
        outputFile.Close()

    SaveCanvases(config["input_path"])

    print("Done")

def CompareVariations(config):
    files = OpenFiles(config)

    variations = []
    h = DMesonJetUtils.GetObject(files[config["default"]["input_name"]], config["default"]["histogram_name"])
    if not h: exit(1)
    baseline = h.Clone("default")
    baseline.SetTitle(config["default"]["title"])
    for s in config["sources"]:
        if "histogram_name_down" in s:
            h = DMesonJetUtils.GetObject(files[s["input_name"]], s["histogram_name_down"])
            if not h: exit(1)
            h_copy = h.Clone("{0}_down".format(s["name"]))
            h_copy.SetTitle(s["histogram_title_down"])
            variations.append(h_copy)
        if "histogram_name_up" in s:
            h = DMesonJetUtils.GetObject(files[s["input_name"]], s["histogram_name_up"])
            if not h: exit(1)
            h_copy = h.Clone("{0}_up".format(s["name"]))
            h_copy.SetTitle(s["histogram_title_up"])
            variations.append(h_copy)

    globalList.extend(variations)
    globalList.append(baseline)
    comp = DMesonJetCompare.DMesonJetCompare(config["name"])
    comp.fOptRatio = "hist"
    comp.fX1LegRatio = 0.45
    comp.fX1LegSpectrum = 0.45
    comp.fLegLineHeight = 0.05
    comp.fGridyRatio = True
    comp.fLogUpperSpace = 3
    r = comp.CompareSpectra(baseline, variations)
    for obj in r:
        globalList.append(obj)

def OpenFiles(config):
    files = dict()

    fname = "{0}/{1}/{1}.root".format(config["input_path"], config["default"]["input_name"])
    f = ROOT.TFile(fname)
    if not f or f.IsZombie():
        print("Could not open file '{0}'".format(fname))
        exit(1)
    files[config["default"]["input_name"]] = f

    for s in config["sources"]:
        if not s["input_name"] in files:
            fname = "{0}/{1}/{1}.root".format(config["input_path"], s["input_name"])
            f = ROOT.TFile(fname)
            if not f or f.IsZombie():
                print("Could not open file '{0}'".format(fname))
                exit(1)
            files[s["input_name"]] = f

    return files

def PlotFSspectraAndSyst(results):
    stat = GetSpectrum(results["default"], name)
    if not stat: return
    hname = name[name.rfind("/") + 1:]
    syst = results["SystematicUncertainty"][name]["{0}_CentralAsymmSyst".format(hname)]
    canvas = ROOT.TCanvas("{0}_canvas".format(name), "{0}_canvas".format(name))
    canvas.SetLogy()
    canvas.SetLeftMargin(0.13)
    canvas.cd()
    syst_copy = syst.Clone()
    syst_copy.SetFillColor(ROOT.kCyan + 1)
    syst_copy.SetMarkerColor(ROOT.kCyan + 1)
    syst_copy.GetYaxis().SetTitleOffset(1.5)
    syst_copy.Draw("a e2")
    stat_copy = stat.DrawCopy("same p e0 x0")
    stat_copy.Scale(1., "width")
    stat_copy.SetMarkerColor(ROOT.kBlue + 1)
    stat_copy.SetMarkerStyle(ROOT.kFullCircle)
    stat_copy.SetMarkerSize(0.6)
    stat_copy.SetLineColor(ROOT.kBlue + 1)
    leg = ROOT.TLegend(0.35, 0.71, 0.89, 0.89, "NB")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(20)
    leg.AddEntry(stat_copy, "Central values with stat. unc.", "pe")
    leg.AddEntry(syst_copy, "Systematic uncertainty", "f")
    leg.Draw()
    globalList.append(syst_copy)
    globalList.append(syst_copy.GetHistogram())
    globalList.append(stat_copy)
    globalList.append(leg)
    globalList.append(canvas)

def GenerateSystematicUncertainty(baseline, spectra):
    hname = baseline.GetName().replace("_copy", "")
    upperLimitsHist = ROOT.TH1D("{0}_UpperSyst".format(hname), "{0}_UpperSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    lowerLimitsHist = ROOT.TH1D("{0}_LowerSyst".format(hname), "{0}_LowerSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    symmetricLimitsHist = ROOT.TH1D("{0}_SymmSyst".format(hname), "{0}_SymmSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    relativeSystHist = ROOT.TH1D("{0}_RelSyst".format(hname), "{0}_RelSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    result = OrderedDict()
    result[upperLimitsHist.GetName()] = upperLimitsHist
    result[lowerLimitsHist.GetName()] = lowerLimitsHist
    result[symmetricLimitsHist.GetName()] = symmetricLimitsHist
    result[relativeSystHist.GetName()] = relativeSystHist
    for ibin in range(1, baseline.GetNbinsX() + 1):
        centralValue = baseline.GetBinContent(ibin)
        for var in spectra:
            diff = var.GetBinContent(ibin) - centralValue
            if diff > upperLimitsHist.GetBinContent(ibin):
                upperLimitsHist.SetBinContent(ibin, diff)
                print("Bin {0}, upper limit {1}".format(ibin, diff))
            if -diff > lowerLimitsHist.GetBinContent(ibin):
                lowerLimitsHist.SetBinContent(ibin, -diff)
                print("Bin {0}, lower limit {1}".format(ibin, -diff))
        symmetricLimitsHist.SetBinContent(ibin, max(upperLimitsHist.GetBinContent(ibin), lowerLimitsHist.GetBinContent(ibin)))
        relativeSystHist.SetBinContent(ibin, symmetricLimitsHist.GetBinContent(ibin) / baseline.GetBinContent(ibin))
    xArray = numpy.array([baseline.GetXaxis().GetBinCenter(ibin) for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    xArrayErr = numpy.array([baseline.GetXaxis().GetBinWidth(ibin) / 2 for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    yArray = numpy.array([baseline.GetBinContent(ibin) / baseline.GetXaxis().GetBinWidth(ibin) for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrUp = numpy.array([upperLimitsHist.GetBinContent(ibin) / upperLimitsHist.GetXaxis().GetBinWidth(ibin) for ibin in range(1, upperLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrLow = numpy.array([lowerLimitsHist.GetBinContent(ibin) / lowerLimitsHist.GetXaxis().GetBinWidth(ibin) for ibin in range(1, lowerLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrSym = numpy.array([symmetricLimitsHist.GetBinContent(ibin) / symmetricLimitsHist.GetXaxis().GetBinWidth(ibin) for ibin in range(1, symmetricLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    symmetricUncGraph = ROOT.TGraphErrors(baseline.GetNbinsX(), xArray, yArray, xArrayErr, yArrayErrSym)
    symmetricUncGraph.SetName("{0}_CentralSymmSyst".format(hname))
    symmetricUncGraph.GetXaxis().SetTitle(baseline.GetXaxis().GetTitle())
    symmetricUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [(mb) (GeV/#it{c})^{-1}]")
    asymmetricUncGraph = ROOT.TGraphAsymmErrors(baseline.GetNbinsX(), xArray, yArray, xArrayErr, xArrayErr, yArrayErrLow, yArrayErrUp)
    asymmetricUncGraph.SetName("{0}_CentralAsymmSyst".format(hname))
    asymmetricUncGraph.GetXaxis().SetTitle(baseline.GetXaxis().GetTitle())
    asymmetricUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [(mb) (GeV/#it{c})^{-1}]")
    result[symmetricUncGraph.GetName()] = symmetricUncGraph
    result[asymmetricUncGraph.GetName()] = asymmetricUncGraph
    return result

def SaveCanvases(input_path):
    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            oname = obj.GetName().replace("/", "_")
            obj.SaveAs("{0}/{1}.pdf".format(input_path, oname))

def GenerateRootList(pdict, name):
    rlist = ROOT.TList()
    rlist.SetName(name)
    print("Loop in {0}".format(name))
    for nobj, obj in pdict.iteritems():
        if isinstance(obj, ROOT.TObject):
            rlist.Add(obj)
        elif isinstance(obj, dict) or isinstance(obj, OrderedDict):
            print("Recursion in {0}".format(nobj))
            rlist.Add(GenerateRootList(obj, nobj))
        else:
            print("Error: type of object {0} not recognized!".format(obj))
    return rlist

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare B feed-down correction file.')
    parser.add_argument('yaml', metavar='conf.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config)

    IPython.embed()
