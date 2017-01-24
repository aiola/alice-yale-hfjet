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
    ROOT.gStyle.SetOptStat(0)

    results = Start(config)

    PlotSystematicUncertaintySummary(config["name"], results)

    output_file_name = "{0}/{1}.root".format(config["input_path"], config["name"])
    print("Storing results in {0}".format(output_file_name))
    outputFile = ROOT.TFile(output_file_name, "recreate")
    if not outputFile or outputFile.IsZombie():
       print("Could not open output file {0}".format(output_file_name))
       exit(1)
    for name, r in results.iteritems():
       rlist = GenerateRootList(r, name)
       rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)
    outputFile.Close()
    globalList.append(rlist)

    SaveCanvases(config["input_path"])

def Start(config):
    histograms = LoadHistograms(config)
    results = dict()
    results["Variations"] = CompareVariations(config, histograms)
    results["Uncertainties"] = GenerateUncertainties(config, histograms)
    results["FinalSpectrum"] = PlotSpectrumStatAndSyst(config["name"], results)
    return results

def LoadHistograms(config):
    crossSection = 62.2  # mb CINT1
    branchingRatio = 0.0393  # D0->Kpi
    antiPartNorm = 2.0  # D0 / D0bar

    files = OpenFiles(config)
    events = LoadEvents(files, config)

    histograms = dict()

    h = DMesonJetUtils.GetObject(files[config["default"]["input_name"]], config["default"]["histogram_name"])
    h.Scale(crossSection / (events[config["default"]["input_name"]] * branchingRatio * antiPartNorm), "width")
    if not h: exit(1)
    histograms["default"] = h
    for s in config["sources"]:
        if not s["active"]: continue
        histograms[s["name"]] = dict()
        if "input_name" in s:
            input_name_up = s["input_name"]
            input_name_down = s["input_name"]
        if "input_name_up" in s: input_name_up = s["input_name_up"]
        if "input_name_down" in s: input_name_down = s["input_name_down"]
        if "histogram_name_down" in s:
            h = DMesonJetUtils.GetObject(files[input_name_down], s["histogram_name_down"])
            h.Scale(crossSection / (events[input_name_down] * branchingRatio * antiPartNorm), "width")
            if not h: exit(1)
            histograms[s["name"]]["down"] = h
        if "histogram_name_up" in s:
            h = DMesonJetUtils.GetObject(files[input_name_up], s["histogram_name_up"])
            h.Scale(crossSection / (events[input_name_up] * branchingRatio * antiPartNorm), "width")
            if not h: exit(1)
            histograms[s["name"]]["up"] = h

    return histograms

def LoadEvents(files, config):
    events = dict()
    for name, f in files.iteritems():
        slash = config["default"]["histogram_name"].find("/")
        heventsName = config["default"]["histogram_name"][:slash + 1]
        heventsName += "Events"
        hevents = DMesonJetUtils.GetObject(files[name], heventsName)
        events[name] = hevents.GetBinContent(1)
    return events

def OpenFiles(config):
    files = dict()

    fname = "{0}/{1}/{1}.root".format(config["input_path"], config["default"]["input_name"])
    f = ROOT.TFile(fname)
    if not f or f.IsZombie():
        print("Could not open file '{0}'".format(fname))
        exit(1)
    files[config["default"]["input_name"]] = f

    for s in config["sources"]:
        if not s["active"]: continue
        input_names = LoadInputNames(s)
        for input_name in input_names:
            if not input_name in files:
                fname = "{0}/{1}/{1}.root".format(config["input_path"], input_name)
                f = ROOT.TFile(fname)
                if not f or f.IsZombie():
                    print("Could not open file '{0}'".format(fname))
                    exit(1)
                files[input_name] = f
    return files

def LoadInputNames(s):
    input_names = []
    if "input_name" in s: input_names.append(s["input_name"])
    if "input_name_up" in s: input_names.append(s["input_name_up"])
    if "input_name_down" in s: input_names.append(s["input_name_down"])
    return input_names

def CompareVariations(config, histograms):
    result = dict()
    variations = []
    h = histograms["default"]
    baseline = h.Clone("default")
    baseline.SetTitle(config["default"]["title"])
    result[baseline.GetName()] = baseline
    for s in config["sources"]:
        if not s["active"]: continue
        if "down" in histograms[s["name"]]:
            h = histograms[s["name"]]["down"]
            h_copy = h.Clone("{0}_down".format(s["name"]))
            h_copy.SetTitle(s["histogram_title_down"])
            result[h_copy.GetName()] = h_copy
            variations.append(h_copy)
        if "up" in histograms[s["name"]]:
            h = histograms[s["name"]]["up"]
            h_copy = h.Clone("{0}_up".format(s["name"]))
            h_copy.SetTitle(s["histogram_title_up"])
            result[h_copy.GetName()] = h_copy
            variations.append(h_copy)

    globalList.extend(variations)
    globalList.append(baseline)
    comp = DMesonJetCompare.DMesonJetCompare("CompareVariations_{0}".format(config["name"]))
    comp.fOptRatio = "hist"
    comp.fX1LegRatio = 0.45
    comp.fX1LegSpectrum = 0.45
    comp.fLegLineHeight = 0.05
    comp.fGridyRatio = True
    comp.fLogUpperSpace = 3
    r = comp.CompareSpectra(baseline, variations)
    for obj in r:
        globalList.append(obj)
    return result

def CalculateFixSystematicUncertainty(config):
    fixed_unc2 = 0
    print("Source & Uncertainty (\\%) \\\\ \\hline")
    for u in config["fixed_uncertainties"]:
        fixed_unc2 += u["uncertainty"] ** 2
        print("{0} & {1:.1f} \\\\".format(u["title"], u["uncertainty"] * 100))
    fixed_unc = math.sqrt(fixed_unc2)
    print("\\hline")
    print("{0} & {1:.1f}".format("Total", fixed_unc * 100))
    return fixed_unc

def GenerateUncertainties(config, histograms):
    baseline = histograms["default"]
    fixed_syst_unc = CalculateFixSystematicUncertainty(config)
    partialRelSystUnc = [DMesonJetUtils.soft_clone(baseline, s["name"], s["title"], "relative uncertainty") for s in config["sources"] if s["active"]]
    totRelSystUnc = DMesonJetUtils.soft_clone(baseline, "tot_rel_syst_unc", "Total Systematic Uncertainty", "relative uncertainty")
    statUnc = DMesonJetUtils.soft_clone(baseline, "stat_unc", "Statistical Uncertainty", "relative uncertainty")
    totUnc = DMesonJetUtils.soft_clone(baseline, "tot_unc", "Total Uncertainty", "relative uncertainty")
    centralSystUnc = DMesonJetUtils.soft_clone(baseline, "central_syst_unc")

    for ibin in range(1, baseline.GetNbinsX() + 1):
        tot_syst_unc2 = 0
        ivar = 0
        for s in config["sources"]:
            if not s["active"]: continue
            if "down" in histograms[s["name"]]:
                h = histograms[s["name"]]["down"]
                diff_down = abs(baseline.GetBinContent(ibin) - h.GetBinContent(ibin))
            else:
                diff_down = 0
            if "up" in histograms[s["name"]]:
                h = histograms[s["name"]]["up"]
                diff_up = abs(baseline.GetBinContent(ibin) - h.GetBinContent(ibin))
            else:
                diff_up = 0
            part_unc = max(diff_down, diff_up)
            partialRelSystUnc[ivar].SetBinContent(ibin, part_unc / baseline.GetBinContent(ibin))
            tot_syst_unc2 += part_unc ** 2
            ivar += 1
        tot_syst_unc2 += (fixed_syst_unc * baseline.GetBinContent(ibin)) ** 2
        tot_syst_unc = math.sqrt(tot_syst_unc2)
        stat_unc = baseline.GetBinError(ibin)
        tot_unc = math.sqrt(tot_syst_unc2 + stat_unc ** 2)
        totRelSystUnc.SetBinContent(ibin, tot_syst_unc / baseline.GetBinContent(ibin))
        statUnc.SetBinContent(ibin, stat_unc / baseline.GetBinContent(ibin))
        totUnc.SetBinContent(ibin, tot_unc / baseline.GetBinContent(ibin))
        centralSystUnc.SetBinContent(ibin, baseline.GetBinContent(ibin))
        centralSystUnc.SetBinError(ibin, tot_syst_unc)

    result = {"PartialSystematicUncertainties" : partialRelSystUnc, "correlated_uncertainty" : fixed_syst_unc, totRelSystUnc.GetName() : totRelSystUnc, centralSystUnc.GetName() : centralSystUnc, statUnc.GetName() : statUnc, totUnc.GetName() : totUnc}
    return result

def PlotSystematicUncertaintySummary(name, results):
    sources = []
    h = results["Uncertainties"]["tot_unc"]
    baseline = h.Clone("{0}_copy".format(h.GetName()))

    print("Source & \\multicolumn{{{0}}}{{c}}{{Uncertainty (\\%)}} \\\\ \\hline".format(baseline.GetNbinsX()))
    print(" & ".join(["\\ptchjet (\\GeVc)"] + ["{0:.0f} - {1:.0f}".format(baseline.GetXaxis().GetBinLowEdge(ibin), baseline.GetXaxis().GetBinUpEdge(ibin)) for ibin in range(1, baseline.GetNbinsX() + 1)]) + "\\\\ \hline")

    h = results["Uncertainties"]["stat_unc"]
    stat_unc = h.Clone("{0}_copy".format(h.GetName()))
    sources.append(stat_unc)

    h = results["Uncertainties"]["tot_rel_syst_unc"]
    tot_rel_syst_unc = h.Clone("{0}_copy".format(h.GetName()))
    sources.append(tot_rel_syst_unc)
    for h in results["Uncertainties"]["PartialSystematicUncertainties"]:
        print(" & ".join([h.GetTitle()] + ["{0:.1f}".format(h.GetBinContent(ibin) * 100) for ibin in range(1, h.GetNbinsX() + 1)]) + "\\\\")
        h_copy = h.Clone("{0}_copy".format(h.GetName()))
        sources.append(h_copy)

    print("\\hline")
    print("Correlated Uncertainty & \\multicolumn{{{0}}}{{c}}{{{1:.1f}}} \\\\".format(baseline.GetNbinsX(), results["Uncertainties"]["correlated_uncertainty"] * 100))
    print("\\hline")
    print(" & ".join([tot_rel_syst_unc.GetTitle()] + ["{0:.1f}".format(tot_rel_syst_unc.GetBinContent(ibin) * 100) for ibin in range(1, tot_rel_syst_unc.GetNbinsX() + 1)]) + "\\\\")
    print("\\hline")
    print(" & ".join([stat_unc.GetTitle()] + ["{0:.1f}".format(stat_unc.GetBinContent(ibin) * 100) for ibin in range(1, stat_unc.GetNbinsX() + 1)]) + "\\\\")
    print("\\hline")
    print(" & ".join([baseline.GetTitle()] + ["{0:.1f}".format(baseline.GetBinContent(ibin) * 100) for ibin in range(1, baseline.GetNbinsX() + 1)]) + "\\\\")

    globalList.extend(sources)
    globalList.append(baseline)
    comp = DMesonJetCompare.DMesonJetCompare("CompareUncertainties_{0}".format(config["name"]))
    comp.fOptSpectrum = "hist"
    comp.fOptSpectrumBaseline = "hist"
    comp.fX1LegSpectrum = 0.45
    comp.fLegLineHeight = 0.05
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = False
    r = comp.CompareSpectra(baseline, sources)
    for obj in r:
        globalList.append(obj)

def PlotSpectrumStatAndSyst(name, results):
    stat = results["Variations"]["default"]
    syst = results["Uncertainties"]["central_syst_unc"]
    canvas = ROOT.TCanvas(name, name)
    canvas.SetLogy()
    canvas.SetLeftMargin(0.13)
    canvas.cd()
    syst_copy = syst.DrawCopy("e2")
    syst_copy.SetFillColor(ROOT.kGray)
    syst_copy.SetMarkerColor(ROOT.kGray)
    syst_copy.SetLineColor(ROOT.kGray)
    syst_copy.GetYaxis().SetTitleOffset(1.5)
    syst_copy.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]")
    syst_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
    stat_copy = stat.DrawCopy("same p e0 x0")
    stat_copy.SetMarkerColor(ROOT.kRed + 2)
    stat_copy.SetLineColor(ROOT.kRed + 2)
    stat_copy.SetMarkerStyle(ROOT.kFullCircle)
    stat_copy.SetMarkerSize(0.9)
    leg = ROOT.TLegend(0.35, 0.71, 0.89, 0.89, "NB")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(20)
    leg.AddEntry(stat_copy, "ALICE", "pe")
    leg.AddEntry(syst_copy, "Systematic Uncertainty", "f")
    leg.Draw()
    globalList.append(syst_copy)
    globalList.append(stat_copy)
    globalList.append(leg)
    globalList.append(canvas)
    result = dict()
    syst_copy_copy = syst_copy.Clone("CentralPointsSystematicUncertainty")
    stat_copy_copy = stat_copy.Clone("CentralPointsStatisticalUncertainty")
    result["CentralPointsSystematicUncertainty"] = syst_copy_copy
    result["CentralPointsStatisticalUncertainty"] = stat_copy_copy
    result["FinalSpectrumCanvas"] = canvas
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
    if isinstance(pdict, dict) or isinstance(pdict, OrderedDict):
        objects = pdict.iteritems()
    else:
        objects = zip([inner_obj.GetName() for inner_obj in pdict], pdict)
    for nobj, obj in objects:
        if isinstance(obj, ROOT.TObject):
            print("Adding in {0}".format(nobj))
            rlist.Add(obj)
        elif isinstance(obj, dict) or isinstance(obj, OrderedDict) or isinstance(obj, list):
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
