#!/usr/bin/env python
# python script to compare efficiencies

import argparse
import yaml
import IPython
import ROOT
import os
import DMesonJetCompare
import DMesonJetUtils

globalList = []


def main(configs, name):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    print("{} configurations have been loaded".format(len(configs)))

    if len(configs) > 1:
        CompareEfficiency_Trains(configs, name, "JetPtDPtSpectrum", "_JetPt_500_3000")
        CompareEfficiency_Trains(configs, name, "JetPtSpectrum_DPt_30", "")
        CompareEfficiency_Trains(configs, name, "JetZSpectrum_DPt_20_JetPt_5_15", "")
        CompareEfficiency_Trains(configs, name, "JetZSpectrum_DPt_60_JetPt_15_30", "")

    CompareEfficiency_DMeson(configs, name)
    CompareEfficiency_Type(configs, name)

    definitions = [{"spectrum" : "JetPtDPtSpectrum", "suffix" : "_JetPt_500_3000", "title" : "Only efficiency"}, {"spectrum" : "DPtSpectrum_JetPt_5_30", "suffix" : "", "title" : "Eff. w/ corr. 5 < #it{p}_{T,ch jet} < 30 GeV/#it{c}"}]
    CompareEfficiency_Definitions(configs, "Comparison_KineCuts_JetPt_5_30", definitions, 3, 29.9)

    definitions = [{"spectrum" : "JetPtDPtSpectrum", "suffix" : "_JetPt_500_1500", "title" : "Only efficiency"}, {"spectrum" : "DPtSpectrum_JetPt_5_15", "suffix" : "", "title" : "Eff. w/ corr. 5 < #it{p}_{T,ch jet} < 15 GeV/#it{c}"}]
    CompareEfficiency_Definitions(configs, "Comparison_KineCuts_JetPt_5_15", definitions, 2, 14.9)

    definitions = [{"spectrum" : "JetPtDPtSpectrum", "suffix" : "_JetPt_1500_3000", "title" : "Only efficiency"}, {"spectrum" : "DPtSpectrum_JetPt_15_30", "suffix" : "", "title" : "Eff. w/ corr. 15 < #it{p}_{T,ch jet} < 30 GeV/#it{c}"}]
    CompareEfficiency_Definitions(configs, "Comparison_KineCuts_JetPt_15_30", definitions, 6, 29.9)

    definitions = [{"spectrum" : "JetPtDPtSpectrum", "suffix" : "_JetPt_500_1500", "title" : "5 < #it{p}_{T,ch jet} < 15 GeV/#it{c}"}, {"spectrum" : "JetPtDPtSpectrum", "suffix" : "_JetPt_1500_3000", "title" : "15 < #it{p}_{T,ch jet} < 30 GeV/#it{c}"}]
    CompareEfficiency_Definitions(configs, "Comparison_JetPtRanges", definitions, 2, 14.9, "")


def GetEfficiency(config, meson_name, jet_type, jet_radius, spectrum="JetPtDPtSpectrum", suffix="_JetPt_500_3000"):
    return GetHisto("Efficiency", config, meson_name, jet_type, jet_radius, spectrum, suffix)


def GetHisto(histo_name, config, meson_name, jet_type, jet_radius, spectrum="JetPtDPtSpectrum", suffix="_JetPt_500_3000"):
    fileName = "{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"])
    file = ROOT.TFile(fileName)
    if not file or file.IsZombie():
        print("Could not open file '{}'".format(fileName))
        exit(1)

    hname = "{meson_name}_Jet_AKT{jet_type}{jet_radius}_pt_scheme_{spectrum}/{meson_name}_Jet_AKT{jet_type}{jet_radius}_pt_scheme_{spectrum}_{histo_name}{suffix}".format(
        meson_name=meson_name, jet_type=jet_type, jet_radius=jet_radius, spectrum=spectrum, histo_name=histo_name, suffix=suffix)
    h = DMesonJetUtils.GetObject(file, hname)
    return h


def CompareEfficiency_Definitions(configs, name, definitions, minpt, maxpt, optratio="hist"):
    print("CompareEfficiency_Definitions")
    if not name: name = "Comparison_KineCuts"
    efficiency_type = "Prompt"
    print("Working on {}".format(efficiency_type))
    for c in configs:
        input_path = c["input_path"]
        print("Working on {}".format(c["name"]))
        for jet in c["analysis"][0]["jets"]:
            jet_type = jet["type"]
            jet_radius = jet["radius"]
            print("Working on jet {} {}".format(jet_type, jet_radius))
            for meson_cuts in c["analysis"][0]["d_meson_cuts"]:
                meson_name = "{}_{}_{}".format(efficiency_type, c["analysis"][0]["d_meson"][0], meson_cuts)
                print("Working on D meson {}".format(meson_name))
                histos = []
                cname = "{}/{}/{}_{}".format(c["train"], c["name"], meson_cuts, name)
                for definition in definitions:
                    print("Definition '{}'".format(definition["title"]))
                    h = GetEfficiency(c, meson_name, jet_type, jet_radius, definition["spectrum"], definition["suffix"])
                    if not h: continue
                    h.SetTitle(definition["title"])
                    globalList.append(h)
                    histos.append(h)

                if len(histos) < 2: continue
                comp = DMesonJetCompare.DMesonJetCompare(cname)
                comp.fOptRatio = optratio
                comp.fX1LegRatio = 0.15
                comp.fX1LegSpectrum = 0.25
                comp.fLinUpperSpace = 0.4
                comp.fGridyRatio = True
                comp.fDoSpectraPlot = "lineary"
                r = comp.CompareSpectra(histos[0], histos[1:])
                for obj in r:
                    if not obj in globalList:
                        globalList.append(obj)

                comp.fMainHistogram.GetXaxis().SetRangeUser(minpt, maxpt)
                comp.fMainRatioHistogram.GetXaxis().SetRangeUser(minpt, maxpt)
                comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasSpectra.GetName()))
                comp.fCanvasRatio.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasRatio.GetName()))


def CompareEfficiency_Type(configs, name):
    print("CompareEfficiency_DMeson")
    if not name: name = "Comparison"
    efficiency_types = ["Prompt", "NonPrompt"]

    for c in configs:
        input_path = c["input_path"]
        print("Working on {}".format(c["name"]))
        if len(c["analysis"][0]["d_meson_cuts"]) < 2:
            print("Skipping {}, since there aren't enough different D meson cuts".format(c["name"]))
            continue
        for jet in c["analysis"][0]["jets"]:
            jet_type = jet["type"]
            jet_radius = jet["radius"]
            print("Working on jet {} {}".format(jet_type, jet_radius))
            cname = "{}/{}/{}_Prompt_NonPrompt".format(c["train"], c["name"], name)
            comp = DMesonJetCompare.DMesonJetCompare(cname)
            comp.fOptRatio = "hist"
            comp.fX1LegRatio = 0.15
            comp.fX1LegSpectrum = 0.25
            comp.fLinUpperSpace = 0.45
            comp.fGridyRatio = True
            comp.fDoSpectraPlot = "lineary"
            colors = [[ROOT.kOrange + 2, ROOT.kRed + 2], [ROOT.kGreen + 2, ROOT.kBlue + 2], [ROOT.kAzure + 2, ROOT.kCyan + 2]]
            markers = [[ROOT.kOpenCircle, ROOT.kFullCircle], [ROOT.kOpenSquare, ROOT.kFullSquare], [ROOT.kOpenDiamond, ROOT.kFullDiamond]]
            lines = [[2, 1], [9, 5], [7, 10]]

            for meson_cuts, cols, marks, lins in zip(c["analysis"][0]["d_meson_cuts"], colors, markers, lines):
                histos = []
                for efficiency_type in efficiency_types:
                    print("Working on {}".format(efficiency_type))
                    meson_name = "{}_{}_{}".format(efficiency_type, c["analysis"][0]["d_meson"][0], meson_cuts)
                    print("Working on D meson {}".format(meson_name))
                    h = GetEfficiency(c, meson_name, jet_type, jet_radius)
                    if not h: continue
                    h.SetBinContent(1, 0)
                    h.SetBinContent(2, 0)
                    h.SetTitle("{}, {}".format(meson_cuts, efficiency_type))
                    globalList.append(h)
                    histos.append(h)

                if len(histos) < 2: continue
                comp.fColors = cols
                comp.fMarkers = marks
                comp.fLines = lins
                r = comp.CompareSpectra(histos[0], histos[1:])
                for obj in r:
                    if not obj in globalList:
                        globalList.append(obj)
                if not "same" in comp.fOptSpectrumBaseline:
                    comp.fOptSpectrumBaseline += " same"

            comp.fMainHistogram.GetXaxis().SetRangeUser(2, 29.9)
            comp.fMainRatioHistogram.GetXaxis().SetRangeUser(2, 29.9)
            comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasSpectra.GetName()))
            comp.fCanvasRatio.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasRatio.GetName()))


def CompareEfficiency_DMeson(configs, name):
    print("CompareEfficiency_DMeson")
    if not name: name = "Comparison"
    efficiency_types = ["Prompt", "NonPrompt"]
    for efficiency_type in efficiency_types:
        print("Working on {}".format(efficiency_type))
        for c in configs:
            input_path = c["input_path"]
            print("Working on {}".format(c["name"]))
            if len(c["analysis"][0]["d_meson_cuts"]) < 2:
                print("Skipping {}, since there aren't enough different D meson cuts".format(c["name"]))
                continue
            for jet in c["analysis"][0]["jets"]:
                jet_type = jet["type"]
                jet_radius = jet["radius"]
                print("Working on jet {} {}".format(jet_type, jet_radius))
                histos = []
                cname = ""
                for meson_cuts in c["analysis"][0]["d_meson_cuts"]:
                    meson_name = "{}_{}_{}".format(efficiency_type, c["analysis"][0]["d_meson"][0], meson_cuts)
                    print("Working on D meson {}".format(meson_name))
                    h = GetEfficiency(c, meson_name, jet_type, jet_radius)
                    if not h: continue
                    h.SetTitle(meson_cuts)
                    if cname: cname = "{}_{}".format(cname, meson_cuts)
                    else: cname = meson_cuts
                    globalList.append(h)
                    histos.append(h)

                if len(histos) < 2: continue
                cname = "{}/{}/{}_{}_{}".format(c["train"], c["name"], name, efficiency_type, cname)
                comp = DMesonJetCompare.DMesonJetCompare(cname)
                comp.fOptRatio = "hist"
                comp.fX1LegRatio = 0.15
                comp.fX1LegSpectrum = 0.25
                comp.fLinUpperSpace = 0.4
                comp.fGridyRatio = True
                comp.fDoSpectraPlot = "lineary"
                r = comp.CompareSpectra(histos[0], histos[1:])
                for obj in r:
                    if not obj in globalList:
                        globalList.append(obj)

                comp.fMainHistogram.GetXaxis().SetRangeUser(2, 29.9)
                comp.fMainRatioHistogram.GetXaxis().SetRangeUser(2, 29.9)
                comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasSpectra.GetName()))
                comp.fCanvasRatio.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasRatio.GetName()))


def CompareEfficiency_Trains(configs, name, spectrum, suffix):
    print("CompareEfficiency_Trains")
    if not name: name = "Comparison"
    input_path = configs[0]["input_path"]

    efficiency_types = ["Prompt", "NonPrompt"]
    histo_names = ["Efficiency", "Truth", "Measured"]

    for histo_name in histo_names:
        print("Working on {}".format(histo_name))
        for efficiency_type in efficiency_types:
            print("Working on {}".format(efficiency_type))
            for meson_cuts in configs[0]["analysis"][0]["d_meson_cuts"]:
                meson_name = "{}_{}_{}".format(efficiency_type, configs[0]["analysis"][0]["d_meson"][0], meson_cuts)
                print("Working on {}".format(meson_name))
                for jet in configs[0]["analysis"][0]["jets"]:
                    jet_type = jet["type"]
                    jet_radius = jet["radius"]
                    print("Working on {}{}".format(jet_type, jet_radius))
                    histos = []
                    cname = ""
                    for c in configs:
                        print("Working on {}".format(c["name"]))
                        h = GetHisto(histo_name, c, meson_name, jet_type, jet_radius, spectrum, suffix)
                        if not h:
                            print("Histogram not found, skipping")
                            continue
                        h.SetTitle(c["name"])
                        if cname: cname = "{}_{}".format(cname, c["name"])
                        else: cname = c["name"]
                        globalList.append(h)
                        histos.append(h)

                    if len(histos) < 2:
                        print("Too few histograms, skipping")
                        continue
                    cname = "{}_{}_{}{}_{}_{}".format(name, histo_name, spectrum, suffix, efficiency_type, cname)
                    comp = DMesonJetCompare.DMesonJetCompare(cname)
                    if histo_name == "Efficiency" or "JetZSpectrum" in spectrum:
                        comp.fDoSpectraPlot = "lineary"
                    else:
                        comp.fDoSpectraPlot = "logy"
                    comp.fOptRatio = "hist"
                    comp.fX1LegRatio = 0.15
                    comp.fX1LegSpectrum = 0.25
                    comp.fLogUpperSpace = 5  # this factor will be used to adjust the y axis in log scale
                    comp.fLogLowerSpace = 1.5  # this factor will be used to adjust the y axis in log scale
                    comp.fLinUpperSpace = 0.45  # this factor will be used to adjust the y axis in linear scale
                    comp.fLinLowerSpace = 0.15  # this factor will be used to adjust the y axis in linear scale
                    comp.fGridyRatio = True
                    r = comp.CompareSpectra(histos[0], histos[1:])
                    for obj in r:
                        if not obj in globalList:
                            globalList.append(obj)

                    comp.fMainHistogram.GetXaxis().SetRangeUser(2, 29.9)
                    comp.fMainRatioHistogram.GetXaxis().SetRangeUser(2, 29.9)
                    comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasSpectra.GetName()))
                    comp.fCanvasRatio.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasRatio.GetName()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compares efficiencies.')
    parser.add_argument('yaml', nargs='*',
                        help='List of YAML configuration files')
    parser.add_argument('--name', metavar='NAME',
                        default=None)
    args = parser.parse_args()

    configs = []

    for fname in args.yaml:
        f = open(fname, 'r')
        c = yaml.load(f)
        f.close()
        configs.append(c)

    main(configs, args.name)

    IPython.embed()
