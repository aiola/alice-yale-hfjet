#!/usr/bin/env python
# python script to compare systematic unc

import argparse
import yaml
import IPython
import ROOT
import os
import DMesonJetCompare
import DMesonJetUtils
import RawYieldSpectrumLoader

globalList = []


def main(configs, name):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    print("{} configurations have been loaded".format(len(configs)))

    CompareSystematic_DMeson(configs, name)
    CompareSystematic_ByGroup(configs, name, "JetPtSpectrum", "DPtCutSideBand")
    CompareSystematic_ByGroup(configs, name, "JetPtSpectrum_DPt_20", "DPtBinWidth")
    CompareSystematic_ByGroup(configs, name, "JetPtSpectrum_DPt_20", "MethodDPt2")
    CompareSystematic_ByGroup(configs, name, "JetZSpectrum_DPt_20_JetPt_5_15", "DPtBinWidth")


def GetSystematicUncertainty(config, meson_name, jet_type, jet_radius, spectrum):
    loader = RawYieldSpectrumLoader.RawYieldSpectrumLoader(config["input_path"], config["train"], config["name"])
    loader.fDMeson = meson_name
    loader.fJetType = jet_type
    loader.fJetRadius = jet_radius
    loader.LoadSpectrumConfig(spectrum)
    h_syst = loader.GetRawYieldSystFromMultiTrial()
    h_yield = loader.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    h_syst.Divide(h_yield)
    h_syst.GetXaxis().SetTitle(h_yield.GetXaxis().GetTitle())
    h_syst.GetYaxis().SetTitle("rel. syst. unc.")
    return h_syst


def CompareSystematic_ByGroup(configs, name, spectrum_name, group_name):
    if not name: name = "ComparisonSystematic_{}_{}".format(spectrum_name, group_name)
    print(name)
    for c in configs:
        input_path = c["input_path"]
        print("Working on {}".format(c["name"]))
        if len(c["analysis"][0]["d_meson"]) < 2:
            print("Skipping {}, since there aren't enough different D meson cuts".format(c["name"]))
            continue
        for jet in c["analysis"][0]["jets"]:
            jet_type = jet["type"]
            jet_radius = jet["radius"]
            print("Working on jet {} {}".format(jet_type, jet_radius))
            for meson_name in c["analysis"][0]["d_meson"]:
                print("Working on D meson {}".format(meson_name))
                histos = []
                cname = "{}/{}/{}_{}".format(c["train"], c["name"], meson_name, name)
                for bin_list in c["analysis"][0]["binLists"]:
                    if not "spectra" in bin_list: continue
                    if not meson_name in bin_list["active_mesons"]: continue
                    for spectrum in bin_list["spectra"]:
                        if not meson_name in spectrum["active_mesons"]: continue
                        if not "multitrial" in spectrum: continue
                        if not meson_name in spectrum["multitrial"]: continue
                        if not spectrum_name in spectrum["name"]: continue
                        if not "compare" in spectrum: continue
                        if not group_name in spectrum["compare"]: continue
                        h = GetSystematicUncertainty(c, meson_name, jet_type, jet_radius, spectrum)
                        if not h: continue
                        h.SetTitle(spectrum["comp_titles"][spectrum["compare"].index(group_name)])
                        globalList.append(h)
                        histos.append(h)

                    if len(histos) < 2: continue

                    comp = DMesonJetCompare.DMesonJetCompare(cname)
                    comp.fOptSpectrumBaseline = "hist"
                    comp.fOptSpectrum = "hist"
                    comp.fX1LegRatio = 0.15
                    comp.fX1LegSpectrum = 0.25
                    comp.fLinUpperSpace = 0.4
                    comp.fGridyRatio = True
                    comp.fDoSpectraPlot = "lineary"
                    comp.fDoRatioPlot = False
                    r = comp.CompareSpectra(histos[0], histos[1:])
                    for obj in r:
                        if not obj in globalList:
                            globalList.append(obj)
                    comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasSpectra.GetName()))


def CompareSystematic_DMeson(configs, name):
    if not name: name = "ComparisonSystematic_DMesonCuts"
    print(name)
    for c in configs:
        input_path = c["input_path"]
        print("Working on {}".format(c["name"]))
        if len(c["analysis"][0]["d_meson"]) < 2:
            print("Skipping {}, since there aren't enough different D meson cuts".format(c["name"]))
            continue
        for jet in c["analysis"][0]["jets"]:
            jet_type = jet["type"]
            jet_radius = jet["radius"]
            print("Working on jet {} {}".format(jet_type, jet_radius))
            for bin_list in c["analysis"][0]["binLists"]:
                if not "spectra" in bin_list: continue
                for spectrum in bin_list["spectra"]:
                    if not "multitrial" in spectrum: continue
                    histos = []
                    cname = "{}/{}/{}_{}_{}".format(c["train"], c["name"], name, spectrum["name"], spectrum["suffix"])
                    for meson_name in c["analysis"][0]["d_meson"]:
                        if not meson_name in spectrum["active_mesons"]: continue
                        if not meson_name in bin_list["active_mesons"]: continue
                        if not meson_name in spectrum["multitrial"]: continue
                        print("Working on D meson {}".format(meson_name))
                        h = GetSystematicUncertainty(c, meson_name, jet_type, jet_radius, spectrum)
                        if not h: continue
                        h.SetTitle(meson_name)
                        globalList.append(h)
                        histos.append(h)

                    if len(histos) < 2: continue

                    comp = DMesonJetCompare.DMesonJetCompare(cname)
                    comp.fOptSpectrumBaseline = "hist"
                    comp.fOptSpectrum = "hist"
                    comp.fX1LegRatio = 0.15
                    comp.fX1LegSpectrum = 0.25
                    comp.fLinUpperSpace = 0.4
                    comp.fGridyRatio = True
                    comp.fDoSpectraPlot = "lineary"
                    comp.fDoRatioPlot = False
                    r = comp.CompareSpectra(histos[0], histos[1:])
                    for obj in r:
                        if not obj in globalList:
                            globalList.append(obj)
                    comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasSpectra.GetName()))


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
