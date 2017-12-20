#!/usr/bin/env python
# python script to compare efficiencies

import argparse
import yaml
import IPython
import ROOT
import os
import DMesonJetUtils
import DMesonJetCompare
import RawYieldSpectrumLoader
import numpy

globalList = []


def main(configEff, configData):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    EfficiencyUncertainty(configEff, "JetPtDPtSpectrum", "_JetPt_500_3000", configData, "DPt", "JetPt_5_30", "InvMassFit", 2, 29.9)


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


def EfficiencyUncertainty(configEff, spectrum, suffix, configData, var, kincuts, method, minpt, maxpt):
    wrap = RawYieldSpectrumLoader.RawYieldSpectrumLoader(configData["input_path"], configData["train"], configData["name"])
    wrap.fVariableName = var
    wrap.fKinematicCuts = kincuts
    wrap.fRawYieldMethod = method
    wrap.fUseReflections = False
    efficiency_type = "Prompt"
    name = "EfficiencyUncertainty"
    print("Working on {}".format(efficiency_type))
    input_path = configEff["input_path"]
    print("Working on {}".format(configEff["name"]))
    for jet in configEff["analysis"][0]["jets"]:
        jet_type = jet["type"]
        jet_radius = jet["radius"]
        wrap.fJetType = jet_type
        wrap.fJetRadius = jet_radius
        print("Working on jet {} {}".format(jet_type, jet_radius))
        for meson_cuts in configEff["analysis"][0]["d_meson_cuts"]:
            meson_name = "{}_{}_{}".format(efficiency_type, configEff["analysis"][0]["d_meson"][0], meson_cuts)
            wrap.fDMeson = "{}_{}".format(configEff["analysis"][0]["d_meson"][0], meson_cuts)
            print("Working on D meson {}".format(meson_name))
            hEff = GetEfficiency(configEff, meson_name, jet_type, jet_radius, spectrum, suffix)
            if not hEff: continue

            hEffStat = DMesonJetUtils.GetRelativeUncertaintyHistogram(hEff)
            hEffStat.GetYaxis().SetTitle("relative statistical uncertainty")
            hEffStat.SetTitle("Efficiency")
            globalList.append(hEffStat)

            hData = wrap.GetDefaultSpectrumFromDMesonJetAnalysis(0, 0, 0)
            hDataStat = DMesonJetUtils.GetRelativeUncertaintyHistogram(hData)
            hDataStat.GetYaxis().SetTitle("relative statistical uncertainty")
            hDataStat.SetTitle("Data")

            # this is a bad hack to fix a different binning between data and efficiency
            # need to go from [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30]
            # to [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 30]
            dptbins = [2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 30]
            hEffStat_copy = ROOT.TH1F("hEffStat_copy", "Efficiency", len(dptbins) - 1, numpy.array(dptbins, dtype=numpy.float64))
            hEffStat_copy.GetYaxis().SetTitle("relative statistical uncertainty")
            skip = 2
            for ibin in range(1, hEffStat_copy.GetNbinsX() + 1):
                if ibin == 7 or ibin == 10:
                    effValue = (hEffStat.GetBinContent(ibin + skip) + hEffStat.GetBinContent(ibin + skip + 1)) / 2
                    skip += 1
                else:
                    effValue = hEffStat.GetBinContent(ibin + skip)
                hEffStat_copy.SetBinContent(ibin, effValue)
            globalList.append(hEffStat_copy)

            globalList.append(hDataStat)

            cname = "{}/{}/{}_{}{}_{}".format(configEff["train"], configEff["name"], meson_cuts, spectrum, suffix, name)

            comp = DMesonJetCompare.DMesonJetCompare(cname)
            comp.fDoRatioPlot = "lineary"
            comp.fDoSpectraPlot = "lineary"
            comp.fLinUpperSpace = 0.05
            comp.fLinLowerSpace = 0.05
            comp.fOptSpectrumBaseline = "hist"
            comp.fOptSpectrum = "hist"
            comp.fGridySpectrum = True
            r = comp.CompareSpectra(hEffStat_copy, [hDataStat])
            for obj in r:
                if not obj in globalList:
                    globalList.append(obj)

            comp.fMainHistogram.GetXaxis().SetRangeUser(minpt, maxpt)
            comp.fMainRatioHistogram.GetXaxis().SetRangeUser(minpt, maxpt)
            comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasSpectra.GetName()))
            comp.fCanvasRatio.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasRatio.GetName()))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compares efficiencies.')
    parser.add_argument('--eff',
                        help='YAML configuration file for efficiency')
    parser.add_argument('--data',
                        help='YAML configuration file for data')
    args = parser.parse_args()

    configs = []

    f = open(args.eff, 'r')
    cEff = yaml.load(f)
    f.close()

    f = open(args.data, 'r')
    cData = yaml.load(f)
    f.close()

    main(cEff, cData)

    IPython.embed()
