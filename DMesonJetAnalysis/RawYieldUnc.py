#!/usr/bin/env python
# python script to do extract B feed down correction factors

import argparse
import yaml
import IPython
import ROOT
import os
import DMesonJetCompare
import RawYieldSpectrumLoader

globalList = []

wrap = RawYieldSpectrumLoader.RawYieldSpectrumLoader()

xaxisTitle = ""
yaxisTitle = ""

def main(config, meson_name, jet_type, jet_radius, spectrum):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    wrap.fInputPath = config["input_path"]
    wrap.fTrain = config["train"]
    wrap.fAnalysisName = config["name"]

    comp = DMesonJetCompare.DMesonJetCompare("AverageRawYieldVsDefault")
    comp.fOptSpectrum = "p x0"
    comp.fOptRatio = ""
    comp.fSeparateBaselineUncertainty = True
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.25
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.25  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = 0.05  # this factor will be used to adjust the y axis in linear scale
    comp.fColors = [ROOT.kBlue + 2] * 2
    comp.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle]
    comp.fFills = [3245] * 2
    comp.fOptSpectrumBaseline = "e2"
    default_vs_average_raw_yield(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, spectrum)
    comp.fColors = [ROOT.kRed + 2] * 2
    comp.fMarkers = [ROOT.kOpenSquare, ROOT.kFullSquare]
    comp.fFills = [3254] * 2
    comp.fOptSpectrumBaseline = "e2 same"
    default_vs_average_raw_yield(comp, "SideBand", config, meson_name, jet_type, jet_radius, spectrum)

    comp = DMesonJetCompare.DMesonJetCompare("DefaultMTRawYieldVsDefault")
    comp.fOptRatio = "hist"
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.20
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.3  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = 0.1  # this factor will be used to adjust the y axis in linear scale
    comp.fNoErrorInBaseline = True
    comp.fColors = [ROOT.kGreen + 2, ROOT.kBlue + 2]
    comp.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle]
    default_vs_default_mt(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, spectrum)
    comp.fColors = [ROOT.kOrange + 2, ROOT.kRed + 2]
    comp.fMarkers = [ROOT.kOpenSquare, ROOT.kFullSquare]
    comp.fOptSpectrumBaseline = "same"
    default_vs_default_mt(comp, "SideBand", config, meson_name, jet_type, jet_radius, spectrum)

    comp = DMesonJetCompare.DMesonJetCompare("DefaultMTRawYieldVsDefaultUncertainties")
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = False
    comp.fOptSpectrum = "hist"
    comp.fOptSpectrumBaseline = "hist"
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.20
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.3  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = 0.1  # this factor will be used to adjust the y axis in linear scale
    comp.fNoErrorInBaseline = True
    comp.fColors = [ROOT.kBlue + 2, ROOT.kBlue + 2]
    default_vs_default_mt_unc(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, spectrum)
    comp.fColors = [ROOT.kRed + 2, ROOT.kRed + 2]
    comp.fOptSpectrumBaseline = "same hist"
    default_vs_default_mt_unc(comp, "SideBand", config, meson_name, jet_type, jet_radius, spectrum)

    comp = DMesonJetCompare.DMesonJetCompare("ReflectionComparison")
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.20
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.15  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = -0.4  # this factor will be used to adjust the y axis in linear scale
    comp.fGridyRatio = True
    comp.fNoErrorInBaseline = True
    comp.fColors = [ROOT.kGreen + 2, ROOT.kBlue + 2]
    comp.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle]
    reflections_raw_yield(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, spectrum)
    comp.fColors = [ROOT.kOrange + 2, ROOT.kRed + 2]
    comp.fMarkers = [ROOT.kOpenSquare, ROOT.kFullSquare]
    comp.fOptSpectrumBaseline = "same"
    reflections_raw_yield(comp, "SideBand", config, meson_name, jet_type, jet_radius, spectrum)

    outputPath = "{0}/{1}/{2}/RawYieldUnc/pdf".format(config["input_path"], config["train"], config["name"])
    if not os.listdir(outputPath): os.makedirs(outputPath)
    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            fname = "{0}/{1}.pdf".format(outputPath, obj.GetName())
            obj.SaveAs(fname)

def default_vs_default_mt(comp, method, config, meson_name, jet_type, jet_radius, spectrum):
    default_spectrum = GetDefaultSpectrum(config, method, meson_name, jet_type, jet_radius, spectrum)
    default_spectrum.SetTitle("Def DMesonJetAnalysis, {0}".format(method))
    wrap.fUseReflections = False
    default_spectrum_from_mt = wrap.GetDefaultSpectrumFromMultiTrial(method)
    default_spectrum_from_mt.SetTitle("Trial Expo Free Sigma, {0}".format(method))
    default_spectrum_from_mt.GetXaxis().SetTitle(xaxisTitle)
    default_spectrum_from_mt.GetYaxis().SetTitle(yaxisTitle)
    globalList.append(default_spectrum)
    globalList.append(default_spectrum_from_mt)
    r = comp.CompareSpectra(default_spectrum_from_mt, [default_spectrum])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

def default_vs_default_mt_unc(comp, method, config, meson_name, jet_type, jet_radius, spectrum):
    default_spectrum = GetDefaultSpectrum(config, method, meson_name, jet_type, jet_radius, spectrum)
    default_spectrum.SetTitle("Def DMesonJetAnalysis, {0}".format(method))
    wrap.fUseReflections = False
    default_spectrum_from_mt = wrap.GetDefaultSpectrumFromMultiTrial(method)
    default_spectrum_from_mt.SetTitle("Trial Expo Free Sigma, {0}".format(method))
    default_spectrum_from_mt.GetXaxis().SetTitle(xaxisTitle)
    default_spectrum_from_mt.GetYaxis().SetTitle(yaxisTitle)
    globalList.append(default_spectrum)
    globalList.append(default_spectrum_from_mt)
    r = comp.CompareUncertainties(default_spectrum_from_mt, [default_spectrum])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

def default_vs_average_raw_yield(comp, method, config, meson_name, jet_type, jet_radius, spectrum):
    default_spectrum = wrap.GetDefaultSpectrumFromMultiTrial(method)
    default_spectrum.SetTitle("Trial Expo Free Sigma, {0}".format(method))
    wrap.fUseReflections = False
    average_spectrum = wrap.GetAverageSpectrumFromMultiTrial(method)
    average_spectrum.SetTitle("Avg Raw Yield Extr Trials, {0}".format(method))
    average_spectrum.GetXaxis().SetTitle(xaxisTitle)
    average_spectrum.GetYaxis().SetTitle(yaxisTitle)
    globalList.append(default_spectrum)
    globalList.append(average_spectrum)
    comp.fRatioRelativeUncertaintyTitle = "Rel. Syst. Unc., {0}".format(method)
    r = comp.CompareSpectra(average_spectrum, [default_spectrum])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

def reflections_raw_yield(comp, method, config, meson_name, jet_type, jet_radius, spectrum):
    wrap.fUseReflections = True
    default_spectrum_wrefl = wrap.GetDefaultSpectrumFromMultiTrial(method)
    default_spectrum_wrefl.SetTitle("Trial Expo Free Sigma w/ refl, {0}".format(method))
    wrap.fUseReflections = False
    default_spectrum_worefl = wrap.GetAverageSpectrumFromMultiTrial(method)
    default_spectrum_worefl.SetTitle("Trial Expo Free Sigma w/o refl, {0}".format(method))
    default_spectrum_worefl.GetXaxis().SetTitle(xaxisTitle)
    default_spectrum_worefl.GetYaxis().SetTitle(yaxisTitle)
    globalList.append(default_spectrum_wrefl)
    globalList.append(default_spectrum_worefl)
    comp.fRatioRelativeUncertaintyTitle = "Rel. Stat. Unc., {0}".format(method)
    r = comp.CompareSpectra(default_spectrum_wrefl, [default_spectrum_worefl])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

def GetDefaultSpectrum(config, method, meson_name, jet_type, jet_radius, spectrum):
    wrap.fUseReflections = False
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fSpectrumName = "_".join([spectrum, method])
    wrap.fDataList = None
    wrap.fDataFile = None
    h = wrap.GetDefaultSpectrumFromDMesonJetAnalysis(method)
    h_copy = h.Clone("{0}_copy".format(spectrum))
    xaxisTitle = h_copy.GetXaxis().GetTitle()
    yaxisTitle = h_copy.GetYaxis().GetTitle()
    return h_copy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='B feed-down.')
    parser.add_argument('yaml', metavar='config.yaml')
    parser.add_argument('--meson', metavar='MESON',
                        default="D0")
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    parser.add_argument('--spectrum', metavar='spectrum',
                        default="JetPtSpectrum_DPt_30")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.meson, args.jet_type, args.jet_radius, args.spectrum)

    IPython.embed()
