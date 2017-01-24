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

def main(config, meson_name, jet_type, jet_radius):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    wrap.fInputPath = config["input_path"]
    wrap.fTrain = config["train"]
    wrap.fAnalysisName = config["name"]
    spectrum = "JetPtSpectrum"
    kincuts = "DPt_30"

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
    default_vs_average_raw_yield(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, spectrum, kincuts)
    comp.fColors = [ROOT.kRed + 2] * 2
    comp.fMarkers = [ROOT.kOpenSquare, ROOT.kFullSquare]
    comp.fFills = [3254] * 2
    comp.fOptSpectrumBaseline = "e2 same"
    default_vs_average_raw_yield(comp, "SideBand", config, meson_name, jet_type, jet_radius, spectrum, kincuts)

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
    default_vs_default_mt(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, spectrum, kincuts)
    comp.fColors = [ROOT.kOrange + 2, ROOT.kRed + 2]
    comp.fMarkers = [ROOT.kOpenSquare, ROOT.kFullSquare]
    comp.fOptSpectrumBaseline = "same"
    default_vs_default_mt(comp, "SideBand", config, meson_name, jet_type, jet_radius, spectrum, kincuts)

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
    default_vs_default_mt_unc(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, spectrum, kincuts)
    comp.fColors = [ROOT.kRed + 2, ROOT.kRed + 2]
    comp.fOptSpectrumBaseline = "same hist"
    default_vs_default_mt_unc(comp, "SideBand", config, meson_name, jet_type, jet_radius, spectrum, kincuts)

    comp = DMesonJetCompare.DMesonJetCompare("ReflectionVariationComparison")
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.20
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.8  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = 0.4  # this factor will be used to adjust the y axis in linear scale
    comp.fGridyRatio = True
    comp.fNoErrorInBaseline = True
    # comp.fColors = [ROOT.kGreen + 2, ROOT.kBlue + 2]
    # comp.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle]
    # reflections_raw_yield(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, spectrum, kincuts)
    # comp.fColors = [ROOT.kOrange + 2, ROOT.kRed + 2]
    # comp.fMarkers = [ROOT.kOpenSquare, ROOT.kFullSquare]
    # comp.fOptSpectrumBaseline = "same"
    reflections_raw_yield(comp, "SideBand", config, meson_name, jet_type, jet_radius, spectrum, kincuts)

    SideBandRawYieldUnc(config["input_path"], config["train"], config["name"])
    SideBandRawYieldReflUnc(config["input_path"], config["train"], config["name"])
    SideBandFinalRawYieldUnc(config["input_path"], config["train"], config["name"])

    outputPath = "{0}/{1}/{2}/RawYieldUnc_pdf".format(config["input_path"], config["train"], config["name"])
    # if not os.listdir(outputPath): os.makedirs(outputPath)
    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            fname = "{0}/{1}.pdf".format(outputPath, obj.GetName())
            obj.SaveAs(fname)

def default_vs_default_mt(comp, method, config, meson_name, jet_type, jet_radius, spectrum, kincuts):
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fSpectrumName = spectrum
    wrap.fKinematicCuts = kincuts
    wrap.fRawYieldMethod = method

    default_spectrum = GetDefaultSpectrum(config, method, meson_name, jet_type, jet_radius, spectrum, kincuts)
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

def default_vs_default_mt_unc(comp, method, config, meson_name, jet_type, jet_radius, spectrum, kincuts):
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fSpectrumName = spectrum
    wrap.fKinematicCuts = kincuts
    wrap.fRawYieldMethod = method

    default_spectrum = GetDefaultSpectrum(config, method, meson_name, jet_type, jet_radius, spectrum, kincuts)
    default_spectrum.SetTitle("Def DMesonJetAnalysis, {0}".format(method))
    wrap.fUseReflections = False
    default_spectrum_from_mt = wrap.GetDefaultSpectrumFromMultiTrial()
    default_spectrum_from_mt.SetTitle("Trial Expo Free Sigma, {0}".format(method))
    default_spectrum_from_mt.GetXaxis().SetTitle(xaxisTitle)
    default_spectrum_from_mt.GetYaxis().SetTitle(yaxisTitle)
    globalList.append(default_spectrum)
    globalList.append(default_spectrum_from_mt)
    r = comp.CompareUncertainties(default_spectrum_from_mt, [default_spectrum])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

def default_vs_average_raw_yield(comp, method, config, meson_name, jet_type, jet_radius, spectrum, kincuts):
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fSpectrumName = spectrum
    wrap.fKinematicCuts = kincuts
    wrap.fRawYieldMethod = method

    wrap.fUseReflections = True
    wrap.fReflFitFunc = "DoubleGaus"
    wrap.fFixedReflOverSignal = 0
    default_spectrum = wrap.GetDefaultSpectrumFromMultiTrial()
    default_spectrum.SetTitle("Trial Expo Free Sigma, {0}".format(method))
    average_spectrum = wrap.GetAverageSpectrumFromMultiTrial()
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

def reflections_raw_yield(comp, method, config, meson_name, jet_type, jet_radius, spectrum, kincuts):
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fSpectrumName = spectrum
    wrap.fKinematicCuts = kincuts
    wrap.fRawYieldMethod = method

    histos = []
    wrap.fUseReflections = True
    wrap.fReflectionFit = "DoubleGaus"
    wrap.fReflectionRoS = 0
    default_spectrum_wrefl_DoubleGaus = wrap.GetDefaultSpectrumFromMultiTrial()
    default_spectrum_wrefl_DoubleGaus.SetTitle("Trial Expo Free Sigma w/ refl double gaus, {0}".format(method))
    globalList.append(default_spectrum_wrefl_DoubleGaus)
    # histos.append(default_spectrum_wrefl_DoubleGaus)

    wrap.fUseReflections = True
    wrap.fReflectionFit = "DoubleGaus"
    wrap.fReflectionRoS = 5
    default_spectrum_wrefl_DoubleGaus_5 = wrap.GetDefaultSpectrumFromMultiTrial()
    default_spectrum_wrefl_DoubleGaus_5.SetTitle("Trial Expo Free Sigma w/ refl double gaus, 0.5 R/S, {0}".format(method))
    globalList.append(default_spectrum_wrefl_DoubleGaus_5)
    histos.append(default_spectrum_wrefl_DoubleGaus_5)

    wrap.fUseReflections = True
    wrap.fReflectionFit = "DoubleGaus"
    wrap.fReflectionRoS = 15
    default_spectrum_wrefl_DoubleGaus_15 = wrap.GetDefaultSpectrumFromMultiTrial()
    default_spectrum_wrefl_DoubleGaus_15.SetTitle("Trial Expo Free Sigma w/ refl double gaus, 1.5 R/S, {0}".format(method))
    globalList.append(default_spectrum_wrefl_DoubleGaus_15)
    histos.append(default_spectrum_wrefl_DoubleGaus_15)

    wrap.fUseReflections = True
    wrap.fReflectionFit = "gaus"
    wrap.fReflectionRoS = 0
    default_spectrum_wrefl_gaus = wrap.GetDefaultSpectrumFromMultiTrial()
    default_spectrum_wrefl_gaus.SetTitle("Trial Expo Free Sigma w/ refl gaus, {0}".format(method))
    globalList.append(default_spectrum_wrefl_gaus)
    histos.append(default_spectrum_wrefl_gaus)

    wrap.fUseReflections = True
    wrap.fReflectionFit = "pol3"
    wrap.fReflectionRoS = 0
    default_spectrum_wrefl_pol3 = wrap.GetDefaultSpectrumFromMultiTrial()
    default_spectrum_wrefl_pol3.SetTitle("Trial Expo Free Sigma w/ refl pol3, {0}".format(method))
    globalList.append(default_spectrum_wrefl_pol3)
    histos.append(default_spectrum_wrefl_pol3)

    wrap.fUseReflections = True
    wrap.fReflectionFit = "pol6"
    wrap.fReflectionRoS = 0
    default_spectrum_wrefl_pol6 = wrap.GetDefaultSpectrumFromMultiTrial()
    default_spectrum_wrefl_pol6.SetTitle("Trial Expo Free Sigma w/ refl pol6, {0}".format(method))
    globalList.append(default_spectrum_wrefl_pol6)
    histos.append(default_spectrum_wrefl_pol6)

    wrap.fUseReflections = False
    default_spectrum_worefl = wrap.GetAverageSpectrumFromMultiTrial()
    default_spectrum_worefl.SetTitle("Trial Expo Free Sigma w/o refl, {0}".format(method))
    default_spectrum_worefl.GetXaxis().SetTitle(xaxisTitle)
    default_spectrum_worefl.GetYaxis().SetTitle(yaxisTitle)
    globalList.append(default_spectrum_worefl)
    histos.append(default_spectrum_worefl)

    comp.fRatioRelativeUncertaintyTitle = "Rel. Stat. Unc., {0}".format(method)
    comp.fOptRatio = "hist"
    r = comp.CompareSpectra(default_spectrum_wrefl_DoubleGaus, histos)
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

def GetDefaultSpectrum(config, method, meson_name, jet_type, jet_radius, spectrum, kincuts):
    wrap.fUseReflections = False
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fSpectrumName = spectrum
    wrap.fKinematicCuts = kincuts
    wrap.fRawYieldMethod = method
    wrap.fDataSpectrumList = None
    wrap.fDataFile = None
    h = wrap.GetDefaultSpectrumFromDMesonJetAnalysis()
    h_copy = h.Clone("{0}_copy".format(spectrum))
    xaxisTitle = h_copy.GetXaxis().GetTitle()
    yaxisTitle = h_copy.GetYaxis().GetTitle()
    return h_copy

def SideBandFinalRawYieldUnc(input_path, train, ana):
    fname = "{0}/{1}/{2}/RawYieldUnc_refl_DoubleGaus/DistributionOfFinalYields_SBApproach_Dzero_AfterDbinSum.root".format(input_path, train, ana)
    file = ROOT.TFile(fname)
    canvas = file.Get("cDistr_Dzero_SideBand")
    histos = []
    for obj in canvas.GetListOfPrimitives():
        if isinstance(obj, ROOT.TH1):
            n = len(histos)
            h_copy = obj.Clone("var{0}".format(n))
            h_copy.SetTitle("Variation {0}".format(n))
            h_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
            h_copy.GetYaxis().SetTitle("counts")
            histos.append(h_copy)
    file.Close()
    fname = "{0}/{1}/{2}/RawYieldUnc_refl_DoubleGaus/TrialExpoFreeS_Dzero_SideBand.root".format(input_path, train, ana)
    file = ROOT.TFile(fname)
    baseline = file.Get("fJetSpectrSBDef")
    baseline_copy = baseline.Clone("baseline")
    baseline_copy.SetTitle("TrialExpoFreeS")
    baseline_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
    baseline_copy.GetYaxis().SetTitle("counts")
    file.Close()
    globalList.append(baseline_copy)
    globalList.extend(histos)
    cname = "CompareRawYieldUncVariations_AfterDbinSum"
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fOptRatio = "hist"
    comp.fNColsLegRatio = 3
    comp.fNColsLegSpectrum = 3
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.15
    r = comp.CompareSpectra(baseline_copy, histos)
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

def SideBandRawYieldUnc(input_path, train, ana):
    for ibin in range(0, 9):
        fname = "{0}/{1}/{2}/RawYieldUnc_refl_DoubleGaus/DistributionOfFinalYields_SBApproach_Dzero_SideBand_{3}.root".format(input_path, train, ana, ibin)
        file = ROOT.TFile(fname)
        canvas = file.Get("cDistr_Dzero_SideBand_{0}".format(ibin))
        histos = []
        for obj in canvas.GetListOfPrimitives():
            if isinstance(obj, ROOT.TH1):
                n = len(histos)
                h_copy = obj.Clone("bin{0}_var{1}".format(ibin, n))
                h_copy.SetTitle("Variation {0}".format(n))
                h_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
                h_copy.GetYaxis().SetTitle("counts")
                histos.append(h_copy)
        file.Close()
        fname = "{0}/{1}/{2}/RawYieldUnc_refl_DoubleGaus/TrialExpoFreeS_Dzero_SideBand_{3}.root".format(input_path, train, ana, ibin)
        file = ROOT.TFile(fname)
        baseline = file.Get("hjetpt{0}".format(ibin))
        baseline_copy = baseline.Clone("bin{0}_baseline".format(ibin))
        baseline_copy.SetTitle("TrialExpoFreeS")
        baseline_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
        baseline_copy.GetYaxis().SetTitle("counts")
        file.Close()
        globalList.append(baseline_copy)
        globalList.extend(histos)
        cname = "CompareRawYieldUncVariations_{0}".format(ibin)
        comp = DMesonJetCompare.DMesonJetCompare(cname)
        comp.fOptRatio = "hist"
        comp.fNColsLegRatio = 3
        comp.fNColsLegSpectrum = 3
        comp.fX1LegRatio = 0.15
        comp.fX1LegSpectrum = 0.15
        r = comp.CompareSpectra(baseline_copy, histos)
        for obj in r:
            if not obj in globalList:
                globalList.append(obj)

def SideBandRawYieldReflUnc(input_path, train, ana):
    reflVar = ["DoubleGaus_15", "DoubleGaus_5", "gaus", "pol3", "pol6"]
    for ibin in range(0, 9):
        histos = []
        for var in reflVar:
            fname = "{0}/{1}/{2}/RawYieldUnc_refl_{3}/TrialExpoFreeS_Dzero_SideBand_{4}.root".format(input_path, train, ana, var, ibin)
            file = ROOT.TFile(fname)
            h = file.Get("hjetpt{0}".format(ibin))
            h_copy = h.Clone("bin{0}_var{1}".format(ibin, var))
            h_copy.SetTitle("TrialExpoFreeS, refl={0}".format(var))
            h_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
            h_copy.GetYaxis().SetTitle("counts")
            file.Close()
            histos.append(h_copy)
        fname = "{0}/{1}/{2}/RawYieldUnc_refl_DoubleGaus/TrialExpoFreeS_Dzero_SideBand_{3}.root".format(input_path, train, ana, ibin)
        file = ROOT.TFile(fname)
        baseline = file.Get("hjetpt{0}".format(ibin))
        baseline_copy = baseline.Clone("bin{0}_baseline".format(ibin))
        baseline_copy.SetTitle("TrialExpoFreeS, refl=DoubleGaus")
        baseline_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
        baseline_copy.GetYaxis().SetTitle("counts")
        file.Close()
        globalList.append(baseline_copy)
        globalList.extend(histos)
        cname = "CompareRawYieldUncReflVariations_{0}".format(ibin)
        comp = DMesonJetCompare.DMesonJetCompare(cname)
        comp.fOptRatio = "hist"
        comp.fX1LegRatio = 0.30
        comp.fX1LegSpectrum = 0.30
        r = comp.CompareSpectra(baseline_copy, histos)
        for obj in r:
            if not obj in globalList:
                globalList.append(obj)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Raw Yield Uncertainty.')
    parser.add_argument('yaml', metavar='config.yaml')
    parser.add_argument('--meson', metavar='MESON',
                        default="D0")
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.meson, args.jet_type, args.jet_radius)

    IPython.embed()
