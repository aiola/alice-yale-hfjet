#!/usr/bin/env python
# python script to do extract B feed down correction factors

import argparse
import yaml
import IPython
import ROOT
import os
import DMesonJetCompare
import RawYieldSpectrumLoader
import subprocess

globalList = []

wrap = RawYieldSpectrumLoader.RawYieldSpectrumLoader()

xaxisTitle = ""
yaxisTitle = ""
do_spectra_plot = "logy"


def main(config, meson_name, jet_type, jet_radius, var, kincuts):
    global do_spectra_plot

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if var == "JetZ" or var == "Jetz": do_spectra_plot = "lineary"

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    wrap.fInputPath = config["input_path"]
    wrap.fTrain = config["train"]
    wrap.fAnalysisName = config["name"]

    comp = DMesonJetCompare.DMesonJetCompare("AverageRawYieldVsDefault")
    comp.fDoSpectraPlot = do_spectra_plot
    comp.fOptSpectrum = "p x0"
    comp.fOptRatio = ""
    comp.fSeparateBaselineUncertainty = True
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.25
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.25  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = 0.05  # this factor will be used to adjust the y axis in linear scale
    comp.fColors = [ROOT.kRed + 2] * 2
    comp.fMarkers = [ROOT.kOpenSquare, ROOT.kFullSquare]
    comp.fFills = [3254] * 2
    comp.fOptSpectrumBaseline = "e2"
    default_vs_average_raw_yield(comp, "SideBand", config, meson_name, jet_type, jet_radius, var, kincuts)
    if var == "JetPt":
        comp.fColors = [ROOT.kBlue + 2] * 2
        comp.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle]
        comp.fFills = [3245] * 2
        comp.fOptSpectrumBaseline = "e2 same"
        default_vs_average_raw_yield(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, var, kincuts)

    comp = DMesonJetCompare.DMesonJetCompare("DefaultMTRawYieldVsDefault")
    comp.fDoSpectraPlot = do_spectra_plot
    comp.fOptRatio = "hist"
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.20
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.3  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = 0.1  # this factor will be used to adjust the y axis in linear scale
    comp.fNoErrorInBaseline = True
    comp.fColors = [ROOT.kOrange + 2, ROOT.kRed + 2]
    comp.fMarkers = [ROOT.kOpenSquare, ROOT.kFullSquare]
    default_vs_default_mt(comp, "SideBand", config, meson_name, jet_type, jet_radius, var, kincuts, None)
    if var == "JetPt":
        comp.fOptSpectrumBaseline = "same"
        comp.fColors = [ROOT.kGreen + 2, ROOT.kBlue + 2]
        comp.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle]
        default_vs_default_mt(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, var, kincuts, None)

    comp = DMesonJetCompare.DMesonJetCompare("DefaultMTRawYieldVsDefault_DoubleGaus")
    comp.fDoSpectraPlot = do_spectra_plot
    comp.fOptRatio = "hist"
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.20
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.3  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = 0.1  # this factor will be used to adjust the y axis in linear scale
    comp.fNoErrorInBaseline = True
    comp.fColors = [ROOT.kOrange + 2, ROOT.kRed + 2]
    comp.fMarkers = [ROOT.kOpenSquare, ROOT.kFullSquare]
    default_vs_default_mt(comp, "SideBand", config, meson_name, jet_type, jet_radius, var, kincuts, "DoubleGaus")
    if var == "JetPt":
        comp.fOptSpectrumBaseline = "same"
        comp.fColors = [ROOT.kGreen + 2, ROOT.kBlue + 2]
        comp.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle]
        default_vs_default_mt(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, var, kincuts, "DoubleGaus")

    comp = DMesonJetCompare.DMesonJetCompare("DefaultMTRawYieldVsDefaultUncertainties")
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = True
    comp.fOptRatio = "hist"
    comp.fOptSpectrum = "hist"
    comp.fOptSpectrumBaseline = "hist"
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.20
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.3  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = 0.1  # this factor will be used to adjust the y axis in linear scale
    comp.fNoErrorInBaseline = True
    comp.fColors = [ROOT.kRed + 2, ROOT.kRed + 2]
    default_vs_default_mt_unc(comp, "SideBand", config, meson_name, jet_type, jet_radius, var, kincuts)
    if var == "JetPt":
        comp.fOptSpectrumBaseline = "same hist"
        comp.fColors = [ROOT.kBlue + 2, ROOT.kBlue + 2]
        default_vs_default_mt_unc(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, var, kincuts)

    comp = DMesonJetCompare.DMesonJetCompare("SideBandReflectionVariationComparison")
    comp.fDoSpectraPlot = do_spectra_plot
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.20
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.8  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = 0.4  # this factor will be used to adjust the y axis in linear scale
    comp.fGridyRatio = True
    comp.fNoErrorInBaseline = True
    reflections_raw_yield(comp, "SideBand", config, meson_name, jet_type, jet_radius, var, kincuts)

    if var == "JetPt":
        comp = DMesonJetCompare.DMesonJetCompare("InvMassFitReflectionVariationComparison")
        comp.fX1LegRatio = 0.15
        comp.fX1LegSpectrum = 0.20
        comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
        comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
        comp.fLinUpperSpace = 0.8  # this factor will be used to adjust the y axis in linear scale
        comp.fLinLowerSpace = 0.4  # this factor will be used to adjust the y axis in linear scale
        comp.fGridyRatio = True
        comp.fNoErrorInBaseline = True
        reflections_raw_yield(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, var, kincuts)

    SideBandRawYieldUnc(config["input_path"], config["train"], config["name"], var, meson_name, kincuts, "SideBand")
    SideBandRawYieldReflUnc(config["input_path"], config["train"], config["name"], var, meson_name, kincuts, "SideBand")
    SideBandFinalRawYieldUnc(config["input_path"], config["train"], config["name"], var, meson_name, kincuts, "SideBand")

    outputPath = "{}/{}/{}/RawYieldUnc_{}_{}_pdf/{}".format(config["input_path"], config["train"], config["name"], var, kincuts, meson_name)
    if not os.path.isdir(outputPath): os.makedirs(outputPath)
    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            fname = "{0}/{1}.pdf".format(outputPath, obj.GetName())
            print("Saving '{}'".format(obj.GetName()))
            obj.SaveAs(fname)


def default_vs_default_mt(comp, method, config, meson_name, jet_type, jet_radius, var, kincuts, refl):
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fVariableName = var
    wrap.fKinematicCuts = kincuts
    wrap.fRawYieldMethod = method
    if refl:
        wrap.fUseReflections = True
        wrap.fReflectionFit = refl
    else:
        wrap.fUseReflections = False

    wrap.fDataSpectrumList = None
    wrap.fDataFile = None
    default_spectrum = wrap.GetDefaultSpectrumFromDMesonJetAnalysis(False, 0, 0)
    xaxisTitle = default_spectrum.GetXaxis().GetTitle()
    yaxisTitle = default_spectrum.GetYaxis().GetTitle()
    default_spectrum.SetTitle("Def DMesonJetAnalysis, {0}".format(method))

    default_spectrum_from_mt = wrap.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    default_spectrum_from_mt.SetTitle("Trial Expo Free Sigma, {0}".format(method))
    default_spectrum_from_mt.GetXaxis().SetTitle(xaxisTitle)
    default_spectrum_from_mt.GetYaxis().SetTitle(yaxisTitle)

    globalList.append(default_spectrum)
    globalList.append(default_spectrum_from_mt)
    r = comp.CompareSpectra(default_spectrum_from_mt, [default_spectrum])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)


def default_vs_default_mt_unc(comp, method, config, meson_name, jet_type, jet_radius, var, kincuts):
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fVariableName = var
    wrap.fKinematicCuts = kincuts
    wrap.fRawYieldMethod = method
    wrap.fUseReflections = False

    wrap.fDataSpectrumList = None
    wrap.fDataFile = None
    default_spectrum = wrap.GetDefaultSpectrumFromDMesonJetAnalysis(False, 0, 0)
    xaxisTitle = default_spectrum.GetXaxis().GetTitle()
    yaxisTitle = default_spectrum.GetYaxis().GetTitle()
    default_spectrum.SetTitle("Def DMesonJetAnalysis, {0}".format(method))

    default_spectrum_from_mt = wrap.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    default_spectrum_from_mt.SetTitle("Trial Expo Free Sigma, {0}".format(method))
    default_spectrum_from_mt.GetXaxis().SetTitle(xaxisTitle)
    default_spectrum_from_mt.GetYaxis().SetTitle(yaxisTitle)
    globalList.append(default_spectrum)
    globalList.append(default_spectrum_from_mt)
    r = comp.CompareUncertainties(default_spectrum_from_mt, [default_spectrum])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)


def default_vs_average_raw_yield(comp, method, config, meson_name, jet_type, jet_radius, var, kincuts):
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fVariableName = var
    wrap.fKinematicCuts = kincuts
    wrap.fRawYieldMethod = method

    wrap.fUseReflections = True
    wrap.fReflFitFunc = "DoubleGaus"
    wrap.fFixedReflOverSignal = 0
    default_spectrum = wrap.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    default_spectrum.SetTitle("Trial Expo Free Sigma, {0}".format(method))
    average_spectrum = wrap.GetAverageSpectrumFromMultiTrial(False, 0)
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


def reflections_raw_yield(comp, method, config, meson_name, jet_type, jet_radius, var, kincuts):
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fVariableName = var
    wrap.fKinematicCuts = kincuts
    wrap.fRawYieldMethod = method

    histos = []
    wrap.fUseReflections = True
    wrap.fReflectionFit = "DoubleGaus"
    wrap.fReflectionRoS = 0
    default_spectrum_wrefl_DoubleGaus = wrap.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    default_spectrum_wrefl_DoubleGaus.SetTitle("Trial Expo Free Sigma w/ refl double gaus, {0}".format(method))
    globalList.append(default_spectrum_wrefl_DoubleGaus)
    # histos.append(default_spectrum_wrefl_DoubleGaus)

    wrap.fUseReflections = True
    wrap.fReflectionFit = "DoubleGaus"
    wrap.fReflectionRoS = 5
    default_spectrum_wrefl_DoubleGaus_5 = wrap.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    default_spectrum_wrefl_DoubleGaus_5.SetTitle("Trial Expo Free Sigma w/ refl double gaus, 0.5 R/S, {0}".format(method))
    globalList.append(default_spectrum_wrefl_DoubleGaus_5)
    histos.append(default_spectrum_wrefl_DoubleGaus_5)

    wrap.fUseReflections = True
    wrap.fReflectionFit = "DoubleGaus"
    wrap.fReflectionRoS = 15
    default_spectrum_wrefl_DoubleGaus_15 = wrap.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    default_spectrum_wrefl_DoubleGaus_15.SetTitle("Trial Expo Free Sigma w/ refl double gaus, 1.5 R/S, {0}".format(method))
    globalList.append(default_spectrum_wrefl_DoubleGaus_15)
    histos.append(default_spectrum_wrefl_DoubleGaus_15)

    wrap.fUseReflections = True
    wrap.fReflectionFit = "gaus"
    wrap.fReflectionRoS = 0
    default_spectrum_wrefl_gaus = wrap.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    default_spectrum_wrefl_gaus.SetTitle("Trial Expo Free Sigma w/ refl gaus, {0}".format(method))
    globalList.append(default_spectrum_wrefl_gaus)
    histos.append(default_spectrum_wrefl_gaus)

    wrap.fUseReflections = True
    wrap.fReflectionFit = "pol3"
    wrap.fReflectionRoS = 0
    default_spectrum_wrefl_pol3 = wrap.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    default_spectrum_wrefl_pol3.SetTitle("Trial Expo Free Sigma w/ refl pol3, {0}".format(method))
    globalList.append(default_spectrum_wrefl_pol3)
    histos.append(default_spectrum_wrefl_pol3)

    wrap.fUseReflections = True
    wrap.fReflectionFit = "pol6"
    wrap.fReflectionRoS = 0
    default_spectrum_wrefl_pol6 = wrap.GetDefaultSpectrumFromMultiTrial(False, 0, 0)
    default_spectrum_wrefl_pol6.SetTitle("Trial Expo Free Sigma w/ refl pol6, {0}".format(method))
    globalList.append(default_spectrum_wrefl_pol6)
    histos.append(default_spectrum_wrefl_pol6)

    wrap.fUseReflections = False
    default_spectrum_worefl = wrap.GetAverageSpectrumFromMultiTrial(False, 0)
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


def SideBandFinalRawYieldUnc(input_path, train, ana, var, dmeson, kin_cuts, suffix):
    var = var.replace("Z", "z")
    spectrumName = "{}Spectrum".format(var.replace("z", "Z"))
    inputSpectrumName = "_".join([s for s in [dmeson[3:], spectrumName, kin_cuts, suffix] if s])

    fname = "{input_path}/{train}/{ana}/RawYieldUnc_refl_DoubleGaus/{spectrum_name}_DistributionOfFinalYields_SBApproach_{var}_Dzero_AfterDbinSum.root".format(input_path=input_path, train=train, ana=ana, var=var, spectrum_name=inputSpectrumName)
    file = ROOT.TFile(fname)
    cname = "{spectrum_name}_cDistr_{var}_Dzero_SideBand".format(var=var, spectrum_name=inputSpectrumName)
    canvas = file.Get(cname)
    if not isinstance(canvas, ROOT.TCanvas):
        print("Object {} in file {} is not a TCanvas??".format(cname, fname))
        exit(1)
    histos = []
    for obj in canvas.GetListOfPrimitives():
        if isinstance(obj, ROOT.TH1):
            n = len(histos)
            h_copy = obj.Clone("var{0}".format(n))
            h_copy.SetTitle("Variation {0}".format(n))
            if var == "JetPt":
                h_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
            elif var == "Jetz":
                h_copy.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")
            h_copy.GetYaxis().SetTitle("counts")
            histos.append(h_copy)
    file.Close()
    fname = "{input_path}/{train}/{ana}/RawYieldUnc_refl_DoubleGaus/{spectrum_name}_TrialExpoFreeS_{var}_Dzero_SideBand.root".format(input_path=input_path, train=train, ana=ana, var=var, spectrum_name=inputSpectrumName)
    file = ROOT.TFile(fname)
    hname = "f{}SpectrSBDef".format(var)
    baseline = file.Get(hname)
    if not baseline:
        print("Could not find histogram {} in file {}".format(hname, fname))
        exit(1)
    baseline_copy = baseline.Clone("baseline")
    baseline_copy.SetTitle("TrialExpoFreeS")
    if var == "JetPt":
        baseline_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
    elif var == "Jetz":
        baseline_copy.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")
    baseline_copy.GetYaxis().SetTitle("counts")
    file.Close()
    globalList.append(baseline_copy)
    globalList.extend(histos)
    cname = "CompareRawYieldUncVariations_AfterDbinSum"
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fDoSpectraPlot = do_spectra_plot
    comp.fOptRatio = "hist"
    comp.fNColsLegRatio = 3
    comp.fNColsLegSpectrum = 3
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.15
    comp.fDoSpectrumLegend = False
    comp.fDoRatioLegend = False
    r = comp.CompareSpectra(baseline_copy, histos)
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)


def SideBandRawYieldUnc(input_path, train, ana, var, dmeson, kin_cuts, suffix):
    var = var.replace("Z", "z")
    spectrumName = "{}Spectrum".format(var.replace("z", "Z"))
    inputSpectrumName = "_".join([s for s in [dmeson[3:], spectrumName, kin_cuts, suffix] if s])

    for ibin in range(0, 10):
        fname = "{input_path}/{train}/{ana}/RawYieldUnc_refl_DoubleGaus/{spectrum_name}_DistributionOfFinalYields_SBApproach_{var}_Dzero_Bin{ibin}.root".format(input_path=input_path, train=train, ana=ana, var=var, ibin=ibin, spectrum_name=inputSpectrumName)
        if not os.path.isfile(fname): continue
        file = ROOT.TFile(fname)
        canvas = file.Get("cDistr_Dzero_SideBand_{0}".format(ibin))
        file.ls()
        histos = []
        for obj in canvas.GetListOfPrimitives():
            if isinstance(obj, ROOT.TH1):
                n = len(histos)
                h_copy = obj.Clone("bin{0}_var{1}".format(ibin, n))
                h_copy.SetTitle("Variation {0}".format(n))
                if var == "JetPt":
                    h_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
                elif var == "Jetz":
                    h_copy.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")
                h_copy.GetYaxis().SetTitle("counts")
                histos.append(h_copy)
        file.Close()
        fname = "{input_path}/{train}/{ana}/RawYieldUnc_refl_DoubleGaus/{spectrum_name}_TrialExpoFreeS_{var}_Dzero_SideBand_{ibin}.root".format(input_path=input_path, train=train, ana=ana, var=var, ibin=ibin, spectrum_name=inputSpectrumName)
        file = ROOT.TFile(fname)
        hname = "hjet{0}".format(ibin)
        baseline = file.Get(hname)
        if not baseline:
            print("Could not find histogram {} in file {}".format(hname, fname))
            exit(1)
        baseline_copy = baseline.Clone("bin{0}_baseline".format(ibin))
        baseline_copy.SetTitle("TrialExpoFreeS")
        if var == "JetPt":
            baseline_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
        elif var == "Jetz":
            baseline_copy.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")
        baseline_copy.GetYaxis().SetTitle("counts")
        file.Close()
        globalList.append(baseline_copy)
        globalList.extend(histos)
        cname = "CompareRawYieldUncVariations_{0}".format(ibin)
        comp = DMesonJetCompare.DMesonJetCompare(cname)
        comp.fDoSpectraPlot = do_spectra_plot
        comp.fDoSpectrumLegend = False
        comp.fDoRatioLegend = False
        comp.fOptRatio = "hist"
        comp.fNColsLegRatio = 3
        comp.fNColsLegSpectrum = 3
        comp.fX1LegRatio = 0.15
        comp.fX1LegSpectrum = 0.15
        r = comp.CompareSpectra(baseline_copy, histos)
        for obj in r:
            if not obj in globalList:
                globalList.append(obj)


def SideBandRawYieldReflUnc(input_path, train, ana, variable, dmeson, kin_cuts, suffix):
    variable = variable.replace("Z", "z")
    spectrumName = "{}Spectrum".format(variable.replace("z", "Z"))
    inputSpectrumName = "_".join([s for s in [dmeson[3:], spectrumName, kin_cuts, suffix] if s])

    reflVar = ["DoubleGaus_15", "DoubleGaus_5", "gaus", "pol3", "pol6"]
    for ibin in range(0, 10):
        histos = []
        for variation in reflVar:
            fname = "{input_path}/{train}/{ana}/RawYieldUnc_refl_{variation}/{spectrum_name}_TrialExpoFreeS_{variable}_Dzero_SideBand_{ibin}.root".format(input_path=input_path, train=train, ana=ana, variable=variable, ibin=ibin, variation=variation, spectrum_name=inputSpectrumName)
            if not os.path.isfile(fname): continue
            file = ROOT.TFile(fname)
            hname = "hjet{0}".format(ibin)
            h = file.Get(hname)
            if not h:
                print("Could not find histogram {} in file {}".format(hname, fname))
                exit(1)
            h_copy = h.Clone("bin{0}_var{1}".format(ibin, variation))
            h_copy.SetTitle("TrialExpoFreeS, refl={0}".format(variation))
            if variable == "JetPt":
                h_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
            elif variable == "Jetz":
                h_copy.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")
            h_copy.GetYaxis().SetTitle("counts")
            file.Close()
            histos.append(h_copy)
        fname = "{input_path}/{train}/{ana}/RawYieldUnc_refl_DoubleGaus/{spectrum_name}_TrialExpoFreeS_{variable}_Dzero_SideBand_{ibin}.root".format(input_path=input_path, train=train, ana=ana, variable=variable, ibin=ibin, spectrum_name=inputSpectrumName)
        if not os.path.isfile(fname): continue
        file = ROOT.TFile(fname)
        hname = "hjet{0}".format(ibin)
        baseline = file.Get(hname)
        if not baseline:
            print("Could not find histogram {} in file {}".format(hname, fname))
            exit(1)

        baseline_copy = baseline.Clone("bin{0}_baseline".format(ibin))
        baseline_copy.SetTitle("TrialExpoFreeS, refl=DoubleGaus")
        if variable == "JetPt":
            baseline_copy.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
        elif variable == "Jetz":
            baseline_copy.GetXaxis().SetTitle("#it{z}_{||}^{ch jet}")
        baseline_copy.GetYaxis().SetTitle("counts")
        file.Close()
        globalList.append(baseline_copy)
        globalList.extend(histos)
        cname = "CompareRawYieldUncReflVariations_{0}".format(ibin)
        comp = DMesonJetCompare.DMesonJetCompare(cname)
        comp.fDoSpectraPlot = do_spectra_plot
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
                        default="D0_D0toKpiCuts")
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    parser.add_argument('--variable', metavar='VAR',
                        default="JetPt")
    parser.add_argument('--kincuts', metavar='KINCUTS',
                        default="DPt_30")

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.meson, args.jet_type, args.jet_radius, args.variable, args.kincuts)

    IPython.embed()
