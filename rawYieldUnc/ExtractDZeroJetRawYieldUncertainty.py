#!/usr/bin/env python
#  Execute with:
#  ./ExtractDZeroJetRawYieldUncertainty.py config.yaml [--refl]
# For an example of a YAML configuration file, see ../DMesonJetAnalysis/LHC10analysis_Train823.yaml

import argparse
import yaml
import IPython
import numpy
import os
import shutil
import glob
import math
import copy

import ROOT

globalList = []

# To mimic ROOT5 behavior
if ROOT.gROOT.GetVersionInt() >= 60000: ROOT.ROOT.Math.IntegratorOneDimOptions.SetDefaultIntegrator("Gauss")


def EvaluateBinPerBinUncertainty(config, cuts, dpt_bins, jetpt_bins, binlist_name, binlist_axis, sigmafixed, DMesonEff, SBweigth, spectrum_name, spectrum_axis, specie, method, ptmin, ptmax, refl=False, singleTrial=False, debug=2):
    # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    if singleTrial:
        interface = GeneratDzeroJetRawYieldUncSingleTrial(config, cuts, dpt_bins, jetpt_bins, binlist_name, binlist_axis, sigmafixed, DMesonEff, SBweigth, spectrum_axis, specie, method, ptmin, ptmax, refl)
    else:
        interface = GeneratDzeroJetRawYieldUnc(config, cuts, dpt_bins, jetpt_bins, binlist_name, binlist_axis, sigmafixed, DMesonEff, SBweigth, spectrum_axis, specie, method, ptmin, ptmax, refl)
    if singleTrial: interface.SetSaveInvMassFitCanvases(True)
    interface.SetYieldMethod(method)
    print("Min pt = {0}, max pt = {1}".format(ptmin, ptmax))
    interface.SetPtBinEdgesForMassPlot(float(ptmin), float(ptmax))

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    extract = interface.ExtractInputMassPlot()
    if not extract:
        print("Error in extracting the mass plot! Exiting...")
        exit(1)

    multitrial = interface.RunMultiTrial(spectrum_name)
    if not multitrial or not interface.Success():
        print("Error in running the MultiTrial code! Exiting...")
        exit(1)

    interface.ClearObjects()
    globalList.append(multitrial)
    globalList.append(interface)
    return interface


def ExtractDJetRawYieldUncertainty(config, cuts, dpt_bins, jetpt_bins, binlist_name, binlist_axis, sigmafixed, DMesonEff, SBweigth, spectrum_name, spectrum_axis, specie, method, single_trial, nTrials=100, allowRepet=False, debug=2):
    interface = GeneratDzeroJetRawYieldUnc(config, cuts, dpt_bins, jetpt_bins, binlist_name, binlist_axis, sigmafixed, DMesonEff, SBweigth, spectrum_axis, specie, method)  # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface.SetYieldMethod(method)
    # only for SB method: number of random trials for each pT(D) bin to build pT(jet) spectrum variations
    if single_trial: interface.SetMaxNTrialsForSidebandMethod(0)
    else: interface.SetMaxNTrialsForSidebandMethod(nTrials)
    interface.SetAllowRepetitionOfTrialExtraction(allowRepet)

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    evalunc = interface.EvaluateUncertainty(spectrum_name)
    if not evalunc:
        print("Error in evaluating the yield uncertainty! Exiting...")
        exit(1)

    interface.ClearObjects()
    globalList.append(interface)
    return interface


def LoadEfficiency(config, eff_config, dmeson, jetName, dpt_bins):
    if not eff_config:
        print("No efficiency requested!")
        return (dpt_bins, None)

    if dpt_bins:
        print("LoadEfficiency: D pt range provided is {}".format(dpt_bins))

    # This is a temporary hack. The detector response analysis does not have an option for "no jet"
    if not jetName: jetName = "Jet_AKTChargedR040_pt_scheme"

    eff_file_name = "{0}/{1}".format(config["input_path"], eff_config["file_name"])
    eff_list_name = "_".join([obj for obj in ["Prompt", dmeson, jetName, eff_config["list_name"]] if obj])
    eff_obj_name = "_".join([obj for obj in ["Prompt", dmeson, jetName, eff_config["list_name"], eff_config["object_name"]] if obj])

    file = ROOT.TFile(eff_file_name)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(eff_file_name))
        exit(1)
    else:
        print("File {0} successfully open".format(eff_file_name))
    rlist = file.Get(eff_list_name)
    if not rlist:
        print("Could not get list {0}".format(eff_list_name))
        exit(1)
    else:
        print("List {0} successfully open".format(eff_list_name))
    hist = rlist.FindObject(eff_obj_name)
    if not hist:
        print("Could not get histogram {0}".format(eff_obj_name))
        exit(1)
    else:
        print("Histogram {0} successfully open".format(eff_obj_name))
    eff_values = []
    ibinDest = 0
    dpt_bins_dest = []
    for ibin in range(1, hist.GetNbinsX() + 1):
        if dpt_bins:
            if hist.GetXaxis().GetBinLowEdge(ibin) < dpt_bins[0]: continue
            if hist.GetXaxis().GetBinLowEdge(ibin) >= dpt_bins[-1]: continue
        eff_values.append(hist.GetBinContent(ibin))
        dpt_bins_dest.append(hist.GetXaxis().GetBinLowEdge(ibin))
        print("Efficiency {0} for bin {1},{2}".format(hist.GetBinContent(ibin),
                                                      hist.GetXaxis().GetBinLowEdge(ibin),
                                                      hist.GetXaxis().GetBinUpEdge(ibin)))
        ibinDest += 1
    dpt_bins_dest.append(hist.GetXaxis().GetBinUpEdge(ibin))
    return (dpt_bins_dest, eff_values)


def GeneratDzeroJetRawYieldUnc(config, cuts, dpt_bins, jetpt_bins, binlist_name, binlist_axis, sigmafixed, DMesonEff, SBweigth, spectrum_axis, specie, method, ptmin=-1, ptmax=-1, refl=False, reflFitFunc="DoubleGaus"):
    # Dzero cfg
    ana = config["analysis"][0]

    chi2cut = 3000
    meansigmaVar = [True, True, True, True, True, True]  # set mean/sigma variations: fixedS, fixedS+15%, fixedS+15%, freeS&M, freeS/fixedM, fixedS&M
    bkgVar = [True, True, True, False, False, False, False, False]  # set bgk variations: exp, lin, pol2, pol3, pol4, pol5, PowLaw, PowLaw*Exp
    rebinStep = [1, 2]
    minMassStep = [1.715, 1.739, 1.691]
    maxMassStep = [2.015, 1.991, 2.039]
    nSigmasBC = [3.5, 4.0]
    # WARNING! set nmask value to active mean/sigma*active bkg variations!
    # And adjust consequently the following matrix (put an entry for each variation, with value: 0=don't consider it, 1=consider it in the final syst eval)
    mask = [
        1, 1, 1,  # fixed sigma (Expo, Lin, Pol2, Pol3, Pol4, Pol5, PowLaw, PowLaw*Exp)
		1, 1, 1,  # fixed sigma+15%
		1, 1, 1,  # fixed sigma-15%
		1, 1, 1,  # free sigma, free mean
		1, 1, 1,  # free sigma, fixed mean
		1, 1, 1  # fixed mean, fixed sigma
        ]

    reader = ROOT.AliDJetTTreeReader()
    reader.AddInputFileName("{0}/{1}/LHC10b/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    reader.AddInputFileName("{0}/{1}/LHC10c/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    reader.AddInputFileName("{0}/{1}/LHC10d/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    reader.AddInputFileName("{0}/{1}/LHC10e/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    reader.SetInputTreename("{0}_{1}_{2}".format(config["task_name"], ana["trigger"][0], "D0_{}".format(cuts)))
    reader.SetInputDBranchname("DmesonJet")
    reader.SetInputJetBranchname("Jet_AKT{0}{1}_pt_scheme".format(ana["jets"][0]["type"], ana["jets"][0]["radius"]))
    reader.SetMassEdgesAndBinWidthForMassPlot(1.565, 2.165, 0.00006)
    reader.SetMassRebin(100)

    interface = ROOT.AliDJetRawYieldUncertainty()
    interface.SetDJetReader(reader)

    if spectrum_axis[0] == "jet_pt":
        interface.SetJetPtBins(len(spectrum_axis[1]) - 1, numpy.array(spectrum_axis[1], dtype=numpy.float64))
    elif spectrum_axis[0] == "d_z":
        interface.SetJetzBins(len(spectrum_axis[1]) - 1, numpy.array(spectrum_axis[1], dtype=numpy.float64))
    else:
        print("GeneratDzeroJetRawYieldUnc: Axis {} not implemented!".format(spectrum_axis[0]))
        exit(1)
    if binlist_axis[0] == "d_pt":
        interface.SetDmesonPtBins(len(binlist_axis[1]) - 1, numpy.array(binlist_axis[1], dtype=numpy.float64))
        interface.SetSigmaToFixDPtBins(numpy.array(sigmafixed, dtype=numpy.float64))
    elif binlist_axis[0] == "jet_pt":
        interface.SetDmesonPtBins(len(dpt_bins) - 1, numpy.array(dpt_bins, dtype=numpy.float64))
        interface.SetSigmaToFixJetPtBins(numpy.array(sigmafixed, dtype=numpy.float64))
    else:
        print("GeneratDzeroJetRawYieldUnc: Axis {} not implemented!".format(binlist_axis[0]))
        exit(1)

    if spectrum_axis[0] != "jet_pt" and binlist_axis[0] != "jet_pt":
        interface.SetJetPtBins(len(jetpt_bins) - 1, numpy.array(jetpt_bins, dtype=numpy.float64))

    interface.SetDmesonPtBinsForEff(len(dpt_bins) - 1, numpy.array(dpt_bins, dtype=numpy.float64))
    interface.SetDmesonEfficiency(numpy.array(DMesonEff))
    interface.SetUseBkgInBinEdges(False)
    interface.SetSigmaForSignalRegion(2)  # only for SB method: sigma range of signal region (usually 3 sigma, also 2 is fine if low S/B)
    interface.SetSigmaSideBandLeft(8, 4)
    interface.SetSigmaSideBandRight(4, 8)
    interface.SetChi2Cut(chi2cut)
    interface.SetMeanSigmaVariations(numpy.array(meansigmaVar, dtype=bool))
    interface.SetBkgVariations(numpy.array(bkgVar, dtype=bool))
    interface.SetRebinSteps(len(rebinStep), numpy.array(rebinStep, dtype=numpy.int32))
    interface.SetMinMassSteps(len(minMassStep), numpy.array(minMassStep, dtype=numpy.float64))
    interface.SetMaxMassSteps(len(maxMassStep), numpy.array(maxMassStep, dtype=numpy.float64))
    interface.SetSigmaBinCounting(len(nSigmasBC), numpy.array(nSigmasBC, dtype=numpy.float64))
    interface.SetMaskOfVariations(len(mask), numpy.array(mask, dtype=bool))
    interface.SetEfficiencyWeightSB(SBweigth)

    interface.SetDmesonSpecie(specie)
    if sum(DMesonEff) < len(DMesonEff) and (SBweigth or method == ROOT.AliDJetRawYieldUncertainty.kEffScale):
        is_eff_corrected = "_efficiency"
    else:
        is_eff_corrected = ""
    if refl and config["reflection_templates"]:  # ATTENTION: the histograms to be set are pT-dependent!!
        SetReflections(interface, config, cuts, binlist_name, binlist_axis, ptmin, ptmax, reflFitFunc, is_eff_corrected)
    else:
        interface.SetFitReflections(False)
    return interface


def GeneratDzeroJetRawYieldUncSingleTrial(config, cuts, dpt_bins, jetpt_bins, binlist_name, binlist_axis, sigmafixed, DMesonEff, SBweigth, spectrum_axis, specie, method, ptmin=-1, ptmax=-1, refl=False, reflFitFunc="DoubleGaus"):
    # Dzero cfg
    ana = config["analysis"][0]

    chi2cut = 3000
    meansigmaVar = [False, False, False, True, False, False]  # set mean/sigma variations: fixedS, fixedS+15%, fixedS+15%, freeS&M, freeS/fixedM, fixedS&M
    bkgVar = [True, False, False, False, False, False, False, False]  # set bgk variations: exp, lin, pol2, pol3, pol4, pol5, PowLaw, PowLaw*Exp
    rebinStep = [1]
    minMassStep = [1.715]
    maxMassStep = [2.015]
    nSigmasBC = []
    # WARNING! set nmask value to active mean/sigma*active bkg variations!
    # And adjust consequently the following matrix (put an entry for each variation, with value: 0=don't consider it, 1=consider it in the final syst eval)
    mask = [
        1  # free sigma, free mean
        ]

    reader = ROOT.AliDJetTTreeReader()
    reader.AddInputFileName("{0}/{1}/LHC10b/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    reader.AddInputFileName("{0}/{1}/LHC10c/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    reader.AddInputFileName("{0}/{1}/LHC10d/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    reader.AddInputFileName("{0}/{1}/LHC10e/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    reader.SetInputTreename("{0}_{1}_{2}".format(config["task_name"], ana["trigger"][0], "D0_{}".format(cuts)))
    reader.SetInputDBranchname("DmesonJet")
    reader.SetInputJetBranchname("Jet_AKT{0}{1}_pt_scheme".format(ana["jets"][0]["type"], ana["jets"][0]["radius"]))
    reader.SetMassEdgesAndBinWidthForMassPlot(1.565, 2.165, 0.00006)
    reader.SetMassRebin(100)

    interface = ROOT.AliDJetRawYieldUncertainty()
    interface.SetDJetReader(reader)

    if spectrum_axis[0] == "jet_pt":
        interface.SetJetPtBins(len(spectrum_axis[1]) - 1, numpy.array(spectrum_axis[1], dtype=numpy.float64))
    elif spectrum_axis[0] == "d_z":
        interface.SetJetzBins(len(spectrum_axis[1]) - 1, numpy.array(spectrum_axis[1], dtype=numpy.float64))
    else:
        print("GeneratDzeroJetRawYieldUncSingleTrial: Axis {} not implemented!".format(spectrum_axis[0]))
        exit(1)
    if binlist_axis[0] == "d_pt":
        interface.SetDmesonPtBins(len(binlist_axis[1]) - 1, numpy.array(binlist_axis[1], dtype=numpy.float64))
        interface.SetSigmaToFixDPtBins(numpy.array(sigmafixed, dtype=numpy.float64))
    elif binlist_axis[0] == "jet_pt":
        interface.SetDmesonPtBins(len(dpt_bins) - 1, numpy.array(dpt_bins, dtype=numpy.float64))
        interface.SetSigmaToFixJetPtBins(numpy.array(sigmafixed, dtype=numpy.float64))
    else:
        print("GeneratDzeroJetRawYieldUncSingleTrial: Axis {} not implemented!".format(binlist_axis[0]))
        exit(1)

    if spectrum_axis[0] != "jet_pt" and binlist_axis[0] != "jet_pt":
        interface.SetJetPtBins(len(jetpt_bins) - 1, numpy.array(jetpt_bins, dtype=numpy.float64))

    interface.SetDmesonPtBinsForEff(len(dpt_bins) - 1, numpy.array(dpt_bins, dtype=numpy.float64))
    interface.SetDmesonEfficiency(numpy.array(DMesonEff))
    interface.SetUseBkgInBinEdges(False)
    interface.SetSigmaForSignalRegion(2.)  # only for SB method: sigma range of signal region (usually 3 sigma, also 2 is fine if low S/B)
    interface.SetSigmaSideBandLeft(8, 4)
    interface.SetSigmaSideBandRight(4, 8)
    interface.SetChi2Cut(chi2cut)
    interface.SetMeanSigmaVariations(numpy.array(meansigmaVar, dtype=bool))
    interface.SetBkgVariations(numpy.array(bkgVar, dtype=bool))
    interface.SetRebinSteps(len(rebinStep), numpy.array(rebinStep, dtype=numpy.int32))
    interface.SetMinMassSteps(len(minMassStep), numpy.array(minMassStep, dtype=numpy.float64))
    interface.SetMaxMassSteps(len(maxMassStep), numpy.array(maxMassStep, dtype=numpy.float64))
    if len(nSigmasBC) > 0: interface.SetSigmaBinCounting(len(nSigmasBC), numpy.array(nSigmasBC, dtype=numpy.float64))
    interface.SetMaskOfVariations(len(mask), numpy.array(mask, dtype=bool))

    interface.SetDmesonSpecie(specie)
    interface.SetEfficiencyWeightSB(SBweigth)

    if sum(DMesonEff) < len(DMesonEff) and (SBweigth or method == ROOT.AliDJetRawYieldUncertainty.kEffScale):
        is_eff_corrected = "_efficiency"
    else:
        is_eff_corrected = ""

    if refl and config["reflection_templates"]:  # ATTENTION: the histograms to be set are pT-dependent!!
        SetReflections(interface, config, cuts, binlist_name, binlist_axis, ptmin, ptmax, reflFitFunc, is_eff_corrected)
    else:
        interface.SetFitReflections(False)
    return interface


def SetReflections(interface, config, cuts, binlist_name, binlist_axis, ptmin, ptmax, reflFitFunc, is_eff_corrected):
    iBin = binlist_axis[1].index(ptmin)
    if binlist_axis[0] == "d_pt":
        varname = "DPt"
    elif  binlist_axis[0] == "jet_pt":
        varname = "JetPt"
    jet = config["analysis"][0]["jets"][0]
    jet_def = "{}_{}".format(jet["type"], jet["radius"])
    reflFileName = "reflTemp/{refl_name}{is_eff_corrected}_{cuts}_{var}_{jet_def}_{binlist_name}_fitted_{fit}.root".format(refl_name=config["reflection_templates"],
                                                                                                         cuts=cuts, jet_def=jet_def, binlist_name=binlist_name,
                                                                                                         var=varname, fit=reflFitFunc,
                                                                                                         is_eff_corrected=is_eff_corrected)
    interface.SetReflFilename(reflFileName)  # file with refl template histo
    interface.SetMCSigFilename(reflFileName)  # file with MC signal histo
    interface.SetReflHistoname("histRflFitted{fit}_ptBin{bin}".format(fit=reflFitFunc, bin=iBin))  # name of template histo
    interface.SetMCSigHistoname("histSgn_{0}".format(iBin))  # name of template histo
    interface.SetValueOfReflOverSignal(-1, 1.715, 2.015)  # 1st: ratio of refl/MCsignal (set by hand). If <0: 2nd and 3rd are the range for its evaluation from histo ratios
    interface.SetFitReflections(True)


def GetLimits(obj, var, cuts):
    result = None
    for cut in cuts:
        if cut["object"] == obj and cut["variable"] == var:
            if "min" in cut and "max" in cut:
                result = [cut["min"], cut["max"]]
                break
            elif "min" in cut and not "max" in cut:
                if var == "fPt":
                    result = [cut["min"], 30]
                    break
                elif var == "fZ":
                    result = [cut["min"], 1.0001]
                    break
    if not result:
        print("GetLimits called with {}, {}, {}. Error.".format(obj, var, cuts))
        exit(1)

    return result


def main(config, reuse_binbybin, skip_binbybin, skip_combine, single_trial, refl, no_refl, bg, do_not_move, debug):
    if bg: ROOT.gROOT.SetBatch(True)

    if no_refl: refl = None

    # load fastjet libraries 3.x
    ROOT.gSystem.Load("libCGAL")

    ROOT.gSystem.Load("libfastjet")
    ROOT.gSystem.Load("libsiscone")
    ROOT.gSystem.Load("libsiscone_spherical")
    ROOT.gSystem.Load("libfastjetplugins")
    ROOT.gSystem.Load("libfastjetcontribfragile")
    ROOT.gROOT.SetMustClean(False)

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    rawYieldUnc = []

    outputPath = "{0}/{1}/{2}/RawYieldUnc".format(config["input_path"], config["train"], config["name"])
    if refl: outputPath += "_refl_{0}".format(refl)

    if reuse_binbybin: CopyFilesBack(outputPath)

    ana = config["analysis"][0]
    for binlist in ana["binLists"]:
        if not "spectra" in binlist: continue
        if len(binlist["bins"]) != 1:
            print("Cannot process bin lists with a number of axis different from 1 (found {}, bin list name {}).".format(len(binlist["bins"]), binlist["name"]))
            continue
        for spectrum in binlist["spectra"]:
            binlist_axis = copy.deepcopy(binlist["bins"].items()[0])
            if not "multitrial" in spectrum: continue
            for dmeson in spectrum["multitrial"]:
                if not dmeson in ana["d_meson"]: continue
                cuts = dmeson[3:]
                spectrum_name = "{}_{}_{}".format(cuts, spectrum["name"], spectrum["suffix"])
                sigmafixed = binlist["sigma_fits"][dmeson]
                SBweigth = False
                if spectrum["type"] == "inv_mass_fit":
                    if binlist_axis[0] != "jet_pt":
                        print("For the invmassfit method the bin list axis must be jet_pt (it is {}, spectrum {})".format(binlist_axis[0], spectrum["name"]))
                        exit(1)
                    method = ROOT.AliDJetRawYieldUncertainty.kEffScale
                    spectrum_axis = binlist_axis
                    (dpt_bins, DMesonEff) = (GetLimits("d", "fPt", binlist["cuts"]), None)
                    if "efficiency" in binlist:
                        (dpt_bins, DMesonEff) = LoadEfficiency(config, binlist["efficiency"], dmeson, "Jet_AKTChargedR040_pt_scheme", dpt_bins)

                elif spectrum["type"] == "side_band":
                    if len(spectrum["axis"]) != 1:
                        print("Cannot process spectra with a number of axis different from 1 (found {}, spectrum name {}).".format(len(spectrum["axis"]), spectrum["name"]))
                        continue
                    if binlist_axis[0] != "d_pt":
                        print("For the sideband method the bin list axis must be d_pt (it is {}, spectrum {})".format(binlist_axis[0], spectrum["name"]))
                        exit(1)
                    if "skip_bins" in spectrum["side_band"]:
                        ibin_prev = -1
                        for ibin in spectrum["side_band"]["skip_bins"]:
                            if ibin != ibin_prev + 1:
                                print("Can only skip the first x bins!")
                                print(spectrum_name)
                                exit(1)
                            print("Skipping bin {}, {}".format(binlist_axis[1][0], binlist_axis[1][1]))
                            binlist_axis[1].pop(0)
                            ibin_prev = ibin
                    method = ROOT.AliDJetRawYieldUncertainty.kSideband
                    spectrum_axis = spectrum["axis"].items()[0]
                    (dpt_bins, DMesonEff) = (binlist_axis[1], None)
                    if "efficiency" in spectrum:
                        (dpt_bins_eff, DMesonEff) = LoadEfficiency(config, spectrum["efficiency"], dmeson, "Jet_AKTChargedR040_pt_scheme", dpt_bins)
                        for dpt_1, dpt_2 in zip(dpt_bins, dpt_bins_eff):
                            if math.fabs(dpt_1 - dpt_2) > 1e-6:
                                print("Detected a mismatch between the efficiency pt bins and the spectrum pt bins")
                                print(dpt_bins)
                                print(dpt_bins_eff)
                                exit(1)
                    elif "efficiency" in binlist:
                        (dpt_bins, DMesonEff) = LoadEfficiency(config, binlist["efficiency"], dmeson, "Jet_AKTChargedR040_pt_scheme", dpt_bins)
                        SBweigth = True

                else:
                    print("Method '{}' not known!".format(spectrum["type"]))
                    exit(1)

                if binlist_axis[0] == "jet_pt":
                    jetpt_bins = binlist_axis[1]
                elif spectrum_axis[0] == "jet_pt":
                    jetpt_bins = spectrum_axis[1]
                else:
                    jetpt_bins = GetLimits("jet", "fPt", binlist["cuts"])
                if not dpt_bins:
                    dpt_bins = GetLimits("d", "fPt", binlist["cuts"])
                if not DMesonEff:
                    DMesonEff = [1.0] * (len(dpt_bins) - 1)

                print("Efficiency: {0}".format(", ".join([str(v) for v in DMesonEff])))
                print("D pt bins: {0}".format(", ".join([str(v) for v in dpt_bins])))
                print("Jet pt bins: {0}".format(", ".join([str(v) for v in jetpt_bins])))

                if not skip_binbybin and not reuse_binbybin:
                    for minPt, maxPt in zip(binlist_axis[1][:-1], binlist_axis[1][1:]):
                       interface = EvaluateBinPerBinUncertainty(config, cuts, dpt_bins, jetpt_bins, binlist["name"], binlist_axis, sigmafixed, DMesonEff, SBweigth, spectrum_name, spectrum_axis, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, method, minPt, maxPt, refl, single_trial)
                       rawYieldUnc.append(interface)
                if not skip_combine: rawYieldUncSummary = ExtractDJetRawYieldUncertainty(config, cuts, dpt_bins, jetpt_bins, binlist["name"], binlist_axis, sigmafixed, DMesonEff, SBweigth, spectrum_name, spectrum_axis, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, method, single_trial)

    if not do_not_move:
        MoveFiles(outputPath, "root")
        MoveFiles(outputPath, "pdf")


def MoveFiles(outputPath, filetype="root"):
    print("Results will be moved to {0}".format(outputPath))
    if not os.path.isdir(outputPath):
        os.makedirs(outputPath)
    for file in glob.glob("./*.{0}".format(filetype)):
        print("Moving file {0}".format(file))
        shutil.copy(file, outputPath)
        os.remove(file)


def CopyFilesBack(outputPath, filetype="root"):
    print("Results will be copied from {0}".format(outputPath))
    for file in glob.glob("{0}/*.{1}".format(outputPath, filetype)):
        print("Copying file {0}".format(file))
        shutil.copy(file, "./")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract Dzero jet raw yield uncertainties.')
    parser.add_argument('yaml', metavar='file.yaml')
    parser.add_argument('--debug', metavar='debug',
                        default=2)
    parser.add_argument('--skip-binbybin', action='store_const',
                        default=False, const=True)
    parser.add_argument('--reuse-binbybin', action='store_const',
                        default=False, const=True)
    parser.add_argument('--skip-combine', action='store_const',
                        default=False, const=True)
    parser.add_argument('--single-trial', action='store_const',
                        default=False, const=True)
    parser.add_argument('--refl',
                        default="DoubleGaus")
    parser.add_argument('--no-refl', action='store_const',
                        default=False, const=True)
    parser.add_argument('-b', action='store_const',
                        default=False, const=True)
    parser.add_argument('--do-not-move', action='store_const',
                        default=False, const=True)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.reuse_binbybin, args.skip_binbybin, args.skip_combine, args.single_trial, args.refl, args.no_refl, args.b, args.do_not_move, args.debug)

    IPython.embed()
