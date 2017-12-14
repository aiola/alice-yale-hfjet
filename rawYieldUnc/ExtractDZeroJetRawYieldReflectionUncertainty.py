#!/usr/bin/env python
#  Execute with:
#  ./ExtractDZeroJetRawYieldUncertainty.py config.yaml
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

import ExtractDZeroJetRawYieldUncertainty

import ROOT

globalList = []

ptDbins = [3, 4, 5, 6, 7, 8, 10, 12, 16, 30]
ptJetbins = [5, 6, 8, 10, 14, 20, 30]  # used for eff.scale approach, but also in sideband approach to define the bins of the output jet spectrum


def EvaluateBinPerBinReflUncertainty(config, cuts, dpt_bins_for_eff, dpt_bins, jetpt_bins, binlist_name, sigmafixed, DMesonEff, SBweigth, spectrum_name, spectrum_axis, specie, method, ptmin, ptmax, debug=2):
    # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface = ExtractDZeroJetRawYieldUncertainty.GeneratDzeroJetRawYieldUncSingleTrial(config, cuts, dpt_bins_for_eff, dpt_bins, jetpt_bins, binlist_name, sigmafixed, DMesonEff, SBweigth, spectrum_axis, specie, method, ptmin, ptmax)
    interface.SetSaveInvMassFitCanvases(True)
    interface.SetYieldMethod(method)
    print("Min pt = {0}, max pt = {1}".format(ptmin, ptmax))
    interface.SetPtBinEdgesForMassPlot(float(ptmin), float(ptmax))
    interface.SetFitReflections(True)

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    extract = interface.ExtractInputMassPlot()
    if not extract:
        print("Error in extracting the mass plot! Exiting...")
        exit(1)

    if sum(DMesonEff) < len(DMesonEff) and (SBweigth or method == ROOT.AliDJetRawYieldUncertainty.kEffScale):
        is_eff_corrected = "_efficiency"
    else:
        is_eff_corrected = ""

    if method == ROOT.AliDJetRawYieldUncertainty.kSideband:
        interface.SetSigmaToFixDPtBins(numpy.array(sigmafixed, dtype=numpy.float64))
        # For the reflection templates
        try:
            iBin = dpt_bins.index(ptmin)
        except:
            iBin = -1
        varname = "DPt"
    elif method == ROOT.AliDJetRawYieldUncertainty.kEffScale and spectrum_axis[0] == "jet_pt":
        interface.SetSigmaToFixJetPtBins(numpy.array(sigmafixed, dtype=numpy.float64))
        # For the reflection templates
        try:
            iBin = jetpt_bins.index(ptmin)
        except:
            iBin = -1
        varname = "JetPt"
    else:
        print("EvaluateBinPerBinReflUncertainty: Method {} with spectrum axis '{}' not implemented!".format(method, spectrum_axis[0]))
        exit(1)

    reflFitFuncs = ["gaus", "pol3", "pol6"]
    for reflFitFunc in reflFitFuncs:
        ExtractDZeroJetRawYieldUncertainty.SetReflections(interface, config, cuts, binlist_name, iBin, varname, reflFitFunc, is_eff_corrected)

        multitrial = interface.RunMultiTrial(spectrum_name)
        if not multitrial or not interface.Success():
            print("Error in running the MultiTrial code! Exiting...")
            exit(1)
        globalList.append(multitrial)

        outputPath = "{0}/{1}/{2}/RawYieldUnc_refl_{3}".format(config["input_path"], config["train"], config["name"], reflFitFunc)
        ExtractDZeroJetRawYieldUncertainty.MoveFiles(outputPath, "root")
        ExtractDZeroJetRawYieldUncertainty.MoveFiles(outputPath, "pdf")

    reflFitFunc = "DoubleGaus"
    rovers_factors = [0.5, 1.5]
    for rovers_factor in rovers_factors:
        ExtractDZeroJetRawYieldUncertainty.SetReflections(interface, config, cuts, binlist_name, iBin, varname, reflFitFunc, is_eff_corrected)
        interface.SetValueOfReflOverSignal(-rovers_factor, 1.715, 2.015)  # 1st: ratio of refl/MCsignal (set by hand). If <0: 2nd and 3rd are the range for its evaluation from histo ratios

        multitrial = interface.RunMultiTrial(spectrum_name)
        if not multitrial or not interface.Success():
            print("Error in running the MultiTrial code! Exiting...")
            exit(1)
        globalList.append(multitrial)

        outputPath = "{0}/{1}/{2}/RawYieldUnc_refl_{3}_{4:.0f}".format(config["input_path"], config["train"], config["name"], reflFitFunc, rovers_factor * 10)
        ExtractDZeroJetRawYieldUncertainty.MoveFiles(outputPath, "root")
        ExtractDZeroJetRawYieldUncertainty.MoveFiles(outputPath, "pdf")

    interface.ClearObjects()
    globalList.append(interface)
    return interface


def ExtractDJetRawYieldReflUncertainty(config, cuts, dpt_bins_for_eff, dpt_bins, jetpt_bins, binlist_name, sigmafixed, DMesonEff, SBweigth, spectrum_name, spectrum_axis, specie, method, debug=2):
    interface = ExtractDZeroJetRawYieldUncertainty.GeneratDzeroJetRawYieldUnc(config, cuts, dpt_bins_for_eff, dpt_bins, jetpt_bins, binlist_name, sigmafixed, DMesonEff, SBweigth, spectrum_axis, specie, method)  # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface.SetYieldMethod(method)
    interface.SetMaxNTrialsForSidebandMethod(0)  # only for SB method: number of random trials for each pT(D) bin to build pT(jet) spectrum variations

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    reflFitFuncs = ["gaus", "pol3", "pol6"]
    for reflFitFunc in reflFitFuncs:
        outputPath = "{0}/{1}/{2}/RawYieldUnc_refl_{3}".format(config["input_path"], config["train"], config["name"], reflFitFunc)
        ExtractDZeroJetRawYieldUncertainty.CopyFilesBack(outputPath)
        evalunc = interface.EvaluateUncertainty(spectrum_name)
        if not evalunc:
            print("Error in evaluating the yield uncertainty! Exiting...")
            exit(1)
        ExtractDZeroJetRawYieldUncertainty.MoveFiles(outputPath)

    reflFitFunc = "DoubleGaus"
    rovers_factors = [0.5, 1.5]
    for rovers_factor in rovers_factors:
        outputPath = "{0}/{1}/{2}/RawYieldUnc_refl_{3}_{4:.0f}".format(config["input_path"], config["train"], config["name"], reflFitFunc, rovers_factor * 10)
        ExtractDZeroJetRawYieldUncertainty.CopyFilesBack(outputPath)
        evalunc = interface.EvaluateUncertainty(spectrum_name)
        if not evalunc:
            print("Error in evaluating the yield uncertainty! Exiting...")
            exit(1)
        ExtractDZeroJetRawYieldUncertainty.MoveFiles(outputPath)

    interface.ClearObjects()
    globalList.append(interface)
    return interface


def main(config, b, debug):
    # subprocess.call("make")
    # ROOT.gSystem.Load("AliDJetRawYieldUncertainty.so")

    if b: ROOT.gROOT.SetBatch(True)

    ROOT.gInterpreter.AddIncludePath("$ALICE_ROOT/include");
    ROOT.gInterpreter.AddIncludePath("$ALICE_PHYSICS/include");
    ROOT.gInterpreter.AddIncludePath("$FASTJET/include");

    # load fastjet libraries 3.x
    ROOT.gSystem.Load("libCGAL")

    ROOT.gSystem.Load("libfastjet")
    ROOT.gSystem.Load("libsiscone")
    ROOT.gSystem.Load("libsiscone_spherical")
    ROOT.gSystem.Load("libfastjetplugins")
    ROOT.gSystem.Load("libfastjetcontribfragile")

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    rawYieldUnc = []

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
                    (dpt_bins, DMesonEff) = (ExtractDZeroJetRawYieldUncertainty.GetLimits("d", "fPt", binlist["cuts"]), None)
                    if "efficiency" in binlist:
                        (dpt_bins_for_eff, DMesonEff) = ExtractDZeroJetRawYieldUncertainty.LoadEfficiency(config, binlist["efficiency"], dmeson, "Jet_AKTChargedR040_pt_scheme", dpt_bins, True)

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
                            sigmafixed.pop(0)
                            ibin_prev = ibin
                    method = ROOT.AliDJetRawYieldUncertainty.kSideband
                    spectrum_axis = spectrum["axis"].items()[0]
                    (dpt_bins, DMesonEff) = (binlist_axis[1], None)
                    if "efficiency" in spectrum:
                        (dpt_bins_for_eff, DMesonEff) = ExtractDZeroJetRawYieldUncertainty.LoadEfficiency(config, spectrum["efficiency"], dmeson, "Jet_AKTChargedR040_pt_scheme", dpt_bins, False)
                        for dpt_1, dpt_2 in zip(dpt_bins, dpt_bins_for_eff):
                            if math.fabs(dpt_1 - dpt_2) > 1e-6:
                                print("Detected a mismatch between the efficiency pt bins and the spectrum pt bins")
                                print(dpt_bins)
                                print(dpt_bins_for_eff)
                                exit(1)
                    elif "efficiency" in binlist:
                        (dpt_bins_for_eff, DMesonEff) = ExtractDZeroJetRawYieldUncertainty.LoadEfficiency(config, binlist["efficiency"], dmeson, "Jet_AKTChargedR040_pt_scheme", dpt_bins, True)
                        SBweigth = True

                else:
                    print("Method '{}' not known!".format(spectrum["type"]))
                    exit(1)

                if binlist_axis[0] == "jet_pt":
                    jetpt_bins = binlist_axis[1]
                elif spectrum_axis[0] == "jet_pt":
                    jetpt_bins = spectrum_axis[1]
                else:
                    jetpt_bins = ExtractDZeroJetRawYieldUncertainty.GetLimits("jet", "fPt", binlist["cuts"])
                if not dpt_bins:
                    dpt_bins = ExtractDZeroJetRawYieldUncertainty.GetLimits("d", "fPt", binlist["cuts"])
                if not DMesonEff:
                    DMesonEff = [1.0] * (len(dpt_bins) - 1)

                print("Efficiency: {0}".format(", ".join([str(v) for v in DMesonEff])))
                print("D pt bins for efficiency: {0}".format(", ".join([str(v) for v in dpt_bins_for_eff])))
                print("D pt bins: {0}".format(", ".join([str(v) for v in dpt_bins])))
                print("Jet pt bins: {0}".format(", ".join([str(v) for v in jetpt_bins])))

                for minPt, maxPt in zip(binlist_axis[1][:-1], binlist_axis[1][1:]):
                    interface = EvaluateBinPerBinReflUncertainty(config, cuts, dpt_bins_for_eff, dpt_bins, jetpt_bins, binlist["name"], sigmafixed, DMesonEff, SBweigth, spectrum_name, spectrum_axis, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, method, minPt, maxPt)
                    rawYieldUnc.append(interface)
                rawYieldUncSummary = ExtractDJetRawYieldReflUncertainty(config, cuts, dpt_bins_for_eff, dpt_bins, jetpt_bins, binlist["name"], sigmafixed, DMesonEff, SBweigth, spectrum_name, spectrum_axis, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, method)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract Dzero jet raw yield uncertainties.')
    parser.add_argument('yaml', metavar='file.yaml')
    parser.add_argument('--debug', metavar='debug',
                        default=2)
    parser.add_argument('-b', action='store_const',
                        default=False, const=True)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.b, args.debug)

    IPython.embed()
