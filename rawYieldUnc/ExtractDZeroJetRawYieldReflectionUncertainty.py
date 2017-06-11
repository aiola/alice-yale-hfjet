#!/usr/local/bin/python
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

import ExtractDZeroJetRawYieldUncertainty

import ROOT

globalList = []

ptDbins = [3, 4, 5, 6, 7, 8, 10, 12, 16, 30]
ptJetbins = [5, 6, 8, 10, 14, 20, 30]  # used for eff.scale approach, but also in sideband approach to define the bins of the output jet spectrum

def EvaluateBinPerBinReflUncertainty(config, specie, method, ptmin, ptmax, debug=2):
    # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface = ExtractDZeroJetRawYieldUncertainty.GeneratDzeroJetRawYieldUncSingleTrial(config, specie, method, ptmin, ptmax, False)
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

    if method == ROOT.AliDJetRawYieldUncertainty.kEffScale:
        varname = "JetPt"
        iBin = ptJetbins.index(ptmin)
    elif method == ROOT.AliDJetRawYieldUncertainty.kSideband:
        varname = "DPt"
        iBin = ptDbins.index(ptmin)
    interface.SetMCSigHistoname("histSgn_{0}".format(iBin))  # name of template histo

    reflFitFuncs = ["gaus", "pol3", "pol6"]
    for reflFitFunc in reflFitFuncs:
        interface.SetReflFilename("reflTemp/{0}.root".format(config["reflection_templates"].format(var=varname, fit=reflFitFunc)))  # file with refl template histo
        interface.SetMCSigFilename("reflTemp/{0}.root".format(config["reflection_templates"].format(var=varname, fit=reflFitFunc)))  # file with MC signal histo
        interface.SetReflHistoname("histRflFitted{fit}_ptBin{bin}".format(fit=reflFitFunc, bin=iBin))  # name of template histo
        interface.SetValueOfReflOverSignal(-1, 1.715, 2.015)  # 1st: ratio of refl/MCsignal (set by hand). If <0: 2nd and 3rd are the range for its evaluation from histo ratios

        multitrial = interface.RunMultiTrial()
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
        interface.SetReflFilename("reflTemp/{0}.root".format(config["reflection_templates"].format(var=varname, fit=reflFitFunc)))  # file with refl template histo
        interface.SetMCSigFilename("reflTemp/{0}.root".format(config["reflection_templates"].format(var=varname, fit=reflFitFunc)))  # file with MC signal histo
        interface.SetReflHistoname("histRflFitted{fit}_ptBin{bin}".format(fit=reflFitFunc, bin=iBin))  # name of template histo
        interface.SetValueOfReflOverSignal(-rovers_factor, 1.715, 2.015)  # 1st: ratio of refl/MCsignal (set by hand). If <0: 2nd and 3rd are the range for its evaluation from histo ratios

        multitrial = interface.RunMultiTrial()
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

def ExtractDJetRawYieldReflUncertainty(config, specie, method, debug=2):
    interface = ExtractDZeroJetRawYieldUncertainty.GeneratDzeroJetRawYieldUnc(config, specie, method)  # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface.SetYieldMethod(method)
    interface.SetMaxNTrialsForSidebandMethod(0)  # only for SB method: number of random trials for each pT(D) bin to build pT(jet) spectrum variations

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    reflFitFuncs = ["gaus", "pol3", "pol6"]
    for reflFitFunc in reflFitFuncs:
        outputPath = "{0}/{1}/{2}/RawYieldUnc_refl_{3}".format(config["input_path"], config["train"], config["name"], reflFitFunc)
        ExtractDZeroJetRawYieldUncertainty.CopyFilesBack(outputPath)
        evalunc = interface.EvaluateUncertainty()
        if not evalunc:
            print("Error in evaluating the yield uncertainty! Exiting...")
            exit(1)
        ExtractDZeroJetRawYieldUncertainty.MoveFiles(outputPath)

    reflFitFunc = "DoubleGaus"
    rovers_factors = [0.5, 1.5]
    for rovers_factor in rovers_factors:
        outputPath = "{0}/{1}/{2}/RawYieldUnc_refl_{3}_{4:.0f}".format(config["input_path"], config["train"], config["name"], reflFitFunc, rovers_factor * 10)
        ExtractDZeroJetRawYieldUncertainty.CopyFilesBack(outputPath)
        evalunc = interface.EvaluateUncertainty()
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

    rawYieldUncInvMassFit = []
    rawYieldUncSideBand = []

    for minPt, maxPt in zip(ptJetbins[:-1], ptJetbins[1:]):
        interface = EvaluateBinPerBinReflUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kEffScale, minPt, maxPt)
        rawYieldUncInvMassFit.append(interface)
    rawYieldUncSummaryInvMassFit = ExtractDJetRawYieldReflUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kEffScale)

    for minPt, maxPt in zip(ptDbins[:-1], ptDbins[1:]):
        interface = EvaluateBinPerBinReflUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kSideband, minPt, maxPt)
        rawYieldUncSideBand.append(interface)
    rawYieldUncSummarySideBand = ExtractDJetRawYieldReflUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kSideband)

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
