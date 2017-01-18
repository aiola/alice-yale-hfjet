#!/usr/bin/env python
#  Execute with:
#  ROOT.gSystem.SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ANALYSIS/macros -I$ROOTSYS/include");
#  gROOT->LoadMacro("AliDJetRawYieldUncertainty.cxx++")
#  .L ExtractDJetRawYieldUncertainty.C
#  EvaluateBinPerBinUncertainty(...) //to be done for each pT bin in which you have a mass spectrum
#  ExtractDJetRawYieldUncertainty(...) //to build the uncertainty for the various bins of the jet pT spectrum
#

import argparse
import yaml
import IPython
import numpy
import os
import shutil
import glob

import ROOT

globalList = []

ptDbins = [3, 4, 5, 6, 7, 8, 10, 12, 16, 30]
ptJetbins = [5, 6, 8, 10, 14, 20, 30]  # used for eff.scale approach, but also in sideband approach to define the bins of the output jet spectrum

def EvaluateBinPerBinUncertainty(config, specie, method, ptmin, ptmax, refl=False, debug=2):
    interface = GeneratDzeroJetRawYieldUnc(config, specie, method, ptmin, ptmax, refl)  # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface.SetYieldMethod(method)
    print("Min pt = {0}, max pt = {1}".format(ptmin, ptmax))
    interface.SetPtBinEdgesForMassPlot(float(ptmin), float(ptmax))
    interface.SetFitReflections(refl)

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    extract = interface.ExtractInputMassPlot()
    if not extract:
        print("Error in extracting the mass plot! Exiting...")
        exit(1)

    multitrial = interface.RunMultiTrial()
    if not multitrial:
        print("Error in running the MultiTrial code! Exiting...")
        exit(1)

    interface.ClearObjects()
    globalList.append(interface)
    return interface

def ExtractDJetRawYieldUncertainty(config, specie, method, nTrials=10, allowRepet=False, debug=2):
    interface = GeneratDzeroJetRawYieldUnc(config, specie, method)  # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface.SetYieldMethod(method)
    interface.SetMaxNTrialsForSidebandMethod(nTrials)  # only for SB method: number of random trials for each pT(D) bin to build pT(jet) spectrum variations
    interface.SetAllowRepetitionOfTrialExtraction(allowRepet)

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    evalunc = interface.EvaluateUncertainty()
    if not evalunc:
        print("Error in evaluating the yield uncertainty! Exiting...")
        exit(1)

    interface.ClearObjects()
    globalList.append(interface)
    return interface

def ExtractDJetRawYieldUncertainty_FromSB_CoherentTrialChoice(config, specie, nTrials=10, debug=2):
    interface = GeneratDzeroJetRawYieldUnc(config, specie, method)  # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface.SetYieldMethod(ROOT.AliDJetRawYieldUncertainty.kSideband)
    interface.SetMaxNTrialsForSidebandMethod(nTrials)  # only for SB method: number of random trials for each pT(D) bin to build pT(jet) spectrum variations

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    evalunc = interface.EvaluateUncertainty_CoherentTrialChoice()
    if not evalunc:
        print("Error in evaluating the yield uncertainty! Exiting...")
        exit(1)

    interface.ClearObjects()
    globalList.append(interface)
    return interface

def LoadEfficiency(config):
    try:
        eff_config = config["analysis"][0]["binLists"][0]["efficiency"]
    except:
        eff_config = None
    if not eff_config:
        print("No efficiency requested!")
        return [1.0] * (len(ptDbins) - 1)
    fname = "{0}/{1}".format(config["input_path"], eff_config["file_name"])
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    else:
        print("File {0} successfully open".format(fname))
    rlist = file.Get(eff_config["list_name"])
    if not rlist:
        print("Could not get list {0}".format(eff_config["list_name"]))
        exit(1)
    else:
        print("List {0} successfully open".format(eff_config["list_name"]))
    hist = rlist.FindObject(eff_config["object_name"])
    if not hist:
        print("Could not get histogram {0}".format(eff_config["object_name"]))
        exit(1)
    else:
        print("Histogram {0} successfully open".format(eff_config["object_name"]))
    eff_values = []
    ibinDest = 0
    for ibin in range(1, hist.GetNbinsX() + 1):
        if hist.GetXaxis().GetBinLowEdge(ibin) < ptDbins[0]: continue
        eff_values.append(hist.GetBinContent(ibin))
        print("Copying efficiency {0} from bin {1},{2} to bin {3}, {4}".format(hist.GetBinContent(ibin),
                                                                               hist.GetXaxis().GetBinLowEdge(ibin),
                                                                               hist.GetXaxis().GetBinUpEdge(ibin),
                                                                               ptDbins[ibinDest],
                                                                               ptDbins[ibinDest + 1]))
        ibinDest += 1
        if ibinDest + 1 >= len(ptDbins): break
    return eff_values

def GeneratDzeroJetRawYieldUnc(config, specie, method, ptmin=-1, ptmax=-1, refl=False):
    # Dzero cfg
    ana = config["analysis"][0]

    DMesonEff = LoadEfficiency(config)
    print("Efficiency: {0}".format(", ".join([str(v) for v in DMesonEff])))
    sigmafixed_DPtBins = [0.010, 0.014, 0.016, 0.015, 0.016, 0.015, 0.023, 0.023, 0.027]  # chopping 0-1, 1-2, 2-3
    sigmafixed_JetPtBins = [0.012, 0.015, 0.014, 0.016, 0.018, 0.020]

    chi2cut = 3
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

    interface = ROOT.AliDJetRawYieldUncertainty()
    interface.AddInputFileName("{0}/{1}/LHC10b/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    interface.AddInputFileName("{0}/{1}/LHC10c/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    interface.AddInputFileName("{0}/{1}/LHC10d/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    interface.AddInputFileName("{0}/{1}/LHC10e/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    interface.SetInputTreename("{0}_{1}_{2}".format(config["task_name"], ana["trigger"][0], ana["d_meson"][0]))
    interface.SetInputDBranchname("DmesonJet")
    interface.SetInputJetBranchname("Jet_AKT{0}{1}_pt_scheme".format(ana["jets"][0]["type"], ana["jets"][0]["radius"]))

    interface.SetMassEdgesAndBinWidthForMassPlot(1.565, 2.165, 0.006)
    interface.SetDmesonPtBins(len(ptDbins) - 1, numpy.array(ptDbins, dtype=numpy.float64))
    interface.SetJetPtBins(len(ptJetbins) - 1, numpy.array(ptJetbins, dtype=numpy.float64))
    interface.SetDmesonEfficiency(numpy.array(DMesonEff))

    interface.SetSigmaForSignalRegion(2.)  # only for SB method: sigma range of signal region (usually 3 sigma, also 2 is fine if low S/B)
    interface.SetSigmaToFixDPtBins(numpy.array(sigmafixed_DPtBins, dtype=numpy.float64))
    interface.SetSigmaToFixJetPtBins(numpy.array(sigmafixed_JetPtBins, dtype=numpy.float64))
    interface.SetChi2Cut(chi2cut)
    interface.SetMeanSigmaVariations(numpy.array(meansigmaVar, dtype=bool))
    interface.SetBkgVariations(numpy.array(bkgVar, dtype=bool))
    interface.SetRebinSteps(len(rebinStep), numpy.array(rebinStep, dtype=numpy.int32))
    interface.SetMinMassSteps(len(minMassStep), numpy.array(minMassStep, dtype=numpy.float64))
    interface.SetMaxMassSteps(len(maxMassStep), numpy.array(maxMassStep, dtype=numpy.float64))
    interface.SetSigmaBinCounting(len(nSigmasBC), numpy.array(nSigmasBC, dtype=numpy.float64))
    interface.SetMaskOfVariations(len(mask), numpy.array(mask, dtype=bool))

    interface.SetDmesonSpecie(specie)

    if refl:  # ATTENTION: the histograms to be set are pT-dependent!!
        if method == ROOT.AliDJetRawYieldUncertainty.kEffScale:
            varname = "JetPt"
        elif method == ROOT.AliDJetRawYieldUncertainty.kSideband:
            varname = "DPt"
        interface.SetReflFilename("reflTemp/{0}.root".format(config["reflection_templates"]))  # file with refl template histo
        interface.SetMCSigFilename("reflTemp/{0}.root".format(config["reflection_templates"]))  # file with MC signal histo
        interface.SetReflHistoname("histReflection_{0}_{1:.0f}_{2:.0f}".format(varname, ptmin, ptmax))  # name of template histo
        interface.SetMCSigHistoname("histSignal_{0}_{1:.0f}_{2:.0f}".format(varname, ptmin, ptmax))  # name of template histo
        interface.SetValueOfReflOverSignal(-1, 1.715, 2.015)  # 1st: ratio of refl/MCsignal (set by hand). If <0: 2nd and 3rd are the range for its evaluation from histo ratios

    return interface

def main(config, skip_binbybin, refl, debug):
    # subprocess.call("make")
    # ROOT.gSystem.Load("AliDJetRawYieldUncertainty.so")

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

    ROOT.gROOT.LoadMacro("AliDJetRawYieldUncertainty.cxx+g")

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    rawYieldUncInvMassFit = []
    rawYieldUncSideBand = []

    if not skip_binbybin:
        for minPt, maxPt in zip(ptJetbins[:-1], ptJetbins[1:]):
            interface = EvaluateBinPerBinUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kEffScale, minPt, maxPt, refl)
            rawYieldUncInvMassFit.append(interface)
    rawYieldUncSummaryInvMassFit = ExtractDJetRawYieldUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kEffScale)

    if not skip_binbybin:
        for minPt, maxPt in zip(ptDbins[:-1], ptDbins[1:]):
            interface = EvaluateBinPerBinUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kSideband, minPt, maxPt, refl)
            rawYieldUncSideBand.append(interface)
    rawYieldUncSummarySideBand = ExtractDJetRawYieldUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kSideband)

    outputPath = "{0}/{1}/{2}".format(config["input_path"], config["train"], config["name"])
    MoveFiles(outputPath)

def MoveFiles(outputPath):
    outputPath += "/RawYieldUnc"
    print("Results will be moved to {0}".format(outputPath))
    if not os.path.isdir(outputPath):
        os.makedirs(outputPath)
    for file in glob.glob(r'./*.root'):
        print("Moving file {0}".format(file))
        shutil.copy(file, outputPath)
        os.remove(file)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract Dzero jet raw yield uncertainties.')
    parser.add_argument('yaml', metavar='file.yaml')
    parser.add_argument('--debug', metavar='debug',
                        default=2)
    parser.add_argument('--skip-binbybin', action='store_const',
                        default=False, const=True)
    parser.add_argument('--refl', action='store_const',
                        default=False, const=True)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.skip_binbybin, args.refl, args.debug)

    IPython.embed()
