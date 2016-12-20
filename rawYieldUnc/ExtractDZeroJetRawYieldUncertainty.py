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

import ROOT

globalList = []

ptDbins = [3., 4., 5., 6., 7., 8., 10., 16., 30.]
ptJetbins = [5., 6., 8., 10., 14., 20., 30.]  # used for eff.scale approach, but also in sideband approach to define the bins of the output jet spectrum

def EvaluateBinPerBinUncertainty(config, specie, method, ptmin, ptmax, refl=False, debug=2):
    interface = GeneratDzeroJetRawYieldUnc(config, specie)  # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface.SetYieldMethod(method)
    interface.SetPtBinEdgesForMassPlot(ptmin, ptmax)
    interface.SetFitReflections(refl)

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    extract = interface.ExtractInputMassPlot()
    if not extract:
        print("Error in extracting the mass plot! Exiting...")
        return

    multitrial = interface.RunMultiTrial()
    if not multitrial:
        print("Error in running the MultiTrial code! Exiting...")
        return

    interface.ClearObjects()
    globalList.append(interface)
    return interface

def ExtractDJetRawYieldUncertainty(config, specie, method, nTrials=10, allowRepet=False, debug=2):
    interface = GeneratDzeroJetRawYieldUnc(config, specie)  # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface.SetYieldMethod(method)
    interface.SetMaxNTrialsForSidebandMethod(nTrials)  # only for SB method: number of random trials for each pT(D) bin to build pT(jet) spectrum variations
    interface.SetAllowRepetitionOfTrialExtraction(allowRepet)

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    evalunc = interface.EvaluateUncertainty()
    if not evalunc:
        print("Error in evaluating the yield uncertainty! Exiting...")
        return

    interface.ClearObjects()
    globalList.append(interface)
    return interface

def ExtractDJetRawYieldUncertainty_FromSB_CoherentTrialChoice(config, specie, nTrials=10, debug=2):
    interface = GeneratDzeroJetRawYieldUnc(config, specie)  # here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
    interface.SetYieldMethod(ROOT.AliDJetRawYieldUncertainty.kSideband)
    interface.SetMaxNTrialsForSidebandMethod(nTrials)  # only for SB method: number of random trials for each pT(D) bin to build pT(jet) spectrum variations

    interface.SetDebugLevel(debug)  # 0 = just do the job; 1 = additional printout; 2 = print individual fits

    evalunc = interface.EvaluateUncertainty_CoherentTrialChoice()
    if not evalunc:
        print("Error in evaluating the yield uncertainty! Exiting...")
        return

    interface.ClearObjects()
    globalList.append(interface)
    return interface

def GeneratDzeroJetRawYieldUnc(config, specie, refl):
    # Dzero cfg
    ana = config["analysis"][0]

    # DMesonEff = [0.0118323, 0.02011807,  0.03644752, 0.05664352 ,0.07682878 ,0.08783701, 0.09420746, 0.1047988, 0.1338670, 0.2143196, 0.2574591]
    DMesonEff = [0.05664352 , 0.07682878 , 0.08783701, 0.09420746, 0.1047988, 0.1338670, 0.2143196, 0.2574591]  # chopping 0-1, 1-2, 2-3

    sigmafixed = 0.014
    chi2cut = 3
    meansigmaVar = [True, True, True, True, True, True]  # set mean/sigma variations: fixedS, fixedS+15%, fixedS+15%, freeS&M, freeS/fixedM, fixedS&M
    bkgVar = [True, False, True, False, False, False, False, False]  # set bgk variations: exp, lin, pol2, pol3, pol4, pol5, PowLaw, PowLaw*Exp
    rebinStep = [1]
    minMassStep = [1.72, 1.74]
    maxMassStep = [2.00, 2.03]
    nSigmasBC = [3.5, 4.0]
    # WARNING! set nmask value to active mean/sigma*active bkg variations!
    # And adjust consequently the following matrix (put an entry for each variation, with value: 0=don't consider it, 1=consider it in the final syst eval)
    mask = [
        1, 1,  # fixed sigma (Expo, Lin, Pol2, Pol3, Pol4, Pol5, PowLaw, PowLaw*Exp)
		1, 1,  # fixed sigma+15%
		1, 1,  # fixed sigma-15%
		1, 1,  # free sigma, free mean
		1, 1,  # free sigma, fixed mean
		1, 1  # fixed mean, fixed sigma
        ]

    interface = ROOT.AliDJetRawYieldUncertainty()
    interface.AddInputFileName("{0}/{1}/LHC10b/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    interface.AddInputFileName("{0}/{1}/LHC10c/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    interface.AddInputFileName("{0}/{1}/LHC10d/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    interface.AddInputFileName("{0}/{1}/LHC10e/merge/{2}".format(config["input_path"], config["train"], config["file_name"]))
    interface.SetInputTreename("{0}_{1}_{2}".format(config["task_name"], ana["trigger"][0], ana["d_meson"][0]))
    interface.SetInputDBranchname("DmesonJet")
    interface.SetInputJetBranchname("Jet_AKT{0}{1}_pt_scheme".format(ana["jets"][0]["type"], ana["jets"][0]["radius"]))

    interface.SetMassEdgesAndBinWidthForMassPlot(1.5664, 2.1664, 0.006)
    interface.SetDmesonPtBins(len(ptDbins) - 1, numpy.array(ptDbins))
    interface.SetJetPtBins(len(ptJetbins) - 1, numpy.array(ptJetbins))
    interface.SetDmesonEfficiency(numpy.array(DMesonEff))

    interface.SetSigmaForSignalRegion(2.)  # only for SB method: sigma range of signal region (usually 3 sigma, also 2 is fine if low S/B)
    interface.SetSigmaToFix(sigmafixed)
    interface.SetChi2Cut(chi2cut)
    interface.SetMeanSigmaVariations(numpy.array(meansigmaVar, dtype=bool))
    interface.SetBkgVariations(numpy.array(bkgVar, dtype=bool))
    interface.SetRebinSteps(len(rebinStep), numpy.array(rebinStep, dtype=numpy.int32))
    interface.SetMinMassSteps(len(minMassStep), numpy.array(minMassStep, dtype=numpy.float32))
    interface.SetMaxMassSteps(len(maxMassStep), numpy.array(maxMassStep, dtype=numpy.float32))
    interface.SetSigmaBinCounting(len(nSigmasBC), numpy.array(nSigmasBC, dtype=numpy.float32))
    interface.SetMaskOfVariations(len(mask), numpy.array(mask, dtype=bool))

    interface.SetDmesonSpecie(specie)

    if refl:  # ATTENTION: the histograms to be set are pT-dependent!!
        interface.SetReflFilename("reflections_fitted_DoubleGaus.root")  # file with refl template histo
        interface.SetMCSigFilename("reflections_fitted_DoubleGaus.root")  # file with MC signal histo
        interface.SetReflHistoname("histRflFittedDoubleGaus_ptBin5")  # name of template histo
        interface.SetMCSigHistoname("histSgn_5")  # name of template histo
        interface.SetValueOfReflOverSignal(-1, 1.72, 2.00)  # 1st: ratio of refl/MCsignal (set by hand). If <0: 2nd and 3rd are the range for its evaluation from histo ratios

    return interface

def main(config, debug):
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

    for minPt, maxPt in zip(ptJetbins[:-1], ptJetbins[1:]):
        interface = EvaluateBinPerBinUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kEffScale, minPt, maxPt)
        rawYieldUncInvMassFit.append(interface)
    rawYieldUncSummaryInvMassFit = ExtractDJetRawYieldUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kEffScale)

    for minPt, maxPt in zip(ptDbins[:-1], ptDbins[1:]):
        interface = EvaluateBinPerBinUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kSideband, minPt, maxPt)
        rawYieldUncSideBand.append(interface)
    rawYieldUncSummarySideBand = ExtractDJetRawYieldUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kSideband)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract Dzero jet raw yield uncertainties.')
    parser.add_argument('yaml', metavar='file.yaml')
    parser.add_argument('--debug', metavar='debug',
                        default=2)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.debug)

    IPython.embed()
