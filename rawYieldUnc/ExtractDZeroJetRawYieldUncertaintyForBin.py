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

import ExtractDZeroJetRawYieldUncertainty

globalList = []

def main(config, method, minPt, maxPt, debug):
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

    if method == "invmassfit":
        interface = ExtractDZeroJetRawYieldUncertainty.EvaluateBinPerBinUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kEffScale, minPt, maxPt)
    elif method == "sideband":
        interface = ExtractDZeroJetRawYieldUncertainty.EvaluateBinPerBinUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kSideband, minPt, maxPt)
    else:
        print("Method {0} unknown!".format(method))
    globalList.append(interface)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract Dzero jet raw yield uncertainties.')
    parser.add_argument('yaml', metavar='file.yaml')
    parser.add_argument('--debug', metavar='debug',
                        default=2)
    parser.add_argument('--ptmin', metavar='pt',
                        default=0)
    parser.add_argument('--ptmax', metavar='pt',
                        default=0)
    parser.add_argument('--method', metavar='methos',
                        default="invmassfit")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.method, args.ptmin, args.ptmax, args.debug)

    IPython.embed()
