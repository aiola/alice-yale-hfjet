#!/usr/bin/env python

import argparse
import yaml
import IPython
import numpy

import ROOT

import ExtractDZeroJetRawYieldUncertainty

globalList = []

def main(config, method, minPt, maxPt, refl, no_refl, singleTrial, debug):
    # subprocess.call("make")
    # ROOT.gSystem.Load("AliDJetRawYieldUncertainty.so")

    if no_refl: refl = None

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

    if method == "invmassfit":
        interface = ExtractDZeroJetRawYieldUncertainty.EvaluateBinPerBinUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kEffScale, minPt, maxPt, refl, singleTrial)
    elif method == "sideband":
        interface = ExtractDZeroJetRawYieldUncertainty.EvaluateBinPerBinUncertainty(config, ROOT.AliDJetRawYieldUncertainty.kD0toKpi, ROOT.AliDJetRawYieldUncertainty.kSideband, minPt, maxPt, refl, singleTrial)
    else:
        print("Method {0} unknown!".format(method))

    globalList.append(interface)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract Dzero jet raw yield uncertainties.')
    parser.add_argument('yaml', metavar='file.yaml')
    parser.add_argument('--debug', metavar='debug',
                        default=2)
    parser.add_argument('--ptmin', metavar='pt',
                        default=0, type=float)
    parser.add_argument('--ptmax', metavar='pt',
                        default=0, type=float)
    parser.add_argument('--method', metavar='method',
                        default="invmassfit")
    parser.add_argument('--refl',
                        default="DoubleGaus")
    parser.add_argument('--no-refl', action='store_const',
                        default=False, const=True)
    parser.add_argument('--single-trial', action='store_const',
                        default=False, const=True)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.method, args.ptmin, args.ptmax, args.refl, args.no_refl, args.single_trial, args.debug)

    IPython.embed()
