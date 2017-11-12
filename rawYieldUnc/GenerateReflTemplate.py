#!/usr/bin/env python

import argparse
import yaml
import IPython
import sys
import ROOT
import os
import glob
import shutil

sys.path.append("../DMesonJetAnalysis")
import DMesonJetUtils
sys.path.remove("../DMesonJetAnalysis")

def main(config):
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

    ROOT.gSystem.Load("../DMesonJetAnalysis/MassFitter.so")

    ROOT.TH1.AddDirectory(False)
    fname = "{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"])
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)

    path = "{0}/{1}/{2}/reflTemp".format(config["input_path"], config["train"], config["name"])
    if not os.path.isdir(path): os.makedirs(path)

    for templ_config in config["analysis"]:
        GenerateReflTemp(path, config["name"], templ_config)

    dest_dir = "./reflTemp"
    for file in glob.glob("{}/*.root".format(path)):
        print("Copying file '{}' to '{}'".format(file, dest_dir))
        shutil.copy(file, dest_dir)

def GenerateReflTemp(path, name, templ_config):
    for templ in templ_config["templates"]:
        if jet_label:
            fname = "{path}/{name}_{cuts}_{variable}_{jet_label}_{bin_list}.root".format(path=path, name=name, variable=templ["variable"], jet_label=templ["jet_label"], bin_list=templ["bin_list"])
        else:
            fname = "{path}/{name}_{cuts}_{variable}_{bin_list}.root".format(path=path, name=name, variable=templ["variable"], bin_list=templ["bin_list"])
        fileOut = ROOT.TFile(fname, "recreate")

        for ibin, (ptmin, ptmax) in enumerate(zip(templ_config["bins"][:-1], templ_config["bins"][1:])):
            if jet_label:
                obj_sig_name = "D0_kSignalOnly_{cuts}/{jet_label}/D0_kSignalOnly_{cuts}_{jet_label}_{bin_list}/InvMass_D0_kSignalOnly_{cuts}_DPt_{minpt:.0f}_{maxpt:.0f}".format(jet_label=templ["jet_label"], bin_list=templ["bin_list"], minpt=ptmin * 100, maxpt=ptmax * 100)
                obj_refl_name = "D0_WrongPID_{cuts}/{jet_label}/D0_WrongPID_{cuts}_{jet_label}_{bin_list}/InvMass_D0_WrongPID_{cuts}_DPt_{minpt:.0f}_{maxpt:.0f}".format(jet_label=templ["jet_label"], bin_list=templ["bin_list"], minpt=ptmin * 100, maxpt=ptmax * 100)
            else:
                obj_sig_name = "D0_kSignalOnly_{cuts}/D0_kSignalOnly_{cuts}_{bin_list}/InvMass_D0_kSignalOnly_{cuts}_DPt_{minpt:.0f}_{maxpt:.0f}".format(bin_list=templ["bin_list"], minpt=ptmin * 100, maxpt=ptmax * 100)
                obj_refl_name = "D0_WrongPID_{cuts}/D0_WrongPID_{cuts}_{bin_list}/InvMass_D0_WrongPID_{cuts}_DPt_{minpt:.0f}_{maxpt:.0f}".format(bin_list=templ["bin_list"], minpt=ptmin * 100, maxpt=ptmax * 100)

            hSig = DMesonJetUtils.GetObject(file, obj_sig_name)
            if not hSig: exit(1)
            hSig.SetName("histSgn_{0}".format(ibin))
            # sigInt = hSig.Integral(hSig.GetXaxis().FindBin(1.715), hSig.GetXaxis().FindBin(2.015))

            hRefl = DMesonJetUtils.GetObject(file, obj_refl_name)
            if not hRefl: exit(1)
            hRefl.SetName("histRfl_{0}".format(ibin))
            # refInt = hSig.Integral(hRefl.GetXaxis().FindBin(1.715), hRefl.GetXaxis().FindBin(2.015))

            hSig.Scale(1. / hRefl.Integral())
            hRefl.Scale(1. / hRefl.Integral())
            fileOutJetPt.cd()
            hSig.Write()
            hRefl.Write()
        fileOut.Close()
        ROOT.AliDJetRawYieldUncertainty.FitReflDistr(len(templ_config["bins"]) - 1, fname, "DoubleGaus");
        ROOT.AliDJetRawYieldUncertainty.FitReflDistr(len(templ_config["bins"]) - 1, fname, "gaus");
        ROOT.AliDJetRawYieldUncertainty.FitReflDistr(len(templ_config["bins"]) - 1, fname, "pol3");
        ROOT.AliDJetRawYieldUncertainty.FitReflDistr(len(templ_config["bins"]) - 1, fname, "pol6");

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract Dzero jet reflection templates.')
    parser.add_argument('yaml', metavar='file.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config)

    IPython.embed()
