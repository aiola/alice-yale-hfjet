#!/usr/bin/env python

import argparse
import yaml
import IPython

import sys
sys.path.append("../DMesonJetAnalysis")
import DMesonJetUtils
sys.path.remove("../DMesonJetAnalysis")
import ROOT

ptDbins = [3, 4, 5, 6, 7, 8, 10, 12, 16, 30]
ptJetbins = [5, 6, 8, 10, 14, 20, 30]  # used for eff.scale approach, but also in sideband approach to define the bins of the output jet spectrum

def main(config):
    ROOT.TH1.AddDirectory(False)
    fname = "{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"])
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)

    fnameJetPt = "reflTemp/{0}_JetPt.root".format(config["name"])
    fileOutJetPt = ROOT.TFile(fnameJetPt, "recreate")
    for ibin, (ptmin, ptmax) in enumerate(zip(ptJetbins[:-1], ptJetbins[1:])):
        hSig = DMesonJetUtils.GetObject(file, "D0_kSignalOnly/Charged_R040/D0_kSignalOnly_Charged_R040_JetPtBins_DPt_30/InvMass_D0_kSignalOnly_JetPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
        if not hSig: exit(1)
        hSig.SetName("histSgn_{0}".format(ibin))
        hRefl = DMesonJetUtils.GetObject(file, "D0_WrongPID/Charged_R040/D0_WrongPID_Charged_R040_JetPtBins_DPt_30/InvMass_D0_WrongPID_JetPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
        if not hRefl: exit(1)
        hRefl.SetName("histRfl_{0}".format(ibin))
        fileOutJetPt.cd()
        hSig.Write()
        hRefl.Write()
    fileOutJetPt.Close()

    fnameDPt = "reflTemp/{0}_DPt.root".format(config["name"])
    fileOutDPt = ROOT.TFile(fnameDPt, "recreate")
    for ibin, (ptmin, ptmax) in enumerate(zip(ptDbins[:-1], ptDbins[1:])):
        hSig = DMesonJetUtils.GetObject(file, "D0_kSignalOnly/Charged_R040/D0_kSignalOnly_Charged_R040_DPtBins_JetPt_5_30/InvMass_D0_kSignalOnly_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
        if not hSig: exit(1)
        hSig.SetName("histSgn_{0}".format(ibin))
        hRefl = DMesonJetUtils.GetObject(file, "D0_WrongPID/Charged_R040/D0_WrongPID_Charged_R040_DPtBins_JetPt_5_30/InvMass_D0_WrongPID_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
        if not hRefl: exit(1)
        hRefl.SetName("histRfl_{0}".format(ibin))
        fileOutDPt.cd()
        hSig.Write()
        hRefl.Write()
    fileOutDPt.Close()

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

    ROOT.AliDJetRawYieldUncertainty.FitReflDistr(len(ptDbins) - 1, fnameDPt, "gaus");
    ROOT.AliDJetRawYieldUncertainty.FitReflDistr(len(ptJetbins) - 1, fnameJetPt, "gaus");

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract Dzero jet reflection templates.')
    parser.add_argument('yaml', metavar='file.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config)

    IPython.embed()
