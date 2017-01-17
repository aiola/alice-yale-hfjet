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

    fname = "reflTemp/{0}.root".format(config["name"])
    fileOut = ROOT.TFile(fname, "recreate")

    for ptmin, ptmax in zip(ptJetbins[:-1], ptJetbins[1:]):
        hSig = DMesonJetUtils.GetObject(file, "D0_kSignalOnly/Charged_R040/D0_kSignalOnly_Charged_R040_JetPtBins_DPt_30/InvMass_D0_kSignalOnly_JetPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
        if not hSig: exit(1)
        hSig.SetName("histSignal_JetPt_{0:.0f}_{1:.0f}".format(ptmin, ptmax))
        hRefl = DMesonJetUtils.GetObject(file, "D0_WrongPID/Charged_R040/D0_WrongPID_Charged_R040_JetPtBins_DPt_30/InvMass_D0_WrongPID_JetPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
        if not hRefl: exit(1)
        hRefl.SetName("histReflection_JetPt_{0:.0f}_{1:.0f}".format(ptmin, ptmax))
        fileOut.cd()
        hSig.Write()
        hRefl.Write()

    for ptmin, ptmax in zip(ptDbins[:-1], ptDbins[1:]):
        hSig = DMesonJetUtils.GetObject(file, "D0_kSignalOnly/Charged_R040/D0_kSignalOnly_Charged_R040_DPtBins_JetPt_5_30/InvMass_D0_kSignalOnly_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
        if not hSig: exit(1)
        hSig.SetName("histSignal_DPt_{0:.0f}_{1:.0f}".format(ptmin, ptmax))
        hRefl = DMesonJetUtils.GetObject(file, "D0_WrongPID/Charged_R040/D0_WrongPID_Charged_R040_DPtBins_JetPt_5_30/InvMass_D0_WrongPID_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
        if not hRefl: exit(1)
        hRefl.SetName("histReflection_DPt_{0:.0f}_{1:.0f}".format(ptmin, ptmax))
        fileOut.cd()
        hSig.Write()
        hRefl.Write()

    fileOut.Close()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract Dzero jet reflection templates.')
    parser.add_argument('yaml', metavar='file.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config)

    IPython.embed()
