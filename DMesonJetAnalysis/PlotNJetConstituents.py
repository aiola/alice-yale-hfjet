#!/usr/bin/env python
# python script to prepare B feed-down correction file

import argparse
import IPython
import ROOT
import array
import numpy
import math
import yaml
import DMesonJetUtils

globalList = []

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    fname = "{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"])
    file = ROOT.TFile(fname)
    ptJetbins = [5, 6, 8, 10, 14, 20, 30]
    canvas = ROOT.TCanvas("{0}/{1}/{2}/NJetConstituents".format(config["input_path"], config["train"], config["name"]))
    canvas.Divide(3, 2)
    globalList.append(canvas)
    for i, (ptmin, ptmax) in enumerate(zip(ptJetbins[:-1], ptJetbins[1:])):
        hname = "D0/Charged_R040/D0_Charged_R040_JetPtBins_DPt_30/NJetConstituents_AnyINT_D0_JetPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100)
        h = DMesonJetUtils.GetObject(file, hname)
        canvas.cd(i + 1)
        h_copy = h.DrawCopy()
        h_copy.GetXaxis().SetRangeUser(0, 15)
        globalList.append(h_copy)
    canvas.SaveAs("{0}.pdf".format(canvas.GetName()))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare B feed-down correction file.')
    parser.add_argument('yaml', metavar='conf.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config)

    IPython.embed()
