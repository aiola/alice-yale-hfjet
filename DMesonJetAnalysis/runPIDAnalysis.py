#!/usr/bin/env python
# python script analyze the PID cuts

import ROOT
print("ROOT {} imported".format(ROOT.gROOT.GetVersionInt()))
ROOT.PyConfig.IgnoreCommandLineOptions = True

import argparse
import yaml
import IPython
import os

import DMesonJetProjectors
import DMesonJetPID

globalList = []


def main(config, maxEvents, format):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    projector = DMesonJetProjectors.DMesonJetProjector(config["input_path"], config["train"], config["file_name"], config["task_name"], config["merging_type"], config["normalization_factor"], maxEvents)

    ana = DMesonJetPID.DMesonJetPIDAnalysis(config["name"], config["title"], config["mc"], config["trigger"], config["d_meson"], config["d_meson_suffix"], config["jet_pt_bins"], config["d_pt_bins"], projector)

    ana.StartAnalysis()

    output_path = "{}/{}/PIDAnalysis".format(config["input_path"], config["train"])
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    ana.SaveRootFile(output_path)
    ana.SavePlots(output_path, format)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Analysis of the topological cuts.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--events', metavar='N',
                        default=-1, type=int)
    parser.add_argument('--format', metavar='pdf',
                        default="pdf")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.events, args.format)

    IPython.embed()
