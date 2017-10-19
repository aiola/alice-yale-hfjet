#!/usr/bin/env python
# python script analyze the topological cuts

import ROOT
print("ROOT {} imported".format(ROOT.gROOT.GetVersionInt()))
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import yaml
import IPython
import DMesonJetProjectors
import DMesonJetTopoAnalysis
import subprocess

globalList = []

# To mimic ROOT5 behavior
if ROOT.gROOT.GetVersionInt() >= 60000: ROOT.ROOT.Math.IntegratorOneDimOptions.SetDefaultIntegrator("Gauss")

def main(config, maxEvents, format):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    name = config["name"]
    train = config["train"]
    input_path = config["input_path"]
    file_name = config["file_name"]
    output_path = input_path

    projector = DMesonJetProjectors.DMesonJetProjector(input_path, train, file_name, config["task_name"], config["merging_type"], maxEvents)
    ana = DMesonJetTopoAnalysis.DMesonJetTopoAnalysisManager("", "D0", projector)
    globalList.append(ana)

    ana.StartAnalysis()

    ana.SaveRootFile("{0}/{1}".format(output_path, train))
    ana.SavePlots("{0}/{1}".format(output_path, train), format)

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
