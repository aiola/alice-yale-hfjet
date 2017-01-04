#!/usr/bin/env python
# python script to run the D meson jet analysis on 2010 pp data

import argparse
import yaml
import IPython
import DMesonJetAnalysis
import DMesonJetProjectors
import subprocess
import ROOT

globalList = []

def main(config, maxEvents, format, gen, proc, ts, stage):

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    suffix = "{0}_{1}_{2}".format(gen, proc, ts)

    if gen == "powheg":
        collision = "POWHEG+PYTHIA6 "
    else:
        collision = ""
    if proc == "charm":
        collision += "(c#bar{c}) "
    elif proc == "beauty":
        collision += "(b#bar{b}) "
    collision += config["collision_system"]

    name = config["name"]
    if "{0}" in name:
        name = name.format(suffix)

    file_name = config["file_name"]
    if "{0}" in file_name:
        file_name = file_name.format(suffix)

    input_path = config["input_path"]
    if "{0}" and "{1}" in input_path:
        input_path = input_path.format(suffix, stage)

    ana = DMesonJetAnalysis.DMesonJetAnalysis(name)
    projector = DMesonJetProjectors.DMesonJetDataProjector(input_path, config["train"], file_name, config["task_name"], config["merging_type"], maxEvents)
    ana.SetProjector(projector)

    for anaConfig in config["analysis"]:
        ana.StartAnalysis(collision, anaConfig)

    ana.SaveRootFile("{0}/{1}".format(input_path, config["train"]))
    ana.SavePlots("{0}/{1}".format(input_path, config["train"]), format)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='D meson jet analysis for 2010 pp data.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--events', metavar='N',
                        default=-1, type=int)
    parser.add_argument('--format', metavar='pdf',
                        default="pdf")
    parser.add_argument('--gen', metavar='GEN',
                        default=None)
    parser.add_argument('--proc', metavar='PROC',
                        default=None)
    parser.add_argument('--ts', metavar='TS',
                        default=None)
    parser.add_argument('--stage', metavar='N',
                        default=1)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.events, args.format, args.gen, args.proc, args.ts, args.stage)

    IPython.embed()
