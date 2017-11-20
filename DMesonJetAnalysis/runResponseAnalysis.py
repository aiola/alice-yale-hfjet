#!/usr/bin/env python
# python script to run the D meson jet analysis on 2010 pp data

import argparse
import yaml
import IPython
import DMesonJetResponse
import DMesonJetProjectors
import subprocess
import ROOT

globalList = []

def main(config, maxEvents, suffix, format, doNotAsk, bg):
    if bg: ROOT.gROOT.SetBatch(True)

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    if suffix:
        name = "{0}_{1}".format(config["name"], suffix)
    else:
        name = config["name"]

    ana = DMesonJetResponse.DMesonJetResponse(name)

    if "normalization_factor" in config:
        norm_factor = config["normalization_factor"]
    else:
        norm_factor = 1

    projector = DMesonJetProjectors.DMesonJetProjector(config["input_path"], config["train"], config["file_name"], config["task_name"], config["merging_type"], norm_factor, maxEvents)
    projector.fDoNotAsk = doNotAsk
    ana.SetProjector(projector)

    for anaConfig in config["analysis"]:
        ana.StartAnalysis(anaConfig)

    ana.SaveRootFile("{0}/{1}".format(config["input_path"], config["train"]))
    ana.SavePlots("{0}/{1}".format(config["input_path"], config["train"]), format)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='D meson jet response for 2010 pp MC.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--events', metavar='N',
                        default=-1, type=int)
    parser.add_argument('--suffix', metavar='suffix',
                        default="")
    parser.add_argument('--format', metavar='pdf',
                        default="pdf")
    parser.add_argument('-y', metavar='y', help='Does not stop before starting each analysis cycle',
                        action='store_const', default=False, const=True)
    parser.add_argument('-b', action='store_const',
                        default=False, const=True)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.events, args.suffix, args.format, args.y, args.b)

    IPython.embed()

