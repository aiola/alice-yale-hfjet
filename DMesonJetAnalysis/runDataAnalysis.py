#!/usr/bin/env python
# python script to run the D meson jet analysis on 2010 pp data

import os
import argparse
import subprocess
import yaml
import IPython
import DMesonJetProjectors
import DMesonJetAnalysis
import ROOT
print("ROOT {} imported".format(ROOT.gROOT.GetVersionInt()))
ROOT.PyConfig.IgnoreCommandLineOptions = True

globalList = []

# To mimic ROOT5 behavior
if ROOT.gROOT.GetVersionInt() >= 60000: ROOT.ROOT.Math.IntegratorOneDimOptions.SetDefaultIntegrator("Gauss")

def main(config, maxEvents, fmt, gen, proc, ts, stage, ask, bg):
    if bg: ROOT.gROOT.SetBatch(True)

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    if config["train"] == "FastSim" or config["train"] == "FastSimOld":
        suffix = "{0}_{1}_{2}".format(gen, proc, ts)

        if gen == "powheg":
            collision = "POWHEG+PYTHIA6"
        else:
            collision = ""
        if proc == "charm":
            collision += " (c#bar{c})"
        elif proc == "beauty":
            collision += " (b#bar{b})"
        collision += " "
        collision += config["collision_system"]

        name = "{0}_{1}".format(config["name"], suffix)
        output_path = "{0}/FastSim_{1}/".format(config["input_path"], suffix)

        many_ts = ts.split("+")
        if isinstance(stage, str) and "+" in stage:
            many_stage = [int(obj) for obj in stage.split("+")]
        elif (isinstance(stage, str) and stage.isdigit()) or isinstance(stage, (int, long)):
            many_stage = [stage] * len(many_ts)
        else:
            print("Could not parse 'stage' option '{}'".format(stage))
            exit(1)
        input_path = []
        file_name = []
        for iter_ts, iter_stage in zip(many_ts, many_stage):
            iter_suffix = "{0}_{1}_{2}".format(gen, proc, iter_ts)
            iter_input_path = "{0}/FastSim_{1}/".format(config["input_path"], iter_suffix)

            if iter_stage >= 0:
                iter_input_path += "stage_{0}/output".format(iter_stage)
                if config["train"] == "FastSimOld":
                    iter_file_name = "AnalysisResults_FastSim_{0}.root".format(iter_suffix)
                else:
                    iter_file_name = "AnalysisResults_FastSim_{0}_{1}.root".format(gen, proc)
            else:
                iter_input_path += "output"
                iter_file_name = "AnalysisResults_FastSim_{0}_{1}.root".format(gen, proc)
            input_path.append(iter_input_path)
            file_name.append(iter_file_name)

        train = ""
        reflection_templates = None
    else:
        collision = config["collision_system"]
        reflection_templates = config["reflection_templates"]
        name = config["name"]
        train = config["train"]
        input_path = config["input_path"]
        file_name = config["file_name"]
        output_path = input_path

    if "normalization_factor" in config:
        norm_factor = config["normalization_factor"]
    else:
        norm_factor = 1

    if "tree_type" in config:
        tree_type = config["tree_type"]
    else:
        tree_type = "simple"

    if "max_pt_hard" in config:
        max_pt_hard = config["max_pt_hard"]
    else:
        max_pt_hard = -1

    if "reject_outliers" in config and config["reject_outliers"]:
        reject_outliers = dict()
        reject_outliers["outlier_pt_hard_jet_factor"] = config["outlier_pt_hard_jet_factor"]
        reject_outliers["outlier_jet_def"] = config["outlier_jet_def"]
    else:
        reject_outliers = False

    print("The output will be stored in '{}'.".format(output_path))
    if not os.path.isdir(output_path):
        os.makedirs(output_path)

    ana = DMesonJetAnalysis.DMesonJetAnalysis(name)
    projector = DMesonJetProjectors.DMesonJetProjector(input_path, train, file_name, config["task_name"], tree_type, config["merging_type"], norm_factor, maxEvents, max_pt_hard, reject_outliers)
    projector.fDoNotAsk = not ask
    ana.SetProjector(projector)

    for anaConfig in config["analysis"]:
        ana.StartAnalysis(collision, reflection_templates, anaConfig)

    ana.SaveRootFile("{0}/{1}".format(output_path, train))
    ana.SavePlots("{0}/{1}".format(output_path, train), fmt)
    ana.SavePlots("{0}/{1}".format(output_path, train), "C")

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
                        default=-1)
    parser.add_argument('-a', metavar='a', help='Ask before starting each analysis cycle',
                        action='store_const', default=False, const=True)
    parser.add_argument('-b', action='store_const',
                        default=False, const=True)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig, args.events, args.format, args.gen, args.proc, args.ts, args.stage, args.a, args.b)

    IPython.embed()
