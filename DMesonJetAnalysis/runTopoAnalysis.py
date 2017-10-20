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

def main(configs, maxEvents, format):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    manager = DMesonJetTopoAnalysis.DMesonJetTopoAnalysisManager("D0")
    globalList.append(manager)

    name_chain = ""

    for config in configs:
        name = config["name"]
        train = config["train"]
        input_path = config["input_path"]
        file_name = config["file_name"]
        output_path = input_path
        name_chain += name

        projector = DMesonJetProjectors.DMesonJetProjector(input_path, train, file_name, config["task_name"], config["merging_type"], maxEvents)

        for topo_ana in config["topo_studies"]:
            manager.AddAnalysis(topo_ana["title"], topo_ana["trigger"], projector, topo_ana["d_meson_suffix"])

    manager.StartAnalysis()

    manager.SaveRootFile("{0}/{1}".format(output_path, name_chain))
    manager.SavePlots("{0}/{1}".format(output_path, name_chain), format)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Analysis of the topological cuts.')
    parser.add_argument('yaml', metavar='config.yaml', nargs='*',
                        help='YAML configuration file')
    parser.add_argument('--events', metavar='N',
                        default=-1, type=int)
    parser.add_argument('--format', metavar='pdf',
                        default="pdf")
    args = parser.parse_args()

    configs = []

    for yaml_filename in args.yaml:
        f = open(yaml_filename, 'r')
        config = yaml.load(f)
        f.close()
        configs.append(config)

    main(configs, args.events, args.format)

    IPython.embed()
