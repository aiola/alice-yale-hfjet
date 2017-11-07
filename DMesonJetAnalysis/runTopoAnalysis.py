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
import os

globalList = []

# To mimic ROOT5 behavior
if ROOT.gROOT.GetVersionInt() >= 60000: ROOT.ROOT.Math.IntegratorOneDimOptions.SetDefaultIntegrator("Gauss")

def main(configs, maxEvents, format):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    manager = DMesonJetTopoAnalysis.DMesonJetTopoAnalysisManager("D0", configs[0]["topo_studies"]["jet_pt_bins"])
    globalList.append(manager)

    name_chain = "_".join([config["name"] for config in configs])

    for config in configs:
        name = config["name"]
        train = config["train"]
        input_path = config["input_path"]
        file_name = config["file_name"]

        projector = DMesonJetProjectors.DMesonJetProjector(input_path, train, file_name, config["task_name"], config["merging_type"], config["normalization_factor"], maxEvents)

        for topo_ana in config["topo_studies"]["analysis"]:
            manager.AddAnalysis(topo_ana["name"], topo_ana["title"], topo_ana["trigger"], projector, topo_ana["d_meson_suffix"], topo_ana["normalization"])

    manager.StartAnalysis()

    output_path = "{0}/TopoAnalysis_{1}".format(input_path, name_chain)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    manager.SaveRootFile(output_path)
    manager.SavePlots(output_path, format)

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
