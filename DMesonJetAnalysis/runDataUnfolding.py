#!/usr/bin/env python
# python script to run the D meson jet unfolding on 2010 pp data

import argparse
import yaml
import IPython
import DMesonJetUnfolding
import subprocess
import ROOT

globalList = []

def main(config, mt, format):
    # subprocess.call("make")
    # ROOT.gSystem.Load("MassFitter.so")

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    ana = DMesonJetUnfolding.DMesonJetUnfolding(config["name"], config["input_path"], config["data_train"], config["data"], config["response_train"], config["response"], mt)

    for anaConfig in config["analysis"]:
        if anaConfig["active"]:
            ana.StartUnfolding(anaConfig, config["efficiency"], config["use_overflow"])

    ana.SaveRootFile()
    ana.SavePlots(format)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='D meson jet unfolding for 2010 pp data.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--format', metavar='pdf',
                        default="pdf")
    parser.add_argument("--mt", action='store_const',
                        default=False, const=True,
                        help='Use results from Multi-Trial analysis.')
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.mt, args.format)

    IPython.embed()
