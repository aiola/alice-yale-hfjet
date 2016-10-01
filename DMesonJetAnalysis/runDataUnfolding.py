#!/usr/bin/env python
#python script to run the D meson jet unfolding on 2010 pp data

import argparse
import yaml
import IPython
import DMesonJetUnfolding
import subprocess
import ROOT

globalList = []

def main(config, format):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    ana = DMesonJetUnfolding.DMesonJetUnfolding(config["name"], config["input_path"], config["data_train"], config["data"], config["response_train"], config["response"])

    for anaConfig in config["analysis"]:
        if anaConfig["active"]:
            ana.StartUnfolding(anaConfig)

    #ana.SaveRootFile("{0}/{1}".format(config["input_path"], config["train"]))
    #ana.SavePlots("{0}/{1}".format(config["input_path"], config["train"]), format)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='D meson jet unfolding for 2010 pp data.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--format', metavar='pdf',
                        default="pdf")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.format)
    
    IPython.embed()
