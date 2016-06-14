#!/usr/bin/env python
#python script to run the D meson jet analysis on 2010 pp data

import argparse
import yaml
import IPython
import DMesonJetAnalysis
import DMesonJetProjectors
import ROOT

globalList = []

def main(config):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    ana = DMesonJetAnalysis.DMesonJetAnalysis(config["name"])
    projector = DMesonJetProjectors.DMesonJetDataProjector(config["input_path"], config["train"], config["file_name"], config["task_name"])
    ana.SetProjector(projector)
    
    for anaConfig in config["analysis"]:
        ana.StartAnalysis(anaConfig)
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='D meson jet analysis for 2010 pp data.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    args = parser.parse_args()
    
    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config)
    
    IPython.embed()
    