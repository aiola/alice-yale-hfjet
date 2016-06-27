#!/usr/bin/env python
#python script to run the D meson jet analysis on 2010 pp data

import argparse
import yaml
import IPython
import DMesonJetResponse
import DMesonJetProjectors
import subprocess
import ROOT

globalList = []

def main(config, maxEvents):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    ana = DMesonJetResponse.DMesonJetResponse(config["name"])
    projector = DMesonJetProjectors.DMesonJetResponseProjector(config["input_path"], config["train"], config["file_name"], config["task_name"], maxEvents)
    ana.SetProjector(projector)
    
    for anaConfig in config["analysis"]:
        ana.StartAnalysis(anaConfig)
        ana.SaveRootFile("{0}/{1}".format(config["input_path"], config["train"]))
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='D meson jet response for 2010 pp MC.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('--events', metavar='N',
                        default=-1, type=int)
    args = parser.parse_args()
    
    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.events)
    
    IPython.embed()
    