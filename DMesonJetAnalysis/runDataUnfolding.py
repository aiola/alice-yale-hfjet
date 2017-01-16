#!/usr/bin/env python
# python script to run the D meson jet unfolding on 2010 pp data

import argparse
import yaml
import IPython
import DMesonJetUnfolding
import subprocess
import ROOT

globalList = []

def main(config, mt, fd_syst, ry_syst, format):
    # subprocess.call("make")
    # ROOT.gSystem.Load("MassFitter.so")

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    ana = DMesonJetUnfolding.DMesonJetUnfolding(config["name"], config["input_path"], config["data_train"], config["data"], config["response_train"], config["response"], mt)

    for anaConfig in config["analysis"]:
        if anaConfig["active"]:
            ana.StartUnfolding(anaConfig, config["efficiency"], config["use_overflow"], 0, 0, True)
            if fd_syst:
                for fd_error_band in [-1, 1]: ana.StartUnfolding(anaConfig, config["efficiency"], config["use_overflow"], fd_error_band, 0, False)
            if ry_syst:
                for ry_error_band in [-1, 1]: ana.StartUnfolding(anaConfig, config["efficiency"], config["use_overflow"], 0, ry_error_band, False)

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
    parser.add_argument("--fd-syst", action='store_const',
                        default=False, const=True,
                        help='Do feed-down systematics.')
    parser.add_argument("--ry-syst", action='store_const',
                        default=False, const=True,
                        help='Do raw-yield systematics.')
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.mt, args.fd_syst, args.ry_syst, args.format)

    IPython.embed()
