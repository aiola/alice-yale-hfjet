#!/usr/bin/env python
# python script to run the D meson jet unfolding on 2010 pp data

import argparse
import yaml
import IPython
import DMesonJetUnfolding
import subprocess
import ROOT

globalList = []


def main(config, no_mt, no_refl, refl_fit, refl_ros, fd_syst, ry_syst, format, bg):
    if bg: ROOT.gROOT.SetBatch(True)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gSystem.Load("libRooUnfold")

    ana = DMesonJetUnfolding.DMesonJetUnfolding(config["name"], config["input_path"], config["data_train"], config["data"], config["response_train"], config["response"], not no_mt, not no_refl, refl_fit, refl_ros)

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
    parser.add_argument("--no-mt", action='store_const',
                        default=False, const=True,
                        help='Use results from DMesonJetAnalysis.')
    parser.add_argument('--refl-fit',
                        default="DoubleGaus",
                        help='Use results from Multi-Trial analysis with reflections with fit.')
    parser.add_argument('--refl-ros',
                        default=0, type=int,
                        help='Use results from Multi-Trial analysis with reflections with R/S.')
    parser.add_argument('--no-refl', action='store_const',
                        default=False, const=True,
                        help='Use results without reflections.')
    parser.add_argument("--fd-syst", action='store_const',
                        default=False, const=True,
                        help='Do feed-down systematics.')
    parser.add_argument("--ry-syst", action='store_const',
                        default=False, const=True,
                        help='Do raw-yield systematics.')
    parser.add_argument('-b', action='store_const',
                        default=False, const=True)
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.no_mt, args.no_refl, args.refl_fit, args.refl_ros, args.fd_syst, args.ry_syst, args.format, args.b)

    IPython.embed()
