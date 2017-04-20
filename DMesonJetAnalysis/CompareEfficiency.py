#!/usr/bin/env python
# python script to do extract B feed down correction factors

import argparse
import yaml
import IPython
import ROOT
import os
import DMesonJetCompare
import DMesonJetUtils

globalList = []


def main(configs, meson_name, jet_type, jet_radius, raw, name):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    CompareEfficiency(configs, meson_name, jet_type, jet_radius, name)

def GetEfficiency(config, meson_name, jet_type, jet_radius):
    fileName = "{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"])
    file = ROOT.TFile(fileName)
    if not file or file.IsZombie():
        print("Could not open file '{}'".format(fileName))
        exit(1)

    hname = "{meson_name}_Jet_AKT{jet_type}{jet_radius}_pt_scheme_JetPtDPtSpectrum_CoarseBins/{meson_name}_Jet_AKT{jet_type}{jet_radius}_pt_scheme_JetPtDPtSpectrum_CoarseBins_Efficiency_JetPt_500_3000".format(
        meson_name=meson_name, jet_type=jet_type, jet_radius=jet_radius)
    h = DMesonJetUtils.GetObject(file, hname)
    if not h:
        print("Could not find histogram '{}'".format(hname))
        exit(1)
    return h

def CompareEfficiency(configs, meson_name, jet_type, jet_radius, name):
    histos = []
    for c in configs:
        input_path = c["input_path"]
        h = GetEfficiency(c, meson_name, jet_type, jet_radius)
        h.SetTitle(c["name"])
        globalList.append(h)
        histos.append(h)

    if name:
        cname = name
    else:
        cname = "MyComparison"

    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fOptRatio = "hist"
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.25
    comp.fLogUpperSpace = 5  # this factor will be used to adjust the y axis in log scale
    comp.fLogLowerSpace = 1.5  # this factor will be used to adjust the y axis in log scale
    comp.fLinUpperSpace = 0.45  # this factor will be used to adjust the y axis in linear scale
    comp.fLinLowerSpace = 0.15  # this factor will be used to adjust the y axis in linear scale
    comp.fGridyRatio = True
    r = comp.CompareSpectra(histos[0], histos[1:])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

    if name:
        comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasSpectra.GetName()))
        comp.fCanvasRatio.SaveAs("{0}/{1}.pdf".format(input_path, comp.fCanvasRatio.GetName()))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Raw Yield Uncertainty.')
    parser.add_argument('yaml', nargs='*',
                        help='List of YAML configuration files')
    parser.add_argument('--meson', metavar='MESON',
                        default="D0")
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    parser.add_argument("--raw", action='store_const',
                        default=False, const=True,
                        help='Use raw yield (not unfolded).')
    parser.add_argument('--name', metavar='NAME',
                        default=None)
    args = parser.parse_args()

    configs = []

    for fname in args.yaml:
        f = open(fname, 'r')
        c = yaml.load(f)
        f.close()
        configs.append(c)

    main(configs, args.meson, args.jet_type, args.jet_radius, args.raw, args.name)

    IPython.embed()
