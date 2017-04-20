#!/usr/bin/env python
# python script to do extract B feed down correction factors

import argparse
import yaml
import IPython
import ROOT
import os
import DMesonJetCompare
import DMesonJetUtils
import RawYieldSpectrumLoader

globalList = []


def main(configs, meson_name, jet_type, jet_radius, raw, name):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    CompareFeedDown(configs, meson_name, jet_type, jet_radius, name)

def GetFeedDown(config, fd, meson_name, jet_type, jet_radius):
    wrap = RawYieldSpectrumLoader.RawYieldSpectrumLoader()
    wrap.fInputPath = config["input_path"]
    wrap.fTrain = config["train"]
    wrap.fAnalysisName = config["name"]
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fSpectrumName = "JetPtSpectrum"
    wrap.fKinematicCuts = "DPt_30"
    wrap.fRawYieldMethod = "SideBand"
    wrap.fUseReflections = True
    wrap.fFDConfig = fd

    hFDCorr = wrap.GetFDCorrection(0)
    h = wrap.GetDefaultSpectrumFromMultiTrial(False)

    hFDCorr.Divide(h)
    hFDCorr.GetYaxis().SetTitle("FD Fraction")

    return hFDCorr

def CompareFeedDown(configs, meson_name, jet_type, jet_radius, name):
    histos = []
    for c in configs:
        print("Looking for spectrum config in '{}'".format(c["name"]))
        for binList in c["analysis"][0]["binLists"]:
            if not binList["name"] == "JetPtBins_DPt_30":
                print("Skipping bin list '{}'".format(binList["name"]))
                continue
            print("Found bin list '{}'".format(binList["name"]))
            for spectraConfig in binList["spectra"]:
                if spectraConfig["name"] == "JetPtSpectrum_DPt_30":
                    print("Found spectrum '{}'".format(spectraConfig["name"]))
                    break
                else:
                    print("Skipping spectrum '{}'".format(spectraConfig["name"]))
            if spectraConfig and spectraConfig["name"] == "JetPtSpectrum_DPt_30": break
        if not spectraConfig or not spectraConfig["name"] == "JetPtSpectrum_DPt_30":
            print("Error: could not find spectrum 'JetPtSpectrum_DPt_30'")
            exit(1)
        fd = spectraConfig["FD"]
        input_path = c["input_path"]
        h = GetFeedDown(c, fd, meson_name, jet_type, jet_radius)
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
    comp.fDoSpectraPlot = "lineary"
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
