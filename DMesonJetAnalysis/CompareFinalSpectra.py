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


def main(config1, config2, meson_name, jet_type, jet_radius, raw, name, no_refl, no_fd, raw_yield_method):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if raw:
        CompareSpectra(config1, config2, meson_name, jet_type, jet_radius, name, no_refl, no_fd, raw_yield_method, GetRawSpectrum)
    else:
        CompareSpectra(config1, config2, meson_name, jet_type, jet_radius, name, no_refl, no_fd, raw_yield_method, GetUnfoldedSpectrum)

def GetUnfoldedSpectrum(config, meson_name, jet_type, jet_radiu, no_refl, no_fd, raw_yield_method):
    if no_refl or no_fd: print("'No refl' and 'No FD' options are ignored!")
    name = "{}_mt_refl_DoubleGaus".format(config["name"])
    fileName = "{0}/{1}/{1}.root".format(config["input_path"], name)
    file = ROOT.TFile(fileName)
    if not file or file.IsZombie():
        print("Could not open file '{}'".format(fileName))
        exit(1)

    spectrumName = "JetPtSpectrum"
    kinematicCuts = "DPt_30"
    rawYieldMethod = raw_yield_method
    analysis_name = "{}_{}".format(rawYieldMethod, kinematicCuts)
    analysis = None
    for analysis in config["analysis"]:
        if analysis["name"] == analysis_name: break
    if not analysis or analysis["name"] != analysis_name:
        print("Could not find analysis '{}'".format(analysis_name))
        exit(1)
    unfoldingMethod = analysis["default_method"]
    prior = analysis["default_prior"]
    reg = analysis["methods"][unfoldingMethod]["default_reg"][prior]
    hname = "{rawYieldMethod}_{kinematicCuts}/{unfoldingMethod}/{rawYieldMethod}_{kinematicCuts}_UnfoldedSpectrum_{unfoldingMethod}_Reg{reg}_Prior{prior}".format(
        rawYieldMethod=rawYieldMethod, kinematicCuts=kinematicCuts, unfoldingMethod=unfoldingMethod, reg=reg, prior=prior)
    h = DMesonJetUtils.GetObject(file, hname)
    if not h:
        print("Could not find histogram '{}'".format(hname))
        exit(1)
    return h

def GetRawSpectrum(config, meson_name, jet_type, jet_radius, no_refl, no_fd, raw_yield_method):
    wrap = RawYieldSpectrumLoader.RawYieldSpectrumLoader()
    wrap.fInputPath = config["input_path"]
    wrap.fTrain = config["train"]
    wrap.fAnalysisName = config["name"]
    wrap.fDMeson = meson_name
    wrap.fJetType = jet_type
    wrap.fJetRadius = jet_radius
    wrap.fSpectrumName = "JetPtSpectrum"
    wrap.fKinematicCuts = "DPt_30"
    wrap.fRawYieldMethod = raw_yield_method
    wrap.fUseReflections = not no_refl
    FDcorr = not no_fd
    return wrap.GetDefaultSpectrumFromMultiTrial(FDcorr)

def CompareSpectra(config1, config2, meson_name, jet_type, jet_radius, name, no_refl, no_fd, raw_yield_method, GetSpectrum):
    h1 = GetSpectrum(config1, meson_name, jet_type, jet_radius, no_refl, no_fd, raw_yield_method)
    h2 = GetSpectrum(config2, meson_name, jet_type, jet_radius, no_refl, no_fd, raw_yield_method)
    globalList.append(h1)
    globalList.append(h2)

    if name:
        cname = name
    else:
        cname = "MyComparison"

    comp = DMesonJetCompare.DMesonJetCompare(cname)
    # comp.fOptSpectrum = "p x0"
    comp.fOptRatio = "hist"
    # comp.fSeparateBaselineUncertainty = True
    # comp.fX1LegRatio = 0.15
    # comp.fX1LegSpectrum = 0.25
    # comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
    # comp.fLogLowerSpace = 2  # this factor will be used to adjust the y axis in log scale
    # comp.fLinUpperSpace = 0.25  # this factor will be used to adjust the y axis in linear scale
    # comp.fLinLowerSpace = 0.05  # this factor will be used to adjust the y axis in linear scale
    # comp.fColors = [ROOT.kBlue + 2] * 2
    # comp.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle]
    # comp.fFills = [3245] * 2
    # comp.fOptSpectrumBaseline = "e2"
    r = comp.CompareSpectra(h1, [h2])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

    if name:
        comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(config["input_path"], comp.fCanvasSpectra.GetName()))
        comp.fCanvasRatio.SaveAs("{0}/{1}.pdf".format(config["input_path"], comp.fCanvasRatio.GetName()))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Raw Yield Uncertainty.')
    parser.add_argument('yaml1',
                        help='First YAML configuration file')
    parser.add_argument('yaml2',
                        help='Second YAML configuration file')
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
    parser.add_argument("--no-refl", action='store_const',
                        default=False, const=True,
                        help='Do not use reflections (only for raw spectra).')
    parser.add_argument("--no-fd", action='store_const',
                        default=False, const=True,
                        help='Do not use B feed-down correction (only for raw spectra).')
    parser.add_argument('--raw-yield-method', metavar='METHOD',
                        default="SideBand",
                        help='Raw yield method')
    args = parser.parse_args()

    f = open(args.yaml1, 'r')
    config1 = yaml.load(f)
    f.close()

    f = open(args.yaml2, 'r')
    config2 = yaml.load(f)
    f.close()

    main(config1, config2, args.meson, args.jet_type, args.jet_radius, args.raw, args.name, args.no_refl, args.no_fd, args.raw_yield_method)

    IPython.embed()
