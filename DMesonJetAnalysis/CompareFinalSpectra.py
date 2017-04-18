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


def main(configs, meson_name, jet_type, jet_radius, raw, name, no_refl, no_fd, raw_yield_method):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if raw:
        CompareSpectra(configs, meson_name, jet_type, jet_radius, name, no_refl, no_fd, raw_yield_method, GetRawSpectrum)
    else:
        CompareSpectra(configs, meson_name, jet_type, jet_radius, name, no_refl, no_fd, raw_yield_method, GetUnfoldedSpectrum)

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
    h = wrap.GetDefaultSpectrumFromMultiTrial(FDcorr)
    if not wrap.fEvents: wrap.LoadNumberOfEvents()
    h.Scale(1. / wrap.fEvents)
    return h

def CompareSpectra(configs, meson_name, jet_type, jet_radius, name, no_refl, no_fd, raw_yield_method, GetSpectrum):
    histos = []
    for c in configs:
        h = GetSpectrum(c, meson_name, jet_type, jet_radius, no_refl, no_fd, raw_yield_method)
        h.SetTitle(c["name"])
        globalList.append(h)
        histos.append(h)

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
    r = comp.CompareSpectra(histos[0], histos[1:])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

    if name:
        comp.fCanvasSpectra.SaveAs("{0}/{1}.pdf".format(config["input_path"], comp.fCanvasSpectra.GetName()))
        comp.fCanvasRatio.SaveAs("{0}/{1}.pdf".format(config["input_path"], comp.fCanvasRatio.GetName()))

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

    configs = []

    for fname in args.yaml:
        f = open(fname, 'r')
        c = yaml.load(f)
        f.close()
        configs.append(c)

    main(configs, args.meson, args.jet_type, args.jet_radius, args.raw, args.name, args.no_refl, args.no_fd, args.raw_yield_method)

    IPython.embed()
