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


def main(config, meson_name, jet_type, jet_radius, name, refl, no_fd, raw_yield_method, inclusive_train):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    centralityBins = [0, 10, 30, 50]
    for cmin, cmax in zip(centralityBins[:-1], centralityBins[1:]):
        CompareSpectra(cmin, cmax, config, meson_name, jet_type, jet_radius, name, refl, no_fd, raw_yield_method, inclusive_train)

def GetRawDJetSpectrum(trigger, config, meson_name, jet_type, jet_radius, refl, no_fd, raw_yield_method):
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
    wrap.fUseReflections = refl
    wrap.fTrigger = trigger
    FDcorr = not no_fd

    if FDcorr:
        print("Looking for spectrum config in '{}'".format(config["name"]))
        for binList in config["analysis"][0]["binLists"]:
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
        wrap.fFDConfig = spectraConfig["FD"]

    h = wrap.GetDefaultSpectrumFromDMesonJetAnalysis(FDcorr)
    if not wrap.fEvents: wrap.LoadNumberOfEvents()
    h.Scale(1. / wrap.fEvents)
    h.GetYaxis().SetTitle("per-event raw yield")
    return h

def GetInclusiveJets(input_path, inclusive_train, cmin, cmax):
    histos = []
    fname = "{}/{}/Jets.root".format(input_path, inclusive_train)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not find file '{}'".format(fname))
        exit(1)
    leadPtBins = [0, 3, 5]
    for leadPt in leadPtBins:
        rlist = file.Get(str(((cmin, cmax), leadPt)))
        if not rlist:
            print("Could not find list '{}'".format(str((cmin, cmax), leadPt)))
            exit(1)
        hname = "JetPt_Cent{}_{}_LeadHadPt{}".format(cmin, cmax, leadPt)
        h = rlist.FindObject(hname)
        if not h:
            print("Could not find histogram '{}'".format(hname))
            exit(1)
        if leadPt > 0:
            h.SetTitle("#it{{p}}_{{T,leading}} > {} GeV/#it{{c}}".format(leadPt))
        else:
            h.SetTitle("Inclusive")
        h.Sumw2()
        h.Scale(1., "width")
        histos.append(h)
    return histos

def CompareSpectra(cmin, cmax, config, meson_name, jet_type, jet_radius, name, no_refl, no_fd, raw_yield_method, inclusive_train):
    histos = []

    histos.extend(GetInclusiveJets(config["input_path"], inclusive_train, cmin, cmax))

    trigger = "INT7_Cent_{}_{}".format(cmin, cmax)
    h = GetRawDJetSpectrum(trigger, config, meson_name, jet_type, jet_radius, no_refl, no_fd, raw_yield_method)
    h.SetTitle("#it{p}_{T,D^{0}} > 3 GeV/#it{c} (x 10^{7})")
    h.Scale(1e7, "width")
    globalList.append(h)
    histos.append(h)

    if name:
        cname = name
    else:
        cname = "MyComparison_{}_{}".format(cmin, cmax)

    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fOptRatio = "hist"
    comp.fX1LegRatio = 0.15
    comp.fX1LegSpectrum = 0.45
    comp.fLegLineHeight = 0.055
    comp.fLogUpperSpace = 2  # this factor will be used to adjust the y axis in log scale
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
    parser.add_argument('yaml',
                        help='YAML configuration file')
    parser.add_argument('--inclusive-train', metavar='TRAIN',
                        help='Inclusive jet train (e.g. Jets_EMC_PbPb_1872)')
    parser.add_argument('--meson', metavar='MESON',
                        default="D0")
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    parser.add_argument('--name', metavar='NAME',
                        default=None)
    parser.add_argument("--refl", action='store_const',
                        default=False, const=True,
                        help='Do not use reflections (only for raw spectra).')
    parser.add_argument("--no-fd", action='store_const',
                        default=False, const=True,
                        help='Do not use B feed-down correction (only for raw spectra).')
    parser.add_argument('--raw-yield-method', metavar='METHOD',
                        default="SideBand",
                        help='Raw yield method')
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    c = yaml.load(f)
    f.close()

    main(c, args.meson, args.jet_type, args.jet_radius, args.name, args.refl, args.no_fd, args.raw_yield_method, args.inclusive_train)

    IPython.embed()
