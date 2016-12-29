#!/usr/bin/env python
# python script to do extract B feed down correction factors

import argparse
import yaml
import IPython
import ROOT
import DMesonJetCompare
import os

globalList = []

def main(config, meson_name, jet_type, jet_radius, spectrum):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    comp = DMesonJetCompare.DMesonJetCompare("AverageRawYieldVsDefault")
    comp.fOptRatio = "hist"
    comp.fColors = [ROOT.kBlue + 2, ROOT.kGreen + 2]
    comp.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle]
    comp.fFills = [3004] * 2
    default_vs_average_raw_yield(comp, "InvMassFit", config, meson_name, jet_type, jet_radius, spectrum)
    comp.fColors = [ROOT.kRed + 2, ROOT.kOrange + 2]
    comp.fMarkers = [ROOT.kOpenSquare, ROOT.kFullSquare]
    comp.fFills = [3005] * 2
    default_vs_average_raw_yield(comp, "SideBand", config, meson_name, jet_type, jet_radius, spectrum)
    outputPath = "{0}/{1}/{2}/RawYieldUnc/pdf".format(config["input_path"], config["train"], config["name"])
    os.makedirs(outputPath)
    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            fname = "{0}/{1}.pdf".format(outputPath, obj.GetName())
            obj.SaveAs(fname)

def default_vs_average_raw_yield(comp, method, config, meson_name, jet_type, jet_radius, spectrum):
    default_spectrum = GetDefaaultSpectrum(config, meson_name, jet_type, jet_radius, "_".join([spectrum, method]))
    default_spectrum.SetTitle("Def Fit, {0}".format(method))
    average_spectrum = GetAverageSpectrum(method, config)
    average_spectrum.SetTitle("Avg Raw Yield Extr Trials, {0}".format(method))
    average_spectrum.GetXaxis().SetTitle(default_spectrum.GetXaxis().GetTitle())
    average_spectrum.GetYaxis().SetTitle(default_spectrum.GetYaxis().GetTitle())
    globalList.append(default_spectrum)
    globalList.append(average_spectrum)
    comp.SetRatioRelativeUncertaintyFromHistogram(average_spectrum)
    comp.fRatioRelativeUncertainty.SetTitle("Rel. Syst. Unc., {0}".format(method))
    r = comp.CompareSpectra(average_spectrum, [default_spectrum])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

def GetAverageSpectrum(method, config):
    fname = "{0}/{1}/{2}/RawYieldUnc/FinalRawYieldCentralPlusSystUncertainty_Dzero_{3}.root".format(config["input_path"], config["train"], config["name"], method)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        file.ls()
        exit(1)
    h = file.Get("JetRawYieldCentral")
    if not h:
        print("Could not find histogram {0} in file {1}".format("JetRawYieldCentral", fname))
        file.ls()
        exit(1)
    h_copy = h.Clone("{0}_copy".format(h.GetName()))
    return h_copy

def GetDefaaultSpectrum(config, meson_name, jet_type, jet_radius, spectrum):
    fname = "{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"])
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        file.ls()
        exit(1)
    mesonlistname = meson_name
    mesonlist = file.Get(mesonlistname)
    if not mesonlist:
        print("Could not get list {0} from file {1}".format(mesonlistname, file.GetName()))
        file.ls()
        exit(1)
    jetlistname = "{0}_{1}".format(jet_type, jet_radius)
    jetlist = mesonlist.FindObject(jetlistname)
    if not jetlist:
        print("Could not get list {0} from list {1} in file {2}".format(jetlistname, mesonlistname, file.GetName()))
        mesonlist.Print()
        exit(1)
    spectrumlistname = "_".join([meson_name, jet_type, jet_radius, spectrum])
    spectrumlist = jetlist.FindObject(spectrumlistname)
    if not spectrumlist:
        print("Could not get list {0} from list {1} in list {2} in file {3}".format(spectrumlistname, jetlistname, mesonlistname, file.GetName()))
        jetlist.Print()
        exit(1)
    spectrumname = "_".join([meson_name, jet_type, jet_radius, spectrum])
    h = spectrumlist.FindObject(spectrumname)
    if not h:
        print("Could not find object {0} in list {1}/{2}/{3} in file {4}".format(spectrumname, mesonlistname, jetlistname, spectrumlistname, file.GetName()))
        jetlist.Print()
        exit(1)
    h_copy = h.Clone("{0}_copy".format(spectrum))
    return h_copy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='B feed-down.')
    parser.add_argument('yaml', metavar='config.yaml')
    parser.add_argument('--meson', metavar='MESON',
                        default="D0")
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    parser.add_argument('--spectrum', metavar='spectrum',
                        default="JetPtSpectrum_DPt_30")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.meson, args.jet_type, args.jet_radius, args.spectrum)

    IPython.embed()
