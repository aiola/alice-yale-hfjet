#!/usr/bin/env python
# python script to compare inclusive jet spectra
# us it with
# InclusiveJetCompare.yaml
# InclusiveJetCompareKTCutSuppFact20.yaml
# InclusiveJetCompareKTCutSuppFact60.yaml
# InclusiveJetCompareSuppFactKT1.yaml
# InclusiveJetCompareSuppFactKT10.yaml

import argparse
import numpy
import yaml
import IPython
import ROOT
import DMesonJetCompare
import DMesonJetUtils
import LoadInclusiveJetSpectrum

globalList = []

def LoadHistograms(config, ts, file_name, prefix, title, jet_type):
    result = dict()
    full_file_name = "{}/{}_{}/{}".format(config["input_path"], prefix, ts, file_name)
    myfile = ROOT.TFile(full_file_name, "read")
    for element in config["histograms"]:
        hname = "{}/{}".format(jet_type, element["name"])
        h = DMesonJetUtils.GetObject(myfile, hname)
        if not h: 
            continue
        if "bins" in element:
            h_orig = h
            h = DMesonJetUtils.Rebin1D_fromBins(h_orig, "{}_rebinned".format(h_orig.GetName()), len(element["bins"])-1, numpy.array(element["bins"], dtype=numpy.float64))
        if "max" in element and "min" in element:
            h.GetXaxis().SetRangeUser(element["min"], element["max"])
        h.SetTitle(title)
        result[hname] = h.Clone(hname)
        globalList.append(result[hname])
    return result

def main(config, measured):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if measured: measured_inclusive_cross_section, _ = LoadInclusiveJetSpectrum.GetCrossSection("original")

    histograms = []
    for element in config["list"]:
        if not element["active"]: 
            continue
        if "file_name" in element:
            file_name = element["file_name"]
        else:
            file_name = config["file_name"]
        if "prefix" in element:
            prefix = element["prefix"]
        else:
            prefix = config["prefix"]
        if "jet_type" in element:
            jet_type = element["jet_type"]
        else:
            jet_type = config["jet_type"]
        histograms.append(LoadHistograms(config, element["ts"], file_name, prefix, element["title"], jet_type))
    for element in config["histograms"]:
        if "jet_type" in element:
            jet_type = element["jet_type"]
        else:
            jet_type = config["jet_type"]
        hname = "{}/{}".format(jet_type, element["name"])
        histo_to_compare = [element[hname] for element in histograms if hname in element]
        globalList.extend(histo_to_compare)
        comp = DMesonJetCompare.DMesonJetCompare(hname)
        comp.fX1LegRatio = 0.40
        comp.fX2LegRatio = 0.90
        comp.fX1LegSpectrum = 0.40
        comp.fX2LegSpectrum = 0.90
        comp.fLogUpperSpace = 50
        comp.fLinUpperSpace = 0.4

        if measured and "JetPtExtended" in hname:
            if not histo_to_compare:
                print("No histograms to compare!")
                continue
            globalList.append(measured_inclusive_cross_section)
            r = comp.CompareSpectra(measured_inclusive_cross_section, histo_to_compare)
        else:
            if len(histo_to_compare) <= 1:
                print("No histograms to compare!")
                continue
            r = comp.CompareSpectra(histo_to_compare[0], histo_to_compare[1:])
        for obj in r:
            globalList.append(obj)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare inclusive jet spectra.')
    parser.add_argument('config',
                        help='YAML file')
    parser.add_argument('--measured', action='store_const',
                        default=False, const=True)
    args = parser.parse_args()

    f = open(args.config, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig, args.measured)

    IPython.embed()