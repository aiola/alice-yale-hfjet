#!/usr/bin/env python
# python script to compare inclusive jet spectra
# us it with
# InclusiveJetCompare.yaml
# InclusiveJetCompareKTCut.yaml
# InclusiveJetCompareSuppFactKT1.yaml
# InclusiveJetCompareSuppFactKT10.yaml

import argparse
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
    file = ROOT.TFile(full_file_name, "read")
    for element in config["histograms"]:
        hname = "{}/{}".format(jet_type, element["name"])
        h = DMesonJetUtils.GetObject(file, hname)
        if not h: continue
        if "max" in element and "min" in element:
            h.GetXaxis().SetRangeUser(element["min"], element["max"])
        h.SetTitle(title)
        result[hname] = h.Clone(hname)
        globalList.append(result[hname])
    return result

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    measured_inclusive_cross_section, _ = LoadInclusiveJetSpectrum.GetCrossSection("original")

    histograms = []
    for element in config["list"]:
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
        print(histo_to_compare)
        globalList.extend(histo_to_compare)
        comp = DMesonJetCompare.DMesonJetCompare(hname)

        if "JetPtExtended" in hname:
            if len(histo_to_compare) == 0:
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
    args = parser.parse_args()

    f = open(args.config, 'r')
    yconfig = yaml.load(f)
    f.close()

    main(yconfig)

    IPython.embed()