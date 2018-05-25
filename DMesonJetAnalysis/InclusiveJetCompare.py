#!/usr/bin/env python
# python script to project TTree containing inclusive jet spectra

import argparse
import IPython
import ROOT
import DMesonJetCompare
import yaml
import LoadInclusiveJetSpectrum

globalList = []

def LoadHistograms(config, ts, file_name, title):
    result = dict()
    full_file_name = "{}/{}_{}/{}".format(config["input_path"], config["prefix"], ts, file_name)
    file = ROOT.TFile(full_file_name, "read")
    for element in config["histograms"]:
        hname = element["name"]
        h = file.Get(hname)
        if not h:
            print("Could not find object '{}' in file '{}'. Skipping.".format(hname, full_file_name))
            continue
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
        histograms.append(LoadHistograms(config, element["ts"], file_name, element["title"]))
    for element in config["histograms"]:
        hname = element["name"]
        histo_to_compare = [element[hname] for element in histograms if hname in element]
        globalList.extend(histo_to_compare)
        comp = DMesonJetCompare.DMesonJetCompare(hname)
        
        if "JetPtExtended" in hname:
            globalList.append(measured_inclusive_cross_section)
            r = comp.CompareSpectra(measured_inclusive_cross_section, histo_to_compare)
        else:
            r = comp.CompareSpectra(histo_to_compare[0], histo_to_compare[1:])
        for obj in r:
            globalList.append(obj)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare inclusive jet spectra.')
    parser.add_argument('config',
                        help='YAML file')
    args = parser.parse_args()

    f = open(args.config, 'r')
    config = yaml.load(f)
    f.close()

    main(config)

    IPython.embed()