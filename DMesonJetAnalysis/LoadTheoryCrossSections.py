#!/usr/bin/env python

import ROOT
import DMesonJetUtils

def GetD0JetTheoryCrossSectionAll(config, axis):
    for t in config["theory"]:
        if not t["active"]: continue
        if "spectrum_name" in t:
            spectrum_name = t["spectrum_name"]
        else:
            spectrum_name = config["theory_spectrum"]
        if "jet_type" in t:
            jet_type = t["jet_type"]
        else:
            jet_type = config["jet_type"]
        if "normalize" in config:
            normalize = config["normalize"]
        else:
            normalize = False
        if "scale" in t and not normalize:
            scale = t["scale"]
        else:
            scale = None

        if normalize:
            if "data_minx" in config and "data_maxx" in config:
                data_minx = config["data_minx"]
                data_maxx = config["data_maxx"]
            else:
                print("Error: requsted normalization, but did not provide minx and maxx via the fields 'data_minx' and 'data_maxx'. Aborting.")
                exit(1)
        else:
            data_minx = None
            data_maxx = None

        h = GetD0JetTheoryCrossSection(config["input_path"], t["gen"], t["proc"], t["ts"], scale, spectrum_name, jet_type, axis, normalize, data_minx, data_maxx)
        t["histogram"] = h

def GetD0JetTheoryCrossSection(input_path, gen, proc, ts, scale, spectrum, jet_type, axis, normalize, data_minx, data_maxx):
    fname = "{input_path}/FastSim_{gen}_{proc}_{ts}/FastSimAnalysis_ccbar_{gen}_{proc}_{ts}.root".format(input_path=input_path, gen=gen, proc=proc, ts=ts)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    h_orig = DMesonJetUtils.GetObject(file, "D0_MCTruth/{jet_type}/D0_MCTruth_{jet_type}_{spectrum}/D0_MCTruth_{jet_type}_{spectrum}".format(spectrum=spectrum, jet_type=jet_type))
    if not h_orig:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)

    h = DMesonJetUtils.Rebin1D(h_orig, axis)
    if normalize:
        h.Scale(1.0 / h.Integral(h.GetXaxis().FindBin(data_minx*1.0001), h.GetXaxis().FindBin(data_maxx*0.9999)), "width")
    else:
        if scale:
            h.Scale(scale * 0.5, "width")  # particle/antiparticle
        else:
            h.Scale(0.5, "width")  # particle/antiparticle

    return h

def GetInclusiveJetTheoryCrossSectionAll(config):
    inclusive_jet_cross_sections = dict()
    for t in config["theory"]:
        if not t["active"]: continue
        if not t["inclusive"]: continue
        if "jet_type" in t["inclusive"]:
            jet_type = t["inclusive"]["jet_type"]
        elif "jet_type" in t:
            jet_type = t["jet_type"]
        else:
            jet_type = config["jet_type"]
        if "scale" in t["inclusive"]:
            scale = t["inclusive"]["scale"]
        else:
            scale = None
        if not "color" in t["inclusive"] and "color" in t:
            t["inclusive"]["color"] = t["color"]
        if not "title" in t["inclusive"] and "title" in t:
            t["inclusive"]["title"] = t["title"]
        if not "line" in t["inclusive"] and "line" in t:
            t["inclusive"]["line"] = t["line"]

        jet_key = (t["inclusive"]["ts"], jet_type)

        if jet_key in inclusive_jet_cross_sections:
            t["inclusive"] = inclusive_jet_cross_sections[jet_key]
            continue

        h = GetInclusiveJetTheoryCrossSection(config["input_path"], t["inclusive"]["gen"], t["inclusive"]["proc"], t["inclusive"]["ts"], scale, t["inclusive"]["file_name"], jet_type)
        t["inclusive"]["histogram"] = h
        inclusive_jet_cross_sections[jet_key] = t["inclusive"]
    return inclusive_jet_cross_sections

def GetInclusiveJetTheoryCrossSection(input_path, gen, proc, ts, scale, file_name, jet_type):
    fname = "{input_path}/FastSim_{gen}_{proc}_{ts}/{file_name}".format(input_path=input_path, gen=gen, proc=proc, ts=ts, file_name=file_name)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    if jet_type:
        objname = "{jet_type}/JetPt".format(jet_type=jet_type)
    else:
        objname = "JetPt"
    h_orig = DMesonJetUtils.GetObject(file, objname)
    if not h_orig:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)

    h = h_orig.Clone()

    if scale:
        h.Scale(scale)

    return h
