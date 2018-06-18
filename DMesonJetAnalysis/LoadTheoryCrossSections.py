#!/usr/bin/env python

import re
import ROOT
import DMesonJetUtils

def GetD0JetTheoryCrossSectionAll(config, axis, no_pt_cut=False):
    for t in config["theory"]:
        if not t["active"]: continue
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
        h = GetD0JetTheoryCrossSection(config["input_path"], t["gen"], t["proc"], t["ts"], scale, config["theory_spectrum"], jet_type, axis, normalize)
        t["histogram"] = h
        if no_pt_cut:
            no_pt_cut_spectrum_name = re.sub("DPt_.", "DPt_0", config["theory_spectrum"])
            h = GetD0JetTheoryCrossSection(config["input_path"], t["gen"], t["proc"], t["ts"], scale, no_pt_cut_spectrum_name, jet_type, axis, normalize)
            t["histogram_no_pt_cut"] = h

def GetD0JetTheoryCrossSection(input_path, gen, proc, ts, scale, spectrum, jet_type, axis, normalize):
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
        h.Scale(1.0 / h.Integral(1, h.GetNbinsX()), "width")
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
        if t["inclusive"]["ts"] in inclusive_jet_cross_sections:
            t["inclusive"] = inclusive_jet_cross_sections[t["inclusive"]["ts"]]
            continue
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
        h = GetInclusiveJetTheoryCrossSection(config["input_path"], t["inclusive"]["gen"], t["inclusive"]["proc"], t["inclusive"]["ts"], scale, t["inclusive"]["file_name"], jet_type)
        t["inclusive"]["histogram"] = h
        inclusive_jet_cross_sections[t["inclusive"]["ts"]] = t["inclusive"]
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
