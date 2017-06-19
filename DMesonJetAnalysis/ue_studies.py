#!/usr/bin/env python
# python script for underlying event studies

import argparse
import yaml
import IPython
import ROOT
import os
import DMesonJetCompare
import DMesonJetUtils
import subprocess

globalList = []

def PlotRhoVsCent(Files):
    hname = "AliAnalysisTaskRho_histos/fHistRhovsCent"
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas("RhoVsCent", "RhoVsCent")
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h.RebinX(2)
    prof = h.ProfileX("RhoVsCentProfile", 1, -1, "s")
    prof.GetYaxis().SetTitle("<#rho> (GeV/#it{c})")
    globalList.append(prof)
    cp = ROOT.TCanvas("RhoVsCentProfile", "RhoVsCentProfile")
    globalList.append(cp)
    cp.cd()
    prof.Draw("e4")
    prof.SetFillStyle(1001)
    prof.SetFillColor(ROOT.kRed)

def PlotRhoVsDPt(Files):
    hname = "AliAnalysisTaskDmesonJets_AnyINT_histos/histosAliAnalysisTaskDmesonJets_AnyINT/D0/Jet_AKTChargedR040_pt_scheme/fHistRhoVsLeadDPt"
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas("RhoVsDPt", "RhoVsDPt")
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h.RebinX(2)
    prof = h.ProfileX("RhoVsDPtProfile", 1, -1, "i")
    prof.GetYaxis().SetTitle("<#rho> (GeV/#it{c})")
    globalList.append(prof)
    cp = ROOT.TCanvas("RhoVsDPtProfile", "RhoVsDPtProfile")
    globalList.append(cp)
    cp.cd()
    prof.Draw("")

def PlotRhoVsJetPt(Files):
    hname = "AliAnalysisTaskDmesonJets_AnyINT_histos/histosAliAnalysisTaskDmesonJets_AnyINT/D0/Jet_AKTChargedR040_pt_scheme/fHistRhoVsLeadJetPt"
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas("RhoVsJetPt", "RhoVsJetPt")
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h.RebinX(2)
    prof = h.ProfileX("RhoVsJetPtProfile", 1, -1, "i")
    prof.GetYaxis().SetTitle("<#rho> (GeV/#it{c})")
    globalList.append(prof)
    cp = ROOT.TCanvas("RhoVsJetPtProfile", "RhoVsJetPtProfile")
    globalList.append(cp)
    cp.cd()
    prof.Draw("")

def main(config, meson_name, jet_type, jet_radius):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    FileNames = []

    path = "{0}/{1}".format(config["input_path"], config["train"])
    print("Looking for file {0} in path {1}".format(config["file_name"], path))
    FileNames.extend(DMesonJetUtils.find_file(path, config["file_name"]))

    Files = []
    for fname in FileNames:
        f = ROOT.TFile(fname, "read")
        Files.append(f)

    PlotRhoVsCent(Files)
    PlotRhoVsDPt(Files)
    PlotRhoVsJetPt(Files)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Underlying event studies.')
    parser.add_argument('yaml', metavar='config.yaml')
    parser.add_argument('--meson', metavar='MESON',
                        default="D0")
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.meson, args.jet_type, args.jet_radius)

    IPython.embed()
