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
import numpy

globalList = []

ptBins = numpy.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150], dtype=numpy.float32)
rhoBins = numpy.array(list(DMesonJetUtils.frange(0, 15, 0.5, True)), dtype=numpy.float32)
centBins = numpy.array(list(DMesonJetUtils.frange(0, 100, 10, True)), dtype=numpy.float32)

def PlotRhoVsCent(Files):
    hname = "AliAnalysisTaskRhoDev_Rho_histos/fHistRhoVsCent"
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas("RhoVsCent", "RhoVsCent")
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, "RhoVsCent", len(centBins) - 1, centBins, len(rhoBins) - 1, rhoBins)
    prof = h_rebin.ProfileX("RhoVsCentProfile", 1, -1, "i")
    prof.GetYaxis().SetTitle("<#rho> (GeV/#it{c})")
    globalList.append(prof)
    cp = ROOT.TCanvas("RhoVsCentProfile", "RhoVsCentProfile")
    globalList.append(cp)
    cp.cd()
    prof.Draw("")

def PlotLeadJetPtVsCent(Files):
    hname = "AliAnalysisTaskRhoDev_Rho_histos/fHistLeadJetPtVsCent"
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas("LeadJetPtVsCent", "LeadJetPtVsCent")
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, "LeadJetPtVsCent", len(centBins) - 1, centBins, len(ptBins) - 1, ptBins)
    prof = h_rebin.ProfileX("LeadJetPtVsCentProfile", 1, -1, "i")
    prof.GetYaxis().SetTitle("#it{p}_{T,jet}^{lead} (GeV/#it{c})")
    globalList.append(prof)
    cp = ROOT.TCanvas("LeadJetPtVsCentProfile", "LeadJetPtVsCentProfile")
    globalList.append(cp)
    cp.cd()
    prof.Draw()

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

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, "RhoVsDPt", len(ptBins) - 1, ptBins, len(rhoBins) - 1, rhoBins)
    prof = h_rebin.ProfileX("RhoVsDPtProfile", 1, -1, "i")
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

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, "RhoVsJetPt", len(ptBins) - 1, ptBins, len(rhoBins) - 1, rhoBins)
    prof = h_rebin.ProfileX("RhoVsJetPtProfile", 1, -1, "i")
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
    PlotLeadJetPtVsCent(Files)

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
