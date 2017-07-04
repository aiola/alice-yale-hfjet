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

ptBins = numpy.array([0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 60, 80, 100], dtype=numpy.float32)
rhoBins = numpy.array(list(DMesonJetUtils.frange(0, 15, 0.5, True)), dtype=numpy.float32)
centBins = numpy.array(list(DMesonJetUtils.frange(0, 100, 5, True)), dtype=numpy.float32)

def PlotRhoVsCent(Files, taskName, label=""):
    hname = "{}/fHist{}RhoVsCent".format(taskName, label)
    cname = hname.replace("AliAnalysisTask", "").replace("_histos", "").replace("fHist", "").replace('/', "_")
    pname = "{}Profile".format(cname)
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas(cname, cname)
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, cname, len(centBins) - 1, centBins, len(rhoBins) - 1, rhoBins)
    prof = h_rebin.ProfileX(pname, 1, -1, "i")
    prof.GetYaxis().SetTitle("<#rho> (GeV/#it{c})")
    globalList.append(prof)

    prof_stddev = h_rebin.ProfileX("{}_profSigma".format(pname), 1, -1, "s")
    std_dev = ROOT.TH1F("{}_Sigma".format(pname), "{}_Sigma".format(pname), h_rebin.GetNbinsX(), h_rebin.GetXaxis().GetXbins().GetArray())
    std_dev.GetYaxis().SetTitle("#sigma(#rho) (GeV/#it{c})")
    for i in range(0, prof_stddev.GetNbinsX() + 2):
        std_dev.SetBinContent(i, prof_stddev.GetBinError(i))
    globalList.append(std_dev)
    return prof, std_dev

def PlotRhoVsLeadJetPt(Files, taskName, jetname, label=""):
    hname = "{}/{}_fHist{}RhoVsLeadJetPt".format(taskName, jetname, label)
    cname = hname.replace("AliAnalysisTask", "").replace("_histos", "").replace("fHist", "").replace('/', "_")
    pname = "{}Profile".format(cname)
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas(cname, cname)
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, cname, len(ptBins) - 1, ptBins, len(rhoBins) - 1, rhoBins)
    prof = h_rebin.ProfileX(pname, 1, -1, "i")
    prof.GetYaxis().SetTitle("<#rho> (GeV/#it{c})")
    globalList.append(prof)

    prof_stddev = h_rebin.ProfileX("{}_profSigma".format(pname), 1, -1, "s")
    std_dev = ROOT.TH1F("{}_Sigma".format(pname), "{}_Sigma".format(pname), h_rebin.GetNbinsX(), h_rebin.GetXaxis().GetXbins().GetArray())
    std_dev.GetYaxis().SetTitle("#sigma(#rho) (GeV/#it{c})")
    for i in range(0, prof_stddev.GetNbinsX() + 2):
        std_dev.SetBinContent(i, prof_stddev.GetBinError(i))
    globalList.append(std_dev)
    return prof, std_dev

def PlotRhoVsTrackPt(Files, taskName, label=""):
    hname = "{}/fHist{}RhoVsLeadTrackPt".format(taskName, label)
    cname = hname.replace("AliAnalysisTask", "").replace("_histos", "").replace("fHist", "").replace('/', "_")
    pname = "{}Profile".format(cname)
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas(cname, cname)
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, cname, len(ptBins) - 1, ptBins, len(rhoBins) - 1, rhoBins)
    prof = h_rebin.ProfileX(pname, 1, -1, "i")
    prof.GetYaxis().SetTitle("<#rho> (GeV/#it{c})")
    globalList.append(prof)

    prof_stddev = h_rebin.ProfileX("{}_profSigma".format(pname), 1, -1, "s")
    std_dev = ROOT.TH1F("{}_Sigma".format(pname), "{}_Sigma".format(pname), h_rebin.GetNbinsX(), h_rebin.GetXaxis().GetXbins().GetArray())
    std_dev.GetYaxis().SetTitle("#sigma(#rho) (GeV/#it{c})")
    for i in range(0, prof_stddev.GetNbinsX() + 2):
        std_dev.SetBinContent(i, prof_stddev.GetBinError(i))
    globalList.append(std_dev)
    return prof, std_dev

def PlotLeadJetPtVsCent(Files, taskName, jetname):
    hname = "{}/{}_fHistLeadJetPtVsCent".format(taskName, jetname)
    cname = hname.replace("AliAnalysisTask", "").replace("_histos", "").replace("fHist", "").replace('/', "_")
    pname = "{}Profile".format(cname)
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas(cname, cname)
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, cname, len(centBins) - 1, centBins, len(ptBins) - 1, ptBins)
    prof = h_rebin.ProfileX(pname, 1, -1, "i")
    prof.GetYaxis().SetTitle("<#it{p}_{T,jet}^{lead}> (GeV/#it{c})")
    globalList.append(prof)

    prof_stddev = h_rebin.ProfileX("{}_profSigma".format(pname), 1, -1, "s")
    std_dev = ROOT.TH1F("{}_Sigma".format(pname), "{}_Sigma".format(pname), h_rebin.GetNbinsX(), h_rebin.GetXaxis().GetXbins().GetArray())
    std_dev.GetYaxis().SetTitle("#sigma(#it{p}_{T,jet}^{lead}) (GeV/#it{c})")
    for i in range(0, prof_stddev.GetNbinsX() + 2):
        std_dev.SetBinContent(i, prof_stddev.GetBinError(i))
    globalList.append(std_dev)
    return prof, std_dev

def PlotLeadTrackPtVsCent(Files, taskName):
    hname = "{}/fHistLeadTrackPtVsCent".format(taskName)
    cname = hname.replace("AliAnalysisTask", "").replace("_histos", "").replace("fHist", "").replace('/', "_")
    pname = "{}Profile".format(cname)
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas(cname, cname)
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, cname, len(centBins) - 1, centBins, len(ptBins) - 1, ptBins)
    prof = h_rebin.ProfileX(pname, 1, -1, "i")
    prof.GetYaxis().SetTitle("<#it{p}_{T,track}^{lead}> (GeV/#it{c})")
    globalList.append(prof)

    prof_stddev = h_rebin.ProfileX("{}_profSigma".format(pname), 1, -1, "s")
    std_dev = ROOT.TH1F("{}_Sigma".format(pname), "{}_Sigma".format(pname), h_rebin.GetNbinsX(), h_rebin.GetXaxis().GetXbins().GetArray())
    std_dev.GetYaxis().SetTitle("#sigma(#it{p}_{T,track}^{lead}) (GeV/#it{c})")
    for i in range(0, prof_stddev.GetNbinsX() + 2):
        std_dev.SetBinContent(i, prof_stddev.GetBinError(i))
    globalList.append(std_dev)
    return prof, std_dev

def PlotRhoVsDPt(Files, taskName, meson_name, jet_type, jet_radius):
    hname = "{taskName}_histos/histos{taskName}/D0/Jet_AKT{jet_type}{jet_radius}_pt_scheme/fHistRhoVsLeadDPt".format(taskName=taskName, meson_name=meson_name, jet_type=jet_type, jet_radius=jet_radius)
    cname = "{}_RhoVsDPt".format(meson_name)
    pname = "{}Profile".format(cname)
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas(cname, cname)
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, cname, len(ptBins) - 1, ptBins, len(rhoBins) - 1, rhoBins)
    prof = h_rebin.ProfileX(pname, 1, -1, "i")
    prof.GetYaxis().SetTitle("<#rho> (GeV/#it{c})")
    globalList.append(prof)

    prof_stddev = h_rebin.ProfileX("{}_profSigma".format(pname), 1, -1, "s")
    std_dev = ROOT.TH1F("{}_Sigma".format(pname), "{}_Sigma".format(pname), h_rebin.GetNbinsX(), h_rebin.GetXaxis().GetXbins().GetArray())
    std_dev.GetYaxis().SetTitle("#sigma(#rho) (GeV/#it{c})")
    for i in range(0, prof_stddev.GetNbinsX() + 2):
        std_dev.SetBinContent(i, prof_stddev.GetBinError(i))
    globalList.append(std_dev)
    return prof, std_dev

def PlotRhoVsJetPt(Files, taskName, meson_name, jet_type, jet_radius):
    hname = "{taskName}_histos/histos{taskName}/D0/Jet_AKT{jet_type}{jet_radius}_pt_scheme/fHistRhoVsLeadJetPt".format(taskName=taskName, meson_name=meson_name, jet_type=jet_type, jet_radius=jet_radius)
    cname = "{}_RhoVsJetPt".format(meson_name)
    pname = "{}Profile".format(cname)
    h = DMesonJetUtils.GetObjectAndMerge(Files, hname)
    globalList.append(h)
    c = ROOT.TCanvas(cname, cname)
    globalList.append(c)
    c.cd()
    c.SetLogz()
    h_copy = h.DrawCopy("colz")
    h_copy.GetYaxis().SetRangeUser(0, 50)
    globalList.append(h_copy)

    h_rebin = DMesonJetUtils.Rebin2D_fromBins(h, cname, len(ptBins) - 1, ptBins, len(rhoBins) - 1, rhoBins)
    prof = h_rebin.ProfileX(pname, 1, -1, "i")
    prof.GetYaxis().SetTitle("<#rho> (GeV/#it{c})")
    globalList.append(prof)

    prof_stddev = h_rebin.ProfileX("{}_profSigma".format(pname), 1, -1, "s")
    std_dev = ROOT.TH1F("{}_Sigma".format(pname), "{}_Sigma".format(pname), h_rebin.GetNbinsX(), h_rebin.GetXaxis().GetXbins().GetArray())
    std_dev.GetYaxis().SetTitle("#sigma(#rho) (GeV/#it{c})")
    for i in range(0, prof_stddev.GetNbinsX() + 2):
        std_dev.SetBinContent(i, prof_stddev.GetBinError(i))
    globalList.append(std_dev)
    return prof, std_dev

def main(config, meson_name, jet_type, jet_radius, gen, cent):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    FileNames = []

    path = "{0}/{1}".format(config["input_path"], config["train"])
    print("Looking for file {0} in path {1}".format(config["file_name"], path))
    FileNames.extend(DMesonJetUtils.find_file(path, config["file_name"]))

    trigger = config["analysis"][0]["trigger"][0]
    if trigger:
        taskName = "{}_{}".format(config["task_name"], trigger)
    else:
        taskName = config["task_name"]

    Files = []
    for fname in FileNames:
        f = ROOT.TFile(fname, "read")
        Files.append(f)

    drawSingle = []
    rhoVScent = []
    rhoVStrackPt = []
    rhoVSjetPt = []
    rhoVScentGen = []
    rhoVStrackPtGen = []
    rhoVSjetPtGen = []

    allCompare = {"rhoVScent" : rhoVScent,
                  "rhoVStrackPt" : rhoVStrackPt,
                  "rhoVSjetPt" : rhoVSjetPt,
                  "rhoVScentGen" : rhoVScentGen,
                  "rhoVStrackPtGen" : rhoVStrackPtGen,
                  "rhoVSjetPtGen" : rhoVSjetPtGen }

    if cent:
        drawSingle.extend(PlotLeadTrackPtVsCent(Files, "AliAnalysisTaskRhoDev_Rho_histos"))
        drawSingle.extend(PlotLeadJetPtVsCent(Files, "AliAnalysisTaskRhoDev_Rho_histos", "Background"))
        rhoVScent.append(PlotRhoVsCent(Files, "AliAnalysisTaskRhoDev_Rho_histos"))
        rhoVScent.append(PlotRhoVsCent(Files, "AliAnalysisTaskRhoTransDev_RhoTrans_histos"))
        rhoVScent.append(PlotRhoVsCent(Files, "AliAnalysisTaskRhoTransDev_RhoTrans_histos", "B2B"))
        if gen:
            drawSingle.extend(PlotLeadTrackPtVsCent(Files, "AliAnalysisTaskRhoDev_RhoGen_histos"))
            drawSingle.extend(PlotLeadJetPtVsCent(Files, "AliAnalysisTaskRhoDev_RhoGen_histos", "Background"))
            rhoVScentGen.append(PlotRhoVsCent(Files, "AliAnalysisTaskRhoDev_RhoGen_histos"))
            rhoVScentGen.append(PlotRhoVsCent(Files, "AliAnalysisTaskRhoTransDev_RhoTransGen_histos"))
            rhoVScentGen.append(PlotRhoVsCent(Files, "AliAnalysisTaskRhoTransDev_RhoTransGen_histos", "B2B"))

    rhoVStrackPt.append(PlotRhoVsTrackPt(Files, "AliAnalysisTaskRhoDev_Rho_histos"))
    rhoVStrackPt.append(PlotRhoVsTrackPt(Files, "AliAnalysisTaskRhoTransDev_RhoTrans_histos"))
    rhoVStrackPt.append(PlotRhoVsTrackPt(Files, "AliAnalysisTaskRhoTransDev_RhoTrans_histos", "B2B"))

    rhoVSjetPt.append(PlotRhoVsLeadJetPt(Files, "AliAnalysisTaskRhoDev_Rho_histos", "Background"))
    rhoVSjetPt.append(PlotRhoVsLeadJetPt(Files, "AliAnalysisTaskRhoTransDev_RhoTrans_histos", "Signal"))
    rhoVSjetPt.append(PlotRhoVsLeadJetPt(Files, "AliAnalysisTaskRhoTransDev_RhoTrans_histos", "Signal", "B2B"))

    drawSingle.extend(PlotRhoVsDPt(Files, taskName, meson_name, jet_type, jet_radius))
    drawSingle.extend(PlotRhoVsJetPt(Files, taskName, meson_name, jet_type, jet_radius))

    if gen:
        rhoVStrackPtGen.append(PlotRhoVsTrackPt(Files, "AliAnalysisTaskRhoDev_RhoGen_histos"))
        rhoVStrackPtGen.append(PlotRhoVsTrackPt(Files, "AliAnalysisTaskRhoTransDev_RhoTransGen_histos"))
        rhoVStrackPtGen.append(PlotRhoVsTrackPt(Files, "AliAnalysisTaskRhoTransDev_RhoTransGen_histos", "B2B"))
        rhoVSjetPtGen.append(PlotRhoVsLeadJetPt(Files, "AliAnalysisTaskRhoDev_RhoGen_histos", "Background"))
        rhoVSjetPtGen.append(PlotRhoVsLeadJetPt(Files, "AliAnalysisTaskRhoTransDev_RhoTransGen_histos", "Signal"))
        rhoVSjetPtGen.append(PlotRhoVsLeadJetPt(Files, "AliAnalysisTaskRhoTransDev_RhoTransGen_histos", "Signal", "B2B"))

        drawSingle.extend(PlotRhoVsDPt(Files, taskName, "{}_MCTruth".format(meson_name), jet_type, jet_radius))
        drawSingle.extend(PlotRhoVsJetPt(Files, taskName, "{}_MCTruth".format(meson_name), jet_type, jet_radius))

    for h in drawSingle:
        cname = h.GetName()
        c = ROOT.TCanvas(cname, cname)
        globalList.append(c)
        c.cd()
        if isinstance(h, ROOT.TH2):
            h.Draw("colz")
        else:
            h.Draw("")

    for name, pair_list in allCompare.iteritems():
        def PlotCompare(name, pair_list, i):
            hlist = [pair[i] for pair in pair_list]
            comp = DMesonJetCompare.DMesonJetCompare(name)
            comp.fDoSpectraPlot = "lineary"
            comp.fDoRatioPlot = False
            if i == 1:
                comp.fOptSpectrumBaseline = "hist"
                comp.fOptSpectrum = "hist"
            if len(hlist) >= 2:
                r = comp.CompareSpectra(hlist[0], hlist[1:])
                globalList.extend(r)
            else:
                print("Skipping {}".format(name))

        name_avg = name
        PlotCompare(name_avg, pair_list, 0)

        name_sigma = "{}_Sigma".format(name)
        PlotCompare(name_sigma, pair_list, 1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Underlying event studies.')
    parser.add_argument('yaml', metavar='config.yaml')
    parser.add_argument('--meson', metavar='MESON',
                        default="D0")
    parser.add_argument('--jet-type', metavar='TYPE',
                        default="Charged")
    parser.add_argument('--jet-radius', metavar='RADIUS',
                        default="R040")
    parser.add_argument('--cent', action='store_const',
                        default=False, const=True,
                        help='Centrality available.')
    parser.add_argument('--gen', action='store_const',
                        default=False, const=True,
                        help='Generator level.')
    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.meson, args.jet_type, args.jet_radius, args.gen, args.cent)

    IPython.embed()
