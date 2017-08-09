#!/usr/bin/env python

import ROOT
import math
import array
import IPython
import argparse
import DMesonJetUtils
import DMesonJetCompare
import functools
import numpy

globalList = []

def rsetattr(obj, attr, val):
    pre, _, post = attr.rpartition('.')
    return setattr(rgetattr(obj, pre) if pre else obj, post, val)

sentinel = object()
def rgetattr(obj, attr, default=sentinel):
    if default is sentinel:
        _getattr = getattr
    else:
        def _getattr(obj, name):
            return getattr(obj, name, default)
    return functools.reduce(_getattr, [obj] + attr.split('.'))

globalList = []

def load_efficiency():
    filename = "/Volumes/DATA/ALICE/JetResults/Jets_EMC_pp_MC_1161_1162_1163_1164/LHC15i2_Train1161_response.root"
    file = ROOT.TFile(filename, "read")
    if not file or file.IsZombie():
        print("Could not open {}".format(filename))
        exit(1)
    hname = "D0_Jet_AKTChargedR040_pt_scheme_JetPtSpectrum_DPt_30_Prompt_Extended/D0_Jet_AKTChargedR040_pt_scheme_JetPtSpectrum_DPt_30_Prompt_Extended_Efficiency"
    h = DMesonJetUtils.GetObject(file, hname)
    return h

class histogram_definition:
    def __init__(self, hname, htitle, xbranch, xaxis, xtitle, cuts):
        self.xbranch = xbranch
        self.is_accepted = cuts
        self.limits = []
        self.h = ROOT.TH1F(hname, "{};{};counts".format(htitle, xtitle), len(xaxis) - 1, array.array('d', xaxis))
        self.h.Sumw2()

    def fill(self, obj):
        if not self.is_accepted(obj, self.limits): return
        xval = rgetattr(obj, self.xbranch)
        self.h.Fill(xval)

def ProjectMyTree1D(tree, hdefs):
    for i, obj in enumerate(tree):
        if i % 10000 == 0: print("Event {}".format(i))
        for hdef in hdefs:
            hdef.fill(obj)

def DMesonJetCuts(obj, limits):
    for limit in limits:
        val = rgetattr(obj, limit[0])
        if limit[1] and val < limit[1]: return False
        if limit[2] and val >= limit[2]: return False
    return True

def InclusiveJetCuts(obj, limits):
    return True

def GetDMesonJets(file, dmeson_name):
    treename = "AliAnalysisTaskDmesonJets_{}_MCTruth".format(dmeson_name)
    tree = DMesonJetUtils.GetObject(file, treename)
    jet_pt_axis = [3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
    # start = 0
    # stop = 100
    # step = 1
    # jet_pt_axis = numpy.linspace(start, stop, (stop - start) / step + 1, True)
    charged_jets = histogram_definition("hDMesonChargedJetsPt", "Charged jets", "Jet_AKTChargedR040_pt_scheme.fPt", jet_pt_axis, "#it{p}_{T,ch.jet} (GeV/#it{c})", DMesonJetCuts)
    charged_jets.h.SetMarkerColor(ROOT.kRed + 1)
    charged_jets.h.SetMarkerStyle(ROOT.kFullCircle)
    charged_jets.h.SetMarkerSize(0.9)
    charged_jets.h.SetLineColor(ROOT.kRed + 1)
    charged_jets.limits.append(("DmesonJet.fPt", 3, None))
    charged_jets.limits.append(("Jet_AKTChargedR040_pt_scheme.fEta", -1, 1))
    charged_jets.acceptance = math.pi * 2
    full_jets = histogram_definition("hDMesonFullJetsPt", "Full jets", "Jet_AKTFullR040_pt_scheme.fPt", jet_pt_axis, "#it{p}_{T,jet} (GeV/#it{c})", DMesonJetCuts)
    full_jets.h.SetMarkerColor(ROOT.kBlue + 1)
    full_jets.h.SetMarkerStyle(ROOT.kFullSquare)
    full_jets.h.SetMarkerSize(0.9)
    full_jets.h.SetLineColor(ROOT.kBlue + 1)
    full_jets.limits.append(("DmesonJet.fPt", 3, None))
    full_jets.limits.append(("Jet_AKTFullR040_pt_scheme.fEta", -1, 1))
    full_jets.acceptance = (math.pi * 2 / 360 * 108 - 0.8) * 0.6
    hdefs = [charged_jets, full_jets]
    ProjectMyTree1D(tree, hdefs)
    return hdefs

def do_analysis(filename, dmeson_name):
    histograms = None

    print("Opening file {0}".format(filename))
    file = ROOT.TFile(filename)

    if not file or file.IsZombie():
        print("Could not open file {}".format(filename))
        exit(1)

    w = get_weight(file)

    hdefs = GetDMesonJets(file, dmeson_name)

    for hdef in hdefs:
        hdef.h.Scale(w / 2)  # divide by 2 for the eta acceptance
        hdef.h.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} #Delta#it{p}_{T} (mb)")

    PlotDMesonJets(hdefs)

def PlotDMesonJets(hdefs):
    # Cross sections
    histograms = []
    for hdef in hdefs:
        h = hdef.h.Clone("{}_CrossSection".format(hdef.h.GetName()))
        h.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb / (GeV/#it{c})]")
        h.Scale(1, "width")
        histograms.append(h)
    comp = DMesonJetCompare.DMesonJetCompare("DMesonJetCrossSections")
    comp.fDoRatioPlot = False
    comp.fColors = comp.fColors[1:]
    comp.fMarkers = comp.fMarkers[1:]
    comp.fX1LegSpectrum = 0.45
    comp.fX2LegSpectrum = 0.80
    r = comp.CompareSpectra(histograms[0], histograms[1:])
    globalList.extend(histograms)
    globalList.extend(r)

    # Expected counts
    lumi_list = [15, 50, 150]
    histograms = []
    for lumi in lumi_list:
        if lumi == 15:
            h = hdefs[0].h.Clone("{}_Lumi{}".format(hdefs[0].h.GetName(), lumi))
            h.Scale(hdefs[0].acceptance / (2 * math.pi) * lumi * 1e6)
            h.SetTitle("TPC jets, #it{{L}}_{{int}} = {} nb^{{-1}}".format(lumi))
            h.GetYaxis().SetTitle("counts")
            histograms.append(h)
        h = hdefs[1].h.Clone("{}_Lumi{}".format(hdefs[1].h.GetName(), lumi))
        h.Scale(hdefs[1].acceptance / (2 * math.pi) * lumi * 1e6)
        h.SetTitle("EMCal jets, #it{{L}}_{{int}} = {} nb^{{-1}}".format(lumi))
        h.GetYaxis().SetTitle("counts")
        histograms.append(h)

    comp = DMesonJetCompare.DMesonJetCompare("DMesonJetCounts")
    comp.fX1LegSpectrum = 0.45
    comp.fX2LegSpectrum = 0.80
    comp.fDoRatioPlot = False
    r = comp.CompareSpectra(histograms[0], histograms[1:])
    globalList.extend(histograms)
    globalList.extend(r)

    histograms_BR = []
    for h in histograms:
        h_new = h.Clone("{}_BR".format(h.GetName()))
        h_new.Scale(0.0393)
        histograms_BR.append(h_new)
    comp = DMesonJetCompare.DMesonJetCompare("DMesonJetCounts_BR")
    comp.fDoRatioPlot = False
    comp.fX1LegSpectrum = 0.45
    comp.fX2LegSpectrum = 0.80
    r = comp.CompareSpectra(histograms_BR[0], histograms_BR[1:])
    globalList.extend(histograms_BR)
    globalList.extend(r)

    eff = load_efficiency()

    histograms_BR_eff = []
    for h in histograms_BR:
        h_new = h.Clone("{}_eff".format(h.GetName()))
        h_new.Multiply(eff)
        histograms_BR_eff.append(h_new)
    comp = DMesonJetCompare.DMesonJetCompare("DMesonJetCounts_BR_eff")
    comp.fDoRatioPlot = False
    comp.fX1LegSpectrum = 0.45
    comp.fX2LegSpectrum = 0.80
    r = comp.CompareSpectra(histograms_BR_eff[0], histograms_BR_eff[1:])
    globalList.extend(histograms_BR_eff)
    globalList.extend(r)

def get_weight(file):
    tlist = file.Get("AliAnalysisTaskDmesonJets_histos")
    hTrials = tlist.FindObject("fHistTrialsVsPtHardNoSel")
    hXsec = tlist.FindObject("fHistXsectionVsPtHardNoSel")
    trials = hTrials.Integral()
    xsec = hXsec.GetMean(2)
    w = xsec / trials
    return w

def main(gen, proc, ts, stage, dmeson_name):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    filename = "/Volumes/DATA/ALICE/JetResults/FastSim_{gen}_{proc}_{ts}/stage_{stage}/output/001/AnalysisResults_FastSim_{gen}_{proc}_{ts}.root".format(gen=gen, proc=proc, ts=ts, stage=stage)
    do_analysis(filename, dmeson_name)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trigger studies.')
    parser.add_argument('ts',
                        help='Timestamp of the fast simulation')
    parser.add_argument('--gen',
                        default="powheg+pythia6",
                        help='Generator used in the fast simulation')
    parser.add_argument('--proc',
                        default="charm",
                        help='Process')
    parser.add_argument('--stage',
                        default=1,
                        help='Merging stage')
    parser.add_argument('--d-meson',
                        default="D0",
                        help='D meson')
    args = parser.parse_args()

    main(args.gen, args.proc, args.ts, args.stage, args.d_meson)

    IPython.embed()
