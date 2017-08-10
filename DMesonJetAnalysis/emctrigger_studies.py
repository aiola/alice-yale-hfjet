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
    hname = "D0_Jet_AKTChargedR040_pt_scheme_JetPtDPtSpectrum_Prompt/D0_Jet_AKTChargedR040_pt_scheme_JetPtDPtSpectrum_Prompt_Efficiency_NoJet"
    h = DMesonJetUtils.GetObject(file, hname)
    file.Close()
    return h

class histogram_definition:
    def __init__(self, hname, htitle, xbranch, xaxis, xtitle, cuts):
        self.xbranch = xbranch
        self.is_accepted = cuts
        self.limits = []
        self.h = ROOT.TH1F(hname, "{};{};counts".format(htitle, xtitle), len(xaxis) - 1, array.array('d', xaxis))
        self.h.Sumw2()
        self.h_eff = ROOT.TH1F(hname, "{};{};counts".format(htitle, xtitle), len(xaxis) - 1, array.array('d', xaxis))
        self.h_eff.Sumw2()
        self.eff = None

    def get_w(self, obj):
        if self.eff:
            w = self.eff.GetBinContent(self.eff.GetXaxis().FindBin(obj.DmesonJet.fPt))
        else:
            w = 1
        return w

    def fill(self, obj):
        if not self.is_accepted(obj, self.limits): return
        xval = rgetattr(obj, self.xbranch)
        self.h.Fill(xval)
        if self.eff: self.h_eff.Fill(xval, self.get_w(obj))

def DMesonJetCuts(obj, limits):
    for limit in limits:
        val = rgetattr(obj, limit[0])
        if limit[1] and val < limit[1]: return False
        if limit[2] and val >= limit[2]: return False
    return True

def InclusiveJetCuts(obj, limits):
    return True

def GetD0MesonJets(file):
    return GetDMesonJets(file, "D0")

def GetDStarMesonJets(file):
    return GetDMesonJets(file, "DStar")

def GetDMesonJets(file, dmeson_name):
    treename = "AliAnalysisTaskDmesonJets_{}_MCTruth".format(dmeson_name)
    tree = DMesonJetUtils.GetObject(file, treename)
    jet_pt_axis = [3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
    # start = 0
    # stop = 100
    # step = 1
    # jet_pt_axis = numpy.linspace(start, stop, (stop - start) / step + 1, True)
    eff = load_efficiency()

    w = get_weight(file, "AliAnalysisTaskDmesonJets_histos")

    charged_jets = histogram_definition("hDMesonChargedJetsPt", "Charged jets", "Jet_AKTChargedR040_pt_scheme.fPt", jet_pt_axis, "#it{p}_{T,jet} (GeV/#it{c})", DMesonJetCuts)
    charged_jets.h.SetMarkerColor(ROOT.kRed + 1)
    charged_jets.h.SetMarkerStyle(ROOT.kFullCircle)
    charged_jets.h.SetMarkerSize(0.9)
    charged_jets.h.SetLineColor(ROOT.kRed + 1)
    charged_jets.limits.append(("DmesonJet.fPt", 3, None))
    charged_jets.limits.append(("Jet_AKTChargedR040_pt_scheme.fEta", -1, 1))
    charged_jets.acceptance = math.pi * 2
    charged_jets.eff = eff
    charged_jets.weight = w

    full_jets = histogram_definition("hDMesonFullJetsPt", "Full jets", "Jet_AKTFullR040_pt_scheme.fPt", jet_pt_axis, "#it{p}_{T,jet} (GeV/#it{c})", DMesonJetCuts)
    full_jets.h.SetMarkerColor(ROOT.kBlue + 1)
    full_jets.h.SetMarkerStyle(ROOT.kFullSquare)
    full_jets.h.SetMarkerSize(0.9)
    full_jets.h.SetLineColor(ROOT.kBlue + 1)
    full_jets.limits.append(("DmesonJet.fPt", 3, None))
    full_jets.limits.append(("Jet_AKTFullR040_pt_scheme.fEta", -1, 1))
    full_jets.acceptance = (math.pi * 2 / 360 * 108 - 0.8) * 0.6
    full_jets.eff = eff
    full_jets.weight = w

    hdefs = [charged_jets, full_jets]

    for i, obj in enumerate(tree):
        if i % 10000 == 0: print("Event {}".format(i))
        # if i > 10000: break
        for hdef in hdefs:
            hdef.fill(obj)

    return hdefs

def GetInclusiveJets(file):
    treename = "AliAnalysisTaskEmcalJetTree_jets"
    tree = DMesonJetUtils.GetObject(file, treename)
    jet_pt_axis = [5, 6, 8, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 120, 130, 140, 150, 160, 180, 190, 200, 210, 220, 230, 240, 250]
    # start = 0
    # stop = 100
    # step = 1
    # jet_pt_axis = numpy.linspace(start, stop, (stop - start) / step + 1, True)

    w = get_weight(file, "AliAnalysisTaskEmcalJetTree_histos")

    charged_jets = histogram_definition("hDMesonChargedJetsPt", "Charged jets", "fPt", jet_pt_axis, "#it{p}_{T,jet} (GeV/#it{c})", DMesonJetCuts)
    charged_jets.h.SetMarkerColor(ROOT.kRed + 1)
    charged_jets.h.SetMarkerStyle(ROOT.kFullCircle)
    charged_jets.h.SetMarkerSize(0.9)
    charged_jets.h.SetLineColor(ROOT.kRed + 1)
    charged_jets.limits.append(("fEta", -1, 1))
    charged_jets.acceptance = math.pi * 2
    charged_jets.weight = w
    charged_jets.vector_branch = "Jet_AKTChargedR040_mcparticles_pT0000_pt_scheme"
    full_jets = histogram_definition("hDMesonFullJetsPt", "Full jets", "fPt", jet_pt_axis, "#it{p}_{T,jet} (GeV/#it{c})", DMesonJetCuts)
    full_jets.h.SetMarkerColor(ROOT.kBlue + 1)
    full_jets.h.SetMarkerStyle(ROOT.kFullSquare)
    full_jets.h.SetMarkerSize(0.9)
    full_jets.h.SetLineColor(ROOT.kBlue + 1)
    full_jets.limits.append(("fEta", -1, 1))
    full_jets.vector_branch = "Jet_AKTFullR040_mcparticles_pT0000_pt_scheme"
    full_jets.acceptance = (math.pi * 2 / 360 * 108 - 0.8) * 0.6
    full_jets.weight = w
    hdefs = [charged_jets, full_jets]

    for i, obj in enumerate(tree):
        if i % 10000 == 0: print("Event {}".format(i))
        # if i > 10000: break
        for hdef in hdefs:
            v = getattr(obj, hdef.vector_branch)
            for jet in v:
                hdef.fill(jet)

    return hdefs

def do_analysis(filename, GetHistograms, PlotHistograms):
    histograms = None

    print("Opening file {0}".format(filename))
    file = ROOT.TFile(filename)

    if not file or file.IsZombie():
        print("Could not open file {}".format(filename))
        exit(1)

    hdefs = GetHistograms(file)

    for hdef in hdefs:
        hdef.h.Scale(hdef.weight / 2)  # divide by 2 for the eta acceptance
        hdef.h.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} #Delta#it{p}_{T} (mb)")

        if hdef.h_eff:
            hdef.h_eff.Scale(hdef.weight / 2)  # divide by 2 for the eta acceptance
            hdef.h_eff.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} #Delta#it{p}_{T} #epsilon (mb)")

    outputPath = filename[0:filename.rfind("/")]
    PlotHistograms(hdefs)

    for c in globalList:
        if isinstance(c, ROOT.TCanvas): c.SaveAs("{}/{}.pdf".format(outputPath, c.GetName()))

def no_correction(hdef, name):
    h = hdef.h.Clone(name)
    return h

def br_correction(hdef, name):
    h = no_correction(hdef, name)
    h.Scale(0.0393)
    return h

def br_eff_correction(hdef, name):
    h = hdef.h_eff.Clone(name)
    h.Scale(0.0393)
    return h

def PlotInclusiveJets(hdefs):
    list_of_corrections = []

    # No corrections
    list_of_corrections.append((no_correction, "raw"))

    PlotJets(hdefs, "InclusiveJets", list_of_corrections)

def PlotDMesonJets(hdefs):
    list_of_corrections = []

    # No corrections
    list_of_corrections.append((no_correction, "raw"))

    # B.R.
    list_of_corrections.append((br_correction, "BR"))

    # B.R. + efficiency
    list_of_corrections.append((br_eff_correction, "BR_eff"))

    PlotJets(hdefs, "DMesonJets", list_of_corrections)

    # Load and plot efficiency
    eff = load_efficiency()
    c = ROOT.TCanvas("DMesonJet_Efficiency", "DMesonJet_Efficiency")
    c.cd()
    eff.Draw("hist")
    globalList.append(c)
    globalList.append(eff)

def PlotJets(hdefs, plot_name, list_of_corrections):
    # Cross sections
    histograms = []
    for hdef in hdefs:
        h = hdef.h.Clone("{}_CrossSection".format(hdef.h.GetName()))
        h.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb / (GeV/#it{c})]")
        h.Scale(1, "width")
        histograms.append(h)
    comp = DMesonJetCompare.DMesonJetCompare("{}CrossSections".format(plot_name))
    comp.fDoRatioPlot = False
    comp.fColors = comp.fColors[1:]
    comp.fMarkers = comp.fMarkers[1:]
    comp.fX1LegSpectrum = 0.25
    comp.fX2LegSpectrum = 0.80
    r = comp.CompareSpectra(histograms[0], histograms[1:])
    globalList.extend(histograms)
    globalList.extend(r)

    def ApplyCorrectionAndPlot(get_hist, suffix="raw"):
        lumi_list = [15, 50, 150]
        histograms = []
        for lumi in lumi_list:
            if lumi == 15:
                h = get_hist(hdefs[0], "{}_Lumi{}_{}".format(hdefs[0].h.GetName(), lumi, suffix))
                h.Scale(hdefs[0].acceptance / (2 * math.pi) * lumi * 1e6)
                h.SetTitle("TPC jets, #it{{L}}_{{int}} = {} nb^{{-1}}".format(lumi))
                h.GetYaxis().SetTitle("counts")
                histograms.append(h)
            h = get_hist(hdefs[1], "{}_Lumi{}_{}".format(hdefs[1].h.GetName(), lumi, suffix))
            h.Scale(hdefs[1].acceptance / (2 * math.pi) * lumi * 1e6)
            h.SetTitle("TPC+EMCal jets, #it{{L}}_{{int}} = {} nb^{{-1}}".format(lumi))
            h.GetYaxis().SetTitle("counts")
            histograms.append(h)

        comp = DMesonJetCompare.DMesonJetCompare("{}Counts_{}".format(plot_name, suffix))
        comp.fDoRatioPlot = False
        comp.fX1LegSpectrum = 0.25
        comp.fX2LegSpectrum = 0.80
        r = comp.CompareSpectra(histograms[0], histograms[1:])
        comp.fCanvasSpectra.SetGridy()
        comp.fCanvasSpectra.SetGridx()
        counts_label = ROOT.TPaveText(0.35, 0.70, 0.80, 0.52, "NB NDC")
        counts_label.SetTextFont(43)
        counts_label.SetTextSize(18)
        counts_label.SetFillStyle(0)
        counts_label.SetBorderSize(0)
        counts_label.SetTextAlign(11)
        for h in [histograms[0], histograms[-1]]:
            counts = h.Integral(h.GetXaxis().FindBin(40), h.GetXaxis().GetNbins())
            counts_label.AddText(h.GetTitle())
            counts_label.AddText("with #it{{p}}_{{T,jet}} > 40 GeV/#it{{c}}: {:.0f} counts".format(counts))
        counts_label.Draw()
        globalList.extend(histograms)
        globalList.extend(r)
        globalList.append(counts_label)

    for f, n in list_of_corrections:
        ApplyCorrectionAndPlot(f, n)

def get_weight(file, list_name):
    tlist = file.Get(list_name)
    hTrials = tlist.FindObject("fHistTrialsVsPtHardNoSel")
    hXsec = tlist.FindObject("fHistXsectionVsPtHardNoSel")
    trials = hTrials.Integral()
    xsec = hXsec.GetMean(2)
    w = xsec / trials
    return w

def main(anatype, test, gen, proc, ts, stage, dmeson_name):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    if test:
        filename = test
    else:
        filename = "/Volumes/DATA/ALICE/JetResults/FastSim_{gen}_{proc}_{ts}/stage_{stage}/output/001/AnalysisResults_FastSim_{gen}_{proc}_{ts}.root".format(gen=gen, proc=proc, ts=ts, stage=stage)
    if anatype == "djets":
        if dmeson_name == "D0":
            do_analysis(filename, GetD0MesonJets, PlotDMesonJets)
        elif dmeson_name == "DStar":
            do_analysis(filename, GetDStarMesonJets, PlotDMesonJets)
    elif anatype == "jets":
        do_analysis(filename, GetInclusiveJets, PlotInclusiveJets)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trigger studies.')
    parser.add_argument('anatype',
                        help='Analysis type: either jets or djets')
    parser.add_argument('--ts',
                        help='Timestamp of the fast simulation')
    parser.add_argument('--test',
                        default=None,
                        help='Test')
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

    main(args.anatype, args.test, args.gen, args.proc, args.ts, args.stage, args.d_meson)

    IPython.embed()
