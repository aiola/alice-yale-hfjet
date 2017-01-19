#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import RawYieldSpectrumLoader
import subprocess

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"

crossSection = 62.3  # mb CINT1
branchingRatio = 0.0393  # D0->Kpi
antiPartNorm = 2.0  # D0 / D0bar
events = 0

def LoadEvents(file, hname):
    global events
    slash = hname.find("/")
    heventsName = hname[:slash + 1]
    heventsName += "Events"
    hevents = DMesonJetUtils.GetObject(file, heventsName)
    events = hevents.GetBinContent(1)
    print("The number of events is {0}".format(events))

def GetJetPtSpectrum():
    global events
    unfolding = "LHC10_Train823_LHC15i2_Train961_efficiency_mt"
    fname = "{0}/{1}/{1}.root".format(input_path, unfolding)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    else:
        print("File {0} open successfully".format(fname))

    hname = "SideBand_DPt_30/Bayes/SideBand_DPt_30_UnfoldedSpectrum_Bayes_Reg4_PriorResponseTruth"

    LoadEvents(file, hname)
    h = DMesonJetUtils.GetObject(file, hname)
    print("The number of events is {0}".format(events))
    h.Scale(crossSection / (events * branchingRatio * antiPartNorm), "width")
    return h

def GetDPtSpectrum(kincuts=None, jet_radius=None, jet_type=None):
    global events
    loader = RawYieldSpectrumLoader.RawYieldSpectrumLoader(input_path, "Jets_EMC_pp_823_824_825_826", "LHC10_Train823_efficiency")
    loader.fUseReflections = False
    loader.fDMeson = "D0"
    loader.fSpectrumName = "DPtSpectrum"
    loader.fJetType = jet_type
    loader.fJetRadius = jet_radius
    loader.fKinematicCuts = kincuts
    h = loader.GetDefaultSpectrumFromDMesonJetAnalysis("InvMassFit")
    h.Scale(crossSection / (events * branchingRatio * antiPartNorm), "width")
    return h

def CompareJetPtvsDPt():
    jetPtSpectrumHist = GetJetPtSpectrum()
    dptSpectrumHist = GetDPtSpectrum()

    canvas = ROOT.TCanvas("Comparison_DPt_JetPt_Spectra", "Comparison_DPt_JetPt_Spectra")
    canvas.SetLogy()
    canvas.SetTicks(1, 1)
    globalList.append(canvas)
    canvas.cd()

    h = dptSpectrumHist.DrawCopy("axis")
    globalList.append(h)
    h.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [mb (GeV/#it{c})^{-1}]")
    h.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")

    dptSpectrumHist_copy = dptSpectrumHist.DrawCopy("same")
    globalList.append(dptSpectrumHist_copy)
    dptSpectrumHist_copy.SetMarkerStyle(ROOT.kOpenCircle)
    dptSpectrumHist_copy.SetMarkerSize(1.2)
    dptSpectrumHist_copy.SetMarkerColor(ROOT.kRed + 1)
    dptSpectrumHist_copy.SetLineColor(ROOT.kRed + 1)

    jetPtSpectrumHist_copy = jetPtSpectrumHist.DrawCopy("same")
    globalList.append(jetPtSpectrumHist_copy)
    jetPtSpectrumHist_copy.SetMarkerStyle(ROOT.kFullCircle)
    jetPtSpectrumHist_copy.SetMarkerSize(1.0)
    jetPtSpectrumHist_copy.SetMarkerColor(ROOT.kBlue + 1)
    jetPtSpectrumHist_copy.SetLineColor(ROOT.kBlue + 1)

    pave = ROOT.TPaveText(0.19, 0.80, 0.60, 0.90, "NB NDC")
    globalList.append(pave)
    pave.SetBorderSize(0)
    pave.SetFillStyle(0)
    pave.SetTextFont(43)
    pave.SetTextSize(20)
    pave.SetTextAlign(13)
    pave.AddText("pp, #sqrt{#it{s}} = 7 TeV")
    pave.AddText("Fully corrected spectra, statistical uncertainty only")
    pave.Draw()

    leg = ROOT.TLegend(0.35, 0.50, 0.81, 0.75, "", "NB NDC")
    globalList.append(leg)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(20)
    leg.SetTextAlign(13)
    leg.SetMargin(0.2)
    leg.AddEntry(jetPtSpectrumHist_copy, "Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4", "pe")
    leg.AddEntry(ROOT.nullptr, "|#eta_{jet}| < 0.5", "")
    leg.AddEntry(ROOT.nullptr, "with D^{0}, #it{p}_{T,D} > 3 GeV/#it{c}", "")
    leg.AddEntry(dptSpectrumHist_copy, "D^{0}, |#eta_{D}| < 0.5", "pe")
    leg.Draw()

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    CompareJetPtvsDPt()

    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{0}/{1}.pdf".format(input_path, obj.GetName()))

if __name__ == '__main__':
    main()

    IPython.embed()
