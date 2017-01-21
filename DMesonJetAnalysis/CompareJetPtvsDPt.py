#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import RawYieldSpectrumLoader
import subprocess
import DMesonJetCompare

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
    loader.fRawYieldMethod = "InvMassFit"
    h = loader.GetDefaultSpectrumFromDMesonJetAnalysis()
    h.Scale(crossSection / (events * branchingRatio * antiPartNorm), "width")
    return h

def GetTheoryJetPtCrossSection():
    fname = "{0}/PromptDJetsPrediction.root".format(input_path)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "default/JetPtSpectrum_DPt_30/GeneratorLevel_JetPtSpectrum")
    if not hStat:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)
    hSystUp = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetPtSpectrum_DPt_30//GeneratorLevel_JetPtSpectrum/GeneratorLevel_JetPtSpectrum_UpperSyst")
    if not hSystUp:
        print("Cannot get theory cross section upper systematic uncertainty!")
        exit(1)
    hSystLow = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetPtSpectrum_DPt_30//GeneratorLevel_JetPtSpectrum/GeneratorLevel_JetPtSpectrum_LowerSyst")
    if not hSystUp:
        print("Cannot get theory cross section lower systematic uncertainty!")
        exit(1)

    # scale for the bin width and the antiparticle factor
    hStat.Scale(0.5, "width")
    hSystUp.Scale(0.5, "width")
    hSystLow.Scale(0.5, "width")
    return hStat, hSystUp, hSystLow

def GetTheoryDPtCrossSection():
    fname = "{0}/PromptDJetsPrediction.root".format(input_path)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "default/DPtSpectrum/GeneratorLevel_DPtSpectrum")
    if not hStat:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)
    hSystUp = DMesonJetUtils.GetObject(file, "SystematicUncertainty/DPtSpectrum//GeneratorLevel_DPtSpectrum/GeneratorLevel_DPtSpectrum_UpperSyst")
    if not hSystUp:
        print("Cannot get theory cross section upper systematic uncertainty!")
        exit(1)
    hSystLow = DMesonJetUtils.GetObject(file, "SystematicUncertainty/DPtSpectrum//GeneratorLevel_DPtSpectrum/GeneratorLevel_DPtSpectrum_LowerSyst")
    if not hSystUp:
        print("Cannot get theory cross section lower systematic uncertainty!")
        exit(1)

    # scale for the bin width and the antiparticle factor
    hStat.Scale(0.5, "width")
    hSystUp.Scale(0.5, "width")
    hSystLow.Scale(0.5, "width")
    return hStat, hSystUp, hSystLow

def CompareJetPtvsDPt():
    jetPtSpectrumHist = GetJetPtSpectrum()
    dptSpectrumHist = GetDPtSpectrum()
    jetPtSpectrumHist.SetTitle("Data, D^{0}-jet")
    dptSpectrumHist.SetTitle("Data, D^{0}")

    jetPtSpectrumHist_theory = GetTheoryJetPtCrossSection()[0]
    dptSpectrumHist_theory = GetTheoryDPtCrossSection()[0]
    jetPtSpectrumHist_theory.SetTitle("POWHEG, D^{0}-jet")
    dptSpectrumHist_theory.SetTitle("POWHEG, D^{0}")

    canvasSpectrum = ROOT.TCanvas("Comparison_DPt_JetPt_Spectra")
    canvasSpectrum.SetTicks(1, 1)
    canvasSpectrum.cd()
    h_axis = dptSpectrumHist.DrawCopy("axis")
    h_axis.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [mb (GeV/#it{c})^{-1}]")
    h_axis.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")

    comp = DMesonJetCompare.DMesonJetCompare("Comparison_DPt_JetPt_Spectra")
    comp.fCanvasSpectra = canvasSpectrum
    comp.fOptSpectrumBaseline = "same"
    comp.fOptSpectrum = "same"
    comp.fX1LegSpectrum = 0.58
    comp.fY1LegSpectrum = 0.52

    comp.fMarkers = [ROOT.kFullCircle, ROOT.kFullSquare]
    comp.fColors = [ROOT.kRed + 2, ROOT.kBlue + 2]
    comp.CompareSpectra(jetPtSpectrumHist, [dptSpectrumHist])

    comp.fMarkers = [ROOT.kOpenCircle, ROOT.kOpenSquare]
    comp.fColors = [ROOT.kOrange + 2, ROOT.kGreen + 2]
    r = comp.CompareSpectra(jetPtSpectrumHist_theory, [dptSpectrumHist_theory])
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

    canvasSpectrum.cd()

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

    pave = ROOT.TPaveText(0.35, 0.55, 0.81, 0.75, "NB NDC")
    globalList.append(pave)
    pave.SetBorderSize(0)
    pave.SetFillStyle(0)
    pave.SetTextFont(43)
    pave.SetTextSize(20)
    pave.SetTextAlign(13)
    pave.AddText("Charged D^{0}-jets: Anti-#it{k}_{T}, #it{R}=0.4")
    pave.AddText("|#eta_{jet}| < 0.5 with D^{0}, #it{p}_{T,D} > 3 GeV/#it{c}")
    pave.AddText("D^{0}: |#eta_{D}| < 0.5")
    pave.Draw()

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
