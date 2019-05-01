#!/usr/bin/env python

import subprocess
import IPython
import ROOT
import RawYieldSpectrumLoader
import SideBandInvMass_Pt_Paper
import SideBandInvMass_Z_Paper

globalList = []

input_path = "../../workdir"

def SideBandPlot():
    loader_pt = RawYieldSpectrumLoader.RawYieldSpectrumLoader(input_path, "Jets_EMC_pp_1116_1117_1118_1119", "LHC10_Train1116_efficiency")
    loader_pt.fDMeson = "D0_D0toKpiCuts"
    loader_pt.fJetType = "Charged"
    loader_pt.fJetRadius = "R040"
    loader_pt.fVariableName = "JetPt"
    loader_pt.fKinematicCuts = "DPt_30"
    loader_pt.fRawYieldMethod = "SideBand"
    loader_pt.fJetRecoScheme = "pt_scheme"
    loader_pt.fUseReflections = True
    loader_pt.LoadDataListFromDMesonJetAnalysis()
    dptbinlist_pt_name = "{}_Charged_R040_pt_scheme_DPtBins_JetPt_5_30".format(loader_pt.fDMeson)
    dptbinList_pt = loader_pt.fDataJetList.FindObject(dptbinlist_pt_name)
    if not dptbinList_pt:
        print("Could not find '{}' in '{}'".format(dptbinlist_pt_name, loader_pt.fDataJetList.GetName()))
        loader_pt.fDataJetList.Print()
        exit(1)
    spectrumList_pt = loader_pt.fDataSpectrumList
    sbList_pt = spectrumList_pt.FindObject("SideBandAnalysis")

    loader_high = RawYieldSpectrumLoader.RawYieldSpectrumLoader(input_path, "Jets_EMC_pp_1116_1117_1118_1119", "LHC10_Train1116_efficiency")
    loader_high.fDMeson = "D0_D0toKpiCuts"
    loader_high.fJetType = "Charged"
    loader_high.fJetRadius = "R040"
    loader_high.fVariableName = "JetZ"
    loader_high.fRawYieldMethod = "SideBand"
    loader_high.fJetRecoScheme = "pt_scheme"
    loader_high.fUseReflections = True
    loader_high.fKinematicCuts = "DPt_60_JetPt_15_30"
    loader_high.LoadDataListFromDMesonJetAnalysis()
    dptbinList_high_name = "{}_Charged_R040_pt_scheme_DPtBins_JetPt_15_30".format(loader_high.fDMeson)
    dptbinList_high = loader_high.fDataJetList.FindObject(dptbinList_high_name)
    if not dptbinList_high:
        print("Could not find '{}' in '{}'".format(dptbinList_high_name, loader_high.fDataJetList.GetName()))
        loader_high.fDataJetList.Print()
        exit(1)
    spectrumList_high = loader_high.fDataSpectrumList
    sbList_high = spectrumList_high.FindObject("SideBandAnalysis")

    cname = "SideBandInvMass_Paper"
    canvas = ROOT.TCanvas(cname, cname, 1200, 800)
    globalList.append(canvas)
    canvas.Divide(3, 2)
    histograms_to_output = []
    histograms_to_output += SideBandInvMass_Pt_Paper.PlotSBInvMass(canvas.cd(1), loader_pt.fDMeson, 4, 5, 5, 30, sbList_pt, dptbinList_pt, "_DoubleGaus", True, False)
    histograms_to_output += SideBandInvMass_Pt_Paper.PlotSBInvMass(canvas.cd(2), loader_pt.fDMeson, 6, 7, 5, 30, sbList_pt, dptbinList_pt, "_DoubleGaus", False, True)
    histograms_to_output += SideBandInvMass_Z_Paper.PlotSBInvMass(canvas.cd(3), loader_high.fDMeson, 6, 12, 15, 30, sbList_high, dptbinList_high, "_DoubleGaus", False, False)
    histograms_to_output += SideBandInvMass_Pt_Paper.PlotSBSpectra(canvas.cd(4), loader_pt.fDMeson, 4, 5, sbList_pt, "_DoubleGaus")
    histograms_to_output += SideBandInvMass_Pt_Paper.PlotSBSpectra(canvas.cd(5), loader_pt.fDMeson, 6, 7, sbList_pt, "_DoubleGaus")
    histograms_to_output += SideBandInvMass_Z_Paper.PlotSBSpectra(canvas.cd(6), loader_high.fDMeson, "DPt_60_JetPt_15_30", 6, 12, sbList_high, "_DoubleGaus", True)

    canvas.cd(5)
    htitle = ROOT.TPaveText(0.44, 0.56, 0.92, 0.95, "NB NDC")
    globalList.append(htitle)
    htitle.SetBorderSize(0)
    htitle.SetFillStyle(0)
    htitle.SetTextFont(43)
    htitle.SetTextSize(19)
    htitle.SetTextAlign(11)
    htitle.AddText("Charged Jets")
    htitle.AddText("Anti-#it{k}_{T}, #it{R} = 0.4")
    htitle.AddText("|#it{#eta}_{jet}| < 0.5")
    htitle.AddText("with D^{0} #rightarrow K^{#font[122]{-}}#pi^{+}")
    htitle.AddText("and charge conj.")
    htitle.Draw()

    canvas.cd(4)
    paveALICE = ROOT.TPaveText(0.30, 0.84, 0.66, 0.95, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(20)
    paveALICE.SetTextAlign(11)
    paveALICE.AddText("ALICE, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.Draw()

    return histograms_to_output

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")
    print("MassFitter loaded")

    histograms_to_output = SideBandPlot()

    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{0}/{1}.pdf".format(input_path, obj.GetName()))
            obj.SaveAs("{0}/{1}.C".format(input_path, obj.GetName()))

    output_file = ROOT.TFile("{}/SideBandInvMass_Paper.root".format(input_path), "recreate")
    output_file.cd()
    for h in histograms_to_output:
        h.Write()
    output_file.Close()

if __name__ == '__main__':
    main()

    IPython.embed()
