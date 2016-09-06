#!/usr/bin/env python
#python script to generate plots requested for approval to be presented at HQ16 conference

import argparse
import yaml
import IPython
import ROOT
import subprocess

globalList = []
canvases = []

def main(actions, output_path, output_type):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    f = open("LHC14j4analysis_Train953.yaml", 'r')
    config = yaml.load(f)
    f.close()
    
    f = open("LHC15i2analysis_Train949.yaml", 'r')
    configRM = yaml.load(f)
    f.close()
    
    f = open("LHC10analysis_Train764.yaml", 'r')
    configData = yaml.load(f)
    f.close()

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    file = OpenFile(config)
    fileW = OpenFile(config, "efficiency")
    fileRM = OpenFile(configRM)
    fileData = OpenFile(configData)
    fileWData = OpenFile(configData, "efficiency")

    if "all" in actions or "invmass" in actions:
        InvMassPlots(file, config)

    if "all" in actions or "compare" in actions:
        CompareMethods(file, config)

    if "all" in actions or "unc" in actions:
        CompareUncertainties(file, fileW, config, "")
        #CompareUncertainties(fileData, fileWData, config, "Data")

    if "all" in actions or "eff" in actions:
        EfficiencyPlots(fileRM, configRM)

    if "all" in actions or "resp" in actions:
        DetectorResponsePlots(fileRM, configRM)

    if "all" in actions or "stat" in actions:
        StatisticalUncertaintyData(fileData, configData)
        
    if "all" in actions or "spectra" in actions:
        PlotSpectra(file)
        PlotSpectra(fileW, "Eff")

    for c in canvases:
        c.SaveAs("{0}/{1}.{2}".format(output_path, c.GetName(), output_type))

def StatisticalUncertaintyData(file, config):
    dataList = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum", file)
    hist = dataList["D0_D_Tagged_Jet_PtD_20_Spectrum_Unc"]
    cname = "HQ16_WorkInProgress_StatisticalUncertainty"
    c = ROOT.TCanvas(cname, cname, 650, 500)
    c.SetTicks(1,1)
    canvases.append(c)
    c.cd()
    h = hist.DrawCopy("hist")

    h.SetLineColor(ROOT.kBlue+2)
    h.SetLineWidth(2)
    h.SetFillColor(ROOT.kBlue-10)
    h.SetFillStyle(1001)
    h.GetYaxis().SetRangeUser(0, 0.25)
    h.GetYaxis().SetTitleOffset(1.2)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetXaxis().SetTitle("#it{p}_{T,ch jet}^{det} (GeV/#it{c})")

    pave = ROOT.TPaveText(0.12, 0.88, 0.4, 0.55, "NB NDC")
    pave.SetFillStyle(0)
    pave.SetBorderSize(0)
    pave.SetTextFont(43)
    pave.SetTextSize(15)
    pave.AddText("ALICE Work In Progress")
    pave.AddText("pp, #sqrt{#it{s}} = 7 TeV")
    pave.AddText("316 M minimum-bias events")
    pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4")
    pave.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
    pave.AddText("#it{p}_{T,D} > 2 GeV/#it{c}")
    pave.Draw()

    globalList.append(c)
    globalList.append(h)
    globalList.append(pave)

def DetectorResponsePlots(file, config):
    simuPlot = ["ALICE Simulation", "PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV",
                "Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4, |#eta_{jet}| < 0.5",
                "with D^{0} #rightarrow K^{-}#pi^{+} and c.c.",
                "2 < #it{p}_{T,D} < 24 GeV/#it{c}"]
    respList = LoadHistograms("D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum", file)
    resp = respList["DetectorResponse"]
    histList = [resp["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_5_6"], resp["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_8_10"], 
                resp["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_18_24"]]
    histList[0].SetTitle("5 < #it{p}_{T,ch jet}^{part} < 6 GeV/#it{c}")
    histList[1].SetTitle("8 < #it{p}_{T,ch jet}^{part} < 10 GeV/#it{c}")
    histList[2].SetTitle("18 < #it{p}_{T,ch jet}^{part} < 24 GeV/#it{c}")
    (blank, c) = PlotMultiHistogram(histList, "HQ16_Simulation_DetectorResponse", "#(){#it{p}_{T,ch jet}^{det} #minus #it{p}_{T,ch jet}^{part}} / #it{p}_{T,ch jet}^{part}", "Probability density",
                               [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2], [ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenCross], [1.0, 1.0, 1.4], simuPlot)
    blank.GetYaxis().SetRangeUser(0, 18)
    blank.GetXaxis().SetRangeUser(-0.8, 0.8)
    blank.GetXaxis().SetTitleOffset(1.4)
    c.SetBottomMargin(0.15)
    
    histList = [respList["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_EnergyScaleShift"], respList["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_EnergyScaleShiftMedian"]]
    histList[0].SetTitle("Mean")
    histList[1].SetTitle("Median")
    (blank, c) = PlotMultiHistogram(histList, "HQ16_Simulation_EnergyScaleShift", "#it{p}_{T,ch jet}^{part} (GeV/#it{c})", "#(){#it{p}_{T,ch jet}^{det} #minus #it{p}_{T,ch jet}^{part}} / #it{p}_{T,ch jet}^{part}",
                               [ROOT.kBlue+2, ROOT.kRed+2], [ROOT.kFullCircle, ROOT.kFullSquare], [1.0, 1.0], simuPlot)
    blank.GetYaxis().SetRangeUser(-0.07, 0.08)
    blank.GetXaxis().SetTitleOffset(1.4)
    blank.GetYaxis().SetTitleOffset(1.4)

    histList = [respList["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_Resolution"]]
    histList[0].SetTitle("Resolution")
    (blank, c) = PlotMultiHistogram(histList, "HQ16_Simulation_Resolution", "#it{p}_{T,ch jet}^{part} (GeV/#it{c})", "#sigma#[]{#(){#it{p}_{T,ch jet}^{det} #minus #it{p}_{T,ch jet}^{part}} / #it{p}_{T,ch jet}^{part}}",
                               [ROOT.kBlue+2], [ROOT.kFullCircle], [1.0], simuPlot)
    blank.GetYaxis().SetRangeUser(0.05, 0.18)
    blank.GetXaxis().SetTitleOffset(1.4)
    blank.GetYaxis().SetTitleOffset(1.4)

def EfficiencyPlots(file, config):
    simuPlot = ["ALICE Simulation", "PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV",
                "Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4, |#eta_{jet}| < 0.5",
                "with D^{0} #rightarrow K^{-}#pi^{+} and c.c."]
    ptEff = LoadHistograms("D0_Jet_AKTChargedR040_pt_scheme_D_Spectra", file)
    histList = [ptEff["D0_Jet_AKTChargedR040_pt_scheme_D_Spectra_Efficiency_JetPt_500_2400"], ptEff["D0_Jet_AKTChargedR040_pt_scheme_D_Spectra_Efficiency_JetPt_500_800"], 
                ptEff["D0_Jet_AKTChargedR040_pt_scheme_D_Spectra_Efficiency_JetPt_800_1300"], ptEff["D0_Jet_AKTChargedR040_pt_scheme_D_Spectra_Efficiency_JetPt_1300_2400"]]
    histList[0].SetTitle("5 < #it{p}_{T,ch jet} < 24 GeV/#it{c}")
    histList[1].SetTitle("5 < #it{p}_{T,ch jet} < 8 GeV/#it{c}")
    histList[2].SetTitle("8 < #it{p}_{T,ch jet} < 13 GeV/#it{c}")
    histList[3].SetTitle("13 < #it{p}_{T,ch jet} < 24 GeV/#it{c}")
    PlotMultiHistogram(histList, "HQ16_Simulation_EfficiencyVsDPt", "#it{p}_{T,D} (GeV/#it{c})", "D-tagged Jet Efficiency",
                       [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2], [ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenCross, ROOT.kOpenStar], [1.0, 1.0, 1.4, 1.5], simuPlot)
    (blank,c) = PlotMultiHistogram(histList, "HQ16_Simulation_EfficiencyVsDPt_LogScale", "#it{p}_{T,D} (GeV/#it{c})", "D-tagged Jet Efficiency",
                       [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2], [ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenCross, ROOT.kOpenStar], [1.0, 1.0, 1.4, 1.5], None,True)
    blank.GetYaxis().SetRangeUser(2e-3, 9)

def PlotSpectra(file, suffix=""):
    if suffix:
        suffix = "_{0}".format(suffix)

    spectrumLSlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign", file)
    spectrumSBlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand", file)

    #bins = ["200_300", "300_400", "400_500", "500_600", "600_700", "700_800", "800_1200", "1200_1600", "1600_2400"]
    bins = ["200_300", "500_600", "700_800", "1200_1600", "1600_2400"]

    hlist = [spectrumSBlist["SideBandAnalysis"]["D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand_SideBandWindow_DPt_{0}".format(bin)] for bin in bins]
    PlotMultiHistogram(hlist, "HQ16_Simulation_SBSpectra{0}".format(suffix), "#it{p}_{T,ch jet}^{det} (GeV/#it{c})", "counts",
                       [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kAzure+2, ROOT.kMagenta+2, ROOT.kCyan+2, ROOT.kPink+2, ROOT.kTeal+2], [ROOT.kOpenCircle]*9, [1.0]*9, None,True)
    
    hlist = [spectrumSBlist["SideBandAnalysis"]["D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand_SignalWindow_DPt_{0}".format(bin)] for bin in bins]
    PlotMultiHistogram(hlist, "HQ16_Simulation_SigSpectra{0}".format(suffix), "#it{p}_{T,ch jet}^{det} (GeV/#it{c})", "counts",
                       [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kAzure+2, ROOT.kMagenta+2, ROOT.kCyan+2, ROOT.kPink+2, ROOT.kTeal+2], [ROOT.kOpenSquare]*9, [1.0]*9, None,True)

    hlist = [spectrumLSlist["LikeSignAnalysis"]["D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign_LikeSign_DPt_{0}".format(bin)] for bin in bins]
    PlotMultiHistogram(hlist, "HQ16_Simulation_LSSpectra{0}".format(suffix), "#it{p}_{T,ch jet}^{det} (GeV/#it{c})", "counts",
                       [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kAzure+2, ROOT.kMagenta+2, ROOT.kCyan+2, ROOT.kPink+2, ROOT.kTeal+2], [ROOT.kOpenCross]*9, [1.4]*9, None,True)

def PlotMultiHistogram(histList, cname, xaxisTitle, yaxisTitle, colors=None, markers=None, markerSizes=None, simuPlot=None, logY=False):
    c = ROOT.TCanvas(cname, cname)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.08)
    c.SetRightMargin(0.08)
    c.SetTicks(1,1)
    c.cd()
    if logY:
        c.SetLogy()
    globalList.append(c)
    canvases.append(c)
    blank = ROOT.TH1D("blankHist", "blankHist;{0};{1}".format(xaxisTitle, yaxisTitle), 100, histList[0].GetXaxis().GetXmin(), histList[0].GetXaxis().GetXmax())
    blank.GetXaxis().SetTitleFont(43)
    blank.GetXaxis().SetTitleOffset(1.2)
    blank.GetXaxis().SetTitleSize(19)
    blank.GetXaxis().SetLabelFont(43)
    blank.GetXaxis().SetLabelOffset(0.009)
    blank.GetXaxis().SetLabelSize(18)
    blank.GetYaxis().SetTitleFont(43)
    blank.GetYaxis().SetTitleOffset(1.2)
    blank.GetYaxis().SetTitleSize(19)
    blank.GetYaxis().SetLabelFont(43)
    blank.GetYaxis().SetLabelOffset(0.009)
    blank.GetYaxis().SetLabelSize(18)
    blank.Draw("AXIS")
    globalList.append(blank)
    if not colors:
        colors = [ROOT.kBlack, ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kMagenta+2]
    if not markers:
        markers = [ROOT.kStar, ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullDiamond, ROOT.kFullCross, ROOT.kFullStar]
    if not markerSizes:
        markerSizes = [1.0]*len(markers)
    max = 0;
    min = 1e15
    leg = ROOT.TLegend(0.55, 0.90-len(histList)*0.075, 0.90, 0.85)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(19)
    leg.SetMargin(0.15)
    for color,marker,markerSize,eff in zip(colors,markers,markerSizes,histList):
        h = eff.Clone()
        globalList.append(h)
        h.SetMarkerStyle(marker)
        h.SetMarkerSize(markerSize)
        h.SetMarkerColor(color)
        h.SetLineColor(color)
        leg.AddEntry(h, h.GetTitle(), "pe")
        h.Draw("same")
        for i in range(1, h.GetNbinsX()+1):
            y = h.GetBinContent(i)
            if y > max:
                max = y
            if y < max:
                min = y

    if len(histList) > 1:
        leg.Draw()
    globalList.append(leg)
    if logY:
        blank.GetYaxis().SetRangeUser(min/10, max*10)
    else:
        blank.GetYaxis().SetRangeUser(0, max*1.8)

    if simuPlot:
        paveALICE = ROOT.TPaveText(0.13, 0.64, 0.52, 0.90, "NB NDC")
        globalList.append(paveALICE)
        paveALICE.SetBorderSize(0)
        paveALICE.SetFillStyle(0)
        paveALICE.SetTextFont(43)
        paveALICE.SetTextSize(17)
        paveALICE.SetTextAlign(13)
        for t in simuPlot:
            paveALICE.AddText(t)
        paveALICE.Draw()
    return blank, c

def CompareUncertainties(file, fileW, config, suffix):
    LSlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign", file)
    SBlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand", file)
    IMlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum", file)
    LSlist_eff = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign", fileW)
    SBlist_eff = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand", fileW)
    IMlist_eff = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum", fileW)

    #unc = [IMlist["D0_D_Tagged_Jet_PtD_20_Spectrum_Unc"], SBlist["D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand_Unc"], LSlist["D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign_Unc"],
    #       IMlist_eff["D0_D_Tagged_Jet_PtD_20_Spectrum_Unc"], SBlist_eff["D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand_Unc"], LSlist_eff["D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign_Unc"]]
    #colors = [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2,
    #          ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2]
    #titles = ["Inv.Mass Fit", "Side-band Method", "Like-sign Method",
    #          "Inv.Mass Fit w/ eff.", "Side-band Method w/ eff.", "Like-sign Method w/ eff."]
    #lineStyles = [1,1,1,
    #              2,2,2]

    unc = [IMlist["D0_D_Tagged_Jet_PtD_20_Spectrum_Unc"], LSlist["D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign_Unc"],
           IMlist_eff["D0_D_Tagged_Jet_PtD_20_Spectrum_Unc"]]
    colors = [ROOT.kBlue+2, ROOT.kGreen+2,
              ROOT.kBlue+2]
    titles = ["Inv.Mass Fit", "Like-sign Method",
              "Inv.Mass Fit w/ eff."]
    lineStyles = [1,1,
                  2]
    if suffix:
        cname = "HQ16_Simulation_UncertaintyComparison_{0}".format(suffix)
    else:
        cname = "HQ16_Simulation_UncertaintyComparison"
    c = ROOT.TCanvas(cname, cname, 650, 500)
    c.SetTicks(1,1)
    c.SetBottomMargin(0.12)
    canvases.append(c)
    globalList.append(c)

    h = unc[0].DrawCopy("AXIS")
    globalList.append(h)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(19)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(18)
    h.GetYaxis().SetTitleOffset(1.1)
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleSize(19)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelSize(18)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetYaxis().SetRangeUser(0, 0.50)
    h.GetXaxis().SetTitle("#it{p}_{T,ch jet}^{det} (GeV/#it{c})")

    leg1 = ROOT.TLegend(0.13, 0.45, 0.55, 0.58, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(19)
    leg1.SetTextAlign(13)
    leg1.SetMargin(0.1)

    for s,color,lineStyle,title in zip(unc,colors,lineStyles,titles):
        h = s.DrawCopy("same hist")
        globalList.append(h)
        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.SetLineStyle(lineStyle)
        h.SetFillStyle(0)
        leg1.AddEntry(h, title, "l")

    paveALICE = ROOT.TPaveText(0.11, 0.69, 0.50, 0.90, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(17)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Simulation")
    paveALICE.AddText("PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("350 M minimum-bias events")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4, |#eta_{jet}| < 0.5 with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
    paveALICE.AddText("2 < #it{p}_{T,D} < 24 GeV/#it{c}")
    paveALICE.Draw()

    leg1.Draw()

def CompareMethods(file, config):
    Tlist = LoadHistograms("D0_kSignalOnly_D_Tagged_Jet_PtD_20_Spectrum", file)
    LSlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign", file)
    SBlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand", file)
    IMlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum", file)

    Tspectrum = Tlist["D0_kSignalOnly_D_Tagged_Jet_PtD_20_Spectrum"]
    spectra = [IMlist["D0_D_Tagged_Jet_PtD_20_Spectrum"], SBlist["D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand"], LSlist["D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign"]]

    for ibin in range(0,Tspectrum.GetNbinsX()+2):
        Tspectrum.SetBinError(ibin,0)

    cname = "HQ16_Simulation_MethodComparison_Diff"
    cDiff = ROOT.TCanvas(cname, cname, 650, 500)
    canvases.append(cDiff)

    cname = "HQ16_Simulation_MethodComparison"
    c = ROOT.TCanvas(cname, cname, 700, 700)
    canvases.append(c)
    c.Divide(1, 2)
    globalList.append(c)
    padMain = c.cd(1)
    padMain.SetPad(0, 0.35, 1, 1)
    padMain.SetBottomMargin(0)
    padMain.SetTicks(1,1)
    hS = Tspectrum.DrawCopy("hist")
    globalList.append(hS)
    hS.SetLineColor(ROOT.kGray)
    hS.SetFillColor(ROOT.kGray)
    hS.SetFillStyle(1001)
    hS.GetYaxis().SetRangeUser(0.001,1100)
    hS.GetYaxis().SetTitleFont(43)
    hS.GetYaxis().SetTitleSize(19)
    hS.GetYaxis().SetLabelFont(43)
    hS.GetYaxis().SetLabelSize(18)
    hS.GetYaxis().SetTitleOffset(1.6)

    padRatio = c.cd(2)
    padRatio.SetPad(0, 0., 1, 0.35)
    padRatio.SetTopMargin(0)
    padRatio.SetBottomMargin(0.25)
    padRatio.SetTicks(1,1)
    hSratio = Tspectrum.DrawCopy("hist")
    globalList.append(hSratio)
    hSratio.Divide(Tspectrum)
    hSratio.SetFillStyle(0)
    hSratio.SetLineColor(ROOT.kBlack)
    hSratio.SetLineStyle(2)
    hSratio.SetLineWidth(2)
    hSratio.GetXaxis().SetTitleFont(43)
    hSratio.GetXaxis().SetTitleSize(19)
    hSratio.GetXaxis().SetLabelFont(43)
    hSratio.GetXaxis().SetLabelSize(18)
    hSratio.GetYaxis().SetTitleFont(43)
    hSratio.GetYaxis().SetTitleSize(19)
    hSratio.GetYaxis().SetLabelFont(43)
    hSratio.GetYaxis().SetLabelSize(18)
    hSratio.GetYaxis().SetTitleOffset(1.6)
    hSratio.GetXaxis().SetTitleOffset(3.4)
    hSratio.GetXaxis().SetTitle("#it{p}_{T,ch jet}^{det} (GeV/#it{c})")
    hSratio.GetYaxis().SetTitle("ratio")
    hSratio.GetYaxis().SetNdivisions(504)

    cDiff.cd()
    hSdiff = Tspectrum.DrawCopy("hist")
    globalList.append(hSdiff)
    hSdiff.Add(Tspectrum, -1)
    hSdiff.SetFillStyle(0)
    hSdiff.SetLineColor(ROOT.kBlack)
    hSdiff.SetLineStyle(2)
    hSdiff.SetLineWidth(2)
    hSdiff.GetYaxis().SetRangeUser(-100,100)
    hSdiff.GetXaxis().SetTitle("#it{p}_{T,ch jet}^{det} (GeV/#it{c})")
    hSdiff.GetYaxis().SetTitle("difference")

    colors = [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2]
    markers = [ROOT.kFullSquare, ROOT.kFullCircle, ROOT.kFullCross]
    sizes = [1.2, 1.1, 1.2]
    fillStyles = [3002, 3245, 3254]
    titles = ["Inv.Mass Fit", "Side-band Method", "Like-sign Method"]
    
    leg1 = ROOT.TLegend(0.52, 0.28, 0.96, 0.49, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(19)
    leg1.SetTextAlign(13)
    leg1.AddEntry(hS, "MC truth", "f")

    for s,color,marker,fillStyle,title,size in zip(spectra,colors,markers,fillStyles,titles,sizes):
        padMain.cd()
        h = s.DrawCopy("same E0 X0")
        globalList.append(h)
        h.SetMarkerColor(color)
        h.SetMarkerStyle(marker)
        h.SetMarkerSize(size)
        h.SetLineColor(color)
        
        padRatio.cd()
        h = s.DrawCopy("same E2")
        globalList.append(h)
        h.Divide(Tspectrum)
        (chi2, chi2Avg) = CalculateChi2(h)
        print("The chi2 for {0} is {1}. For the average is {2}".format(title, chi2, chi2Avg))
        h.SetMarkerColor(color)
        h.SetMarkerStyle(marker)
        h.SetMarkerSize(size)
        h.SetLineColor(color)
        h.SetFillColor(color-9)
        h.SetFillStyle(fillStyle)

        leg1.AddEntry(h, title, "pef")

        cDiff.cd()
        h = s.DrawCopy("same hist")
        globalList.append(h)
        h.Add(Tspectrum, -1)
        h.SetLineColor(color)
        h.SetFillStyle(0)

    padMain.cd()
    paveALICE = ROOT.TPaveText(0.24, 0.55, 0.63, 0.88, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(17)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Simulation")
    paveALICE.AddText("PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4, |#eta_{jet}| < 0.5")
    paveALICE.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
    paveALICE.AddText("2 < #it{p}_{T,D} < 24 GeV/#it{c}")
    paveALICE.Draw()

    leg1.Draw()

def CalculateChi2(h):
    chi2 = 0
    avg = 0
    avgErr2 = 0
    for xbin in range(1, h.GetNbinsX()+1):
        chi2 += (h.GetBinContent(xbin) - 1)**2 / h.GetBinError(xbin)**2
        avg += h.GetBinContent(xbin)
        avgErr2 += h.GetBinError(xbin)**2
    avg /= h.GetNbinsX()
    avgErr2 /= h.GetNbinsX()**2
    chi2Avg = (avg - 1)**2 / avgErr2
    return chi2, chi2Avg

def InvMassPlots(file, config):
    spectrumSB = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand", file)
    spectrumLS = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign", file)
    invMass = LoadHistograms("D0_DPtBins_PtD_20", file)
    invMassLS = LoadHistograms("2ProngLikeSign_DPtBins_PtD_20", file)
    
    PlotInvMass(invMass, invMassLS, spectrumSB, spectrumLS, config)

def PlotInvMass(invMassPlotList, invMassPlotListLS, spectrumPlotListSB, spectrumPlotListLS, config):
    #bins = ["200_300", "600_700", "1200_1600"]
    bins = ["600_700"]
    #binTitles = ["2 < #it{p}_{T,D} < 3 GeV/#it{c}", "6 < #it{p}_{T,D} < 7 GeV/#it{c}", "12 < #it{p}_{T,D} < 16 GeV/#it{c}"]
    binTitles = ["6 < #it{p}_{T,D} < 7 GeV/#it{c}"]
    cname = "HQ16_Simulation_InvMassSB"
    c = ROOT.TCanvas(cname, cname, len(bins)*400, 400)
    canvases.append(c)
    c.Divide(len(bins), 1)
    globalList.append(c)
    for i,(bin,binTitle) in enumerate(zip(bins, binTitles)):
        invMassHisto = invMassPlotList["InvMass_D0_DPt_{0}".format(bin)]
        invMassHistoLS = invMassPlotListLS["InvMass_2ProngLikeSign_DPt_{0}".format(bin)]
        invMassFitter = invMassPlotList["InvMass_D0_DPt_{0}_fitter".format(bin)]
        invMassHistoSB = spectrumPlotListSB["SideBandAnalysis"]["InvMassSBWindow_D0_DPt_{0}".format(bin)]
        invMassHistoSig = spectrumPlotListSB["SideBandAnalysis"]["InvMassSigWindow_D0_DPt_{0}".format(bin)]

        pad = c.cd(i+1)
        pad.SetLeftMargin(0.17)
        pad.SetRightMargin(0.05)
        pad.SetTopMargin(0.09)
        pad.SetBottomMargin(0.13)
        pad.SetTicks(1,1)

        h = invMassHisto.DrawCopy("AXIS")
        h.GetXaxis().SetTitle("#it{M}_{K#pi} (GeV/#it{c}^{2})")
        h.GetYaxis().SetTitle("arb. units")

        (h1, h2, hsig) = PlotInvMassSideBands(invMassHistoSB, invMassHistoSig)

        invMassHisto_copy = invMassHisto.DrawCopy("same")
        globalList.append(invMassHisto_copy)
        
        invMassHistoLS_copy = invMassHistoLS.DrawCopy("same hist")
        invMassHistoLS_copy.SetFillStyle(0)
        invMassHistoLS_copy.SetLineColor(ROOT.kGreen+2)
        invMassHistoLS_copy.SetLineStyle(1)
        invMassHistoLS_copy.SetLineWidth(2)

        h.SetMaximum(invMassHisto_copy.GetMaximum()*1.3)
        h.GetXaxis().SetTitleFont(43)
        h.GetXaxis().SetTitleOffset(1.1)
        h.GetXaxis().SetTitleSize(19)
        h.GetXaxis().SetLabelFont(43)
        h.GetXaxis().SetLabelOffset(0.009)
        h.GetXaxis().SetLabelSize(18)
        h.GetYaxis().SetTitleFont(43)
        h.GetYaxis().SetTitleOffset(1.5)
        h.GetYaxis().SetTitleSize(19)
        h.GetYaxis().SetLabelFont(43)
        h.GetYaxis().SetLabelOffset(0.009)
        h.GetYaxis().SetLabelSize(18)
        htitle = ROOT.TPaveText(0.61, 0.71, 0.96, 0.89, "NB NDC")
        htitle.SetBorderSize(0)
        htitle.SetFillStyle(0)
        htitle.SetTextFont(43)
        htitle.SetTextSize(14)
        htitle.SetTextAlign(13)
        htitle.AddText(binTitle)
        htitle.AddText("|#eta_{jet}| < 0.5")
        htitle.AddText("#it{p}_{T,ch jet} > 6 GeV/#it{c}")
        htitle.Draw()
        globalList.append(htitle)
        invMassFitter.Draw("same");

    c.cd(1)    
    paveALICE = ROOT.TPaveText(0.19, 0.80, 0.59, 0.89, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(14)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Simulation")
    paveALICE.AddText("PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.Draw()

    leg1 = ROOT.TLegend(0.20, 0.47, 0.58, 0.74, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(14)
    leg1.SetTextAlign(13)
    leg1.AddEntry(invMassHisto_copy, "Inv.Mass", "pe")
    leg1.AddEntry(h1, "S-B Window", "f")
    leg1.AddEntry(hsig, "Signal Window", "f")
    leg1.AddEntry(invMassHistoLS_copy, "L-S Inv.Mass", "l")
    leg1.AddEntry(invMassFitter.GetFitFunction(), "Fit Sig+Bkg", "l")
    leg1.AddEntry(invMassFitter.GetBkgFunction(), "Fit Bkg-only", "l")
    leg1.Draw()

def PlotInvMassSideBands(sideBandWindowHisto, signalWindowHisto):
    hsb1 = sideBandWindowHisto.DrawCopy("hist same")
    hsb1.GetXaxis().SetRange(9, 16)
    hsb1.SetFillColor(ROOT.kRed-10)
    hsb1.SetFillStyle(1001)
    hsb1.SetLineColor(ROOT.kRed-10)
    hsb2 = sideBandWindowHisto.DrawCopy("hist same")
    hsb2.GetXaxis().SetRange(34, 41)
    hsb2.SetFillColor(ROOT.kRed-10)
    hsb2.SetFillStyle(1001)
    hsb2.SetLineColor(ROOT.kRed-10)
    hsig = signalWindowHisto.DrawCopy("hist same")
    hsig.GetXaxis().SetRange(22, 29)
    hsig.SetFillColor(ROOT.kBlue-10)
    hsig.SetFillStyle(1001)
    hsig.SetLineColor(ROOT.kBlue-10)
    globalList.append(hsb1)
    globalList.append(hsb2)
    globalList.append(hsig)

    return hsb1, hsb2, hsig

def OpenFile(config, suffix=""):
    ROOT.TH1.AddDirectory(0)
    if suffix:
        file = ROOT.TFile.Open("{0}/{1}/{2}_{3}.root".format(config["input_path"], config["train"], config["name"], suffix))
    else:
        file = ROOT.TFile.Open("{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"]))
    return file

def LoadHistograms(lname, file):
    rlist = file.Get(lname)
    return ExtractRootList(rlist)

def ExtractRootList(input):
    result = dict()
    for h in input:
        if isinstance(h, ROOT.TH1) or isinstance(h, ROOT.MassFitter):
            result[h.GetName()] = h
        elif isinstance(h, ROOT.TCollection):
            result[h.GetName()] = ExtractRootList(h)
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Side Band analysis.')
    parser.add_argument('actions', metavar='action',
                        help='Actions to be taken', nargs='*')
    parser.add_argument('-o', metavar='path',
                        help='Output path', default='../notes/HQ16/img')
    parser.add_argument('-f', metavar='format',
                        help='Format (pdf, eps, png,...)', default='eps')
    args = parser.parse_args()

    main(args.actions, args.o, args.f)
    
    IPython.embed()