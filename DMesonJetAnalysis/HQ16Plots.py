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

    if "all" in actions or "SB" in actions:
        SideBandAnalysisPlots(file, config)

    if "all" in actions or "LS" in actions:
        LikeSignAnalysisPlots(file, config)

    if "all" in actions or "compare" in actions:
        CompareMethods(file, config)

    if "all" in actions or "unc" in actions:
        CompareUncertainties(file, fileW, config)

    if "all" in actions or "eff" in actions:
        EfficiencyPlots(fileRM, configRM)

    if "all" in actions or "resp" in actions:
        DetectorResponsePlots(fileRM, configRM)

    if "all" in actions or "stat" in actions:
        StatisticalUncertaintyData(fileData, configData)

    for c in canvases:
        c.SaveAs("{0}/{1}.{2}".format(output_path, c.GetName(), output_type))

def StatisticalUncertaintyData(file, config):
    dataList = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum", file)
    hist = dataList["D0_D_Tagged_Jet_PtD_20_Spectrum_Unc"]
    cname = "HQ16_WorkInProgress_StatisticalUncertainty"
    c = ROOT.TCanvas(cname, cname, 650, 500)
    canvases.append(c)
    c.cd()
    h = hist.DrawCopy("hist")

    h.SetLineColor(ROOT.kBlue+2)
    h.SetLineWidth(2)
    h.SetFillColorAlpha(ROOT.kBlue+2, 0.25)
    h.SetFillStyle(1001)
    h.GetYaxis().SetRangeUser(0, 0.25)
    h.GetYaxis().SetTitleOffset(1.2)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetXaxis().SetTitle("#it{p}_{T,jet}^{ch,reco} (GeV/#it{c})")

    pave = ROOT.TPaveText(0.12, 0.88, 0.4, 0.55, "NB NDC")
    pave.SetFillStyle(0)
    pave.SetBorderSize(0)
    pave.SetTextFont(43)
    pave.SetTextSize(15)
    pave.AddText("ALICE Work In Progress")
    pave.AddText("pp, #sqrt{#it{s}} = 7 TeV")
    pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4")
    pave.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
    pave.AddText("#it{p}_{T,D} > 2 GeV/#it{c}")
    pave.Draw()

    globalList.append(c)
    globalList.append(h)
    globalList.append(pave)

def DetectorResponsePlots(file, config):
    respList = LoadHistograms("D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum", file)
    resp = respList["DetectorResponse"]
    histList = [resp["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_5_6"], resp["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_8_10"], 
                resp["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_18_24"]]
    histList[0].SetTitle("5 < #it{p}_{T,jet}^{ch,truth} < 6 GeV/#it{c}")
    histList[1].SetTitle("8 < #it{p}_{T,jet}^{ch,truth} < 10 GeV/#it{c}")
    histList[2].SetTitle("18 < #it{p}_{T,jet}^{ch,truth} < 24 GeV/#it{c}")
    blank = PlotMultiHistogram(histList, "HQ16_Simulation_DetectorResponse", "#(){#it{p}_{T,jet}^{ch,reco}-#it{p}_{T,jet}^{ch,truth}} / #it{p}_{T,jet}^{ch,truth}", "Probability density",
                               [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kGreen+1], [ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenTriangleUp])
    blank.GetYaxis().SetRangeUser(0, 18)
    blank.GetXaxis().SetRangeUser(-0.8, 0.8)
    
    histList = [respList["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_EnergyScaleShift"], respList["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_EnergyScaleShiftMedian"]]
    histList[0].SetTitle("Mean")
    histList[1].SetTitle("Median")
    blank = PlotMultiHistogram(histList, "HQ16_Simulation_EnergyScaleShift", "#it{p}_{T,jet}^{ch,truth} (GeV/#it{c})", "#(){#it{p}_{T,jet}^{ch,reco}-#it{p}_{T,jet}^{ch,truth}} / #it{p}_{T,jet}^{ch,truth}")
    blank.GetYaxis().SetRangeUser(-0.07, 0.08)

    histList = [respList["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_Spectrum_DetectorResponse_Resolution"]]
    histList[0].SetTitle("Resolution")
    blank = PlotMultiHistogram(histList, "HQ16_Simulation_Resolution", "#it{p}_{T,jet}^{ch,truth} (GeV/#it{c})", "#sigma#(){#it{p}_{T,jet}^{ch,reco}-#it{p}_{T,jet}^{ch,truth}} / #it{p}_{T,jet}^{ch,truth}")
    blank.GetYaxis().SetRangeUser(0.05, 0.18)

def EfficiencyPlots(file, config):
    etaEff = LoadHistograms("D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_JetEta", file)
    histList = [etaEff["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_JetEta_Efficiency_JetPt_500_2400"], etaEff["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_JetEta_Efficiency_JetPt_500_800"], 
                etaEff["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_JetEta_Efficiency_JetPt_800_1300"], etaEff["D0_Jet_AKTChargedR040_pt_scheme_D_Tagged_Jet_JetEta_Efficiency_JetPt_1300_2400"]]
    histList[0].SetTitle("5 < #it{p}_{T,jet}^{ch} < 24 GeV/#it{c}")
    histList[1].SetTitle("5 < #it{p}_{T,jet}^{ch} < 8 GeV/#it{c}")
    histList[2].SetTitle("8 < #it{p}_{T,jet}^{ch} < 13 GeV/#it{c}")
    histList[3].SetTitle("13 < #it{p}_{T,jet}^{ch} < 24 GeV/#it{c}")
    PlotMultiHistogram(histList, "HQ16_Simulation_EfficiencyVsJetEta", "#eta_{jet}", "Efficiency")
    
    ptEff = LoadHistograms("D0_Jet_AKTChargedR040_pt_scheme_D_Spectra", file)
    histList = [ptEff["D0_Jet_AKTChargedR040_pt_scheme_D_Spectra_Efficiency_JetPt_500_2400"], ptEff["D0_Jet_AKTChargedR040_pt_scheme_D_Spectra_Efficiency_JetPt_500_800"], 
                ptEff["D0_Jet_AKTChargedR040_pt_scheme_D_Spectra_Efficiency_JetPt_800_1300"], ptEff["D0_Jet_AKTChargedR040_pt_scheme_D_Spectra_Efficiency_JetPt_1300_2400"]]
    histList[0].SetTitle("5 < #it{p}_{T,jet}^{ch} < 24 GeV/#it{c}")
    histList[1].SetTitle("5 < #it{p}_{T,jet}^{ch} < 8 GeV/#it{c}")
    histList[2].SetTitle("8 < #it{p}_{T,jet}^{ch} < 13 GeV/#it{c}")
    histList[3].SetTitle("13 < #it{p}_{T,jet}^{ch} < 24 GeV/#it{c}")
    PlotMultiHistogram(histList, "HQ16_Simulation_EfficiencyVsDPt", "#it{p}_{T,D} (GeV/#it{c})", "Efficiency")

def PlotMultiHistogram(histList, cname, xaxisTitle, yaxisTitle, colors=None, markers=None):
    c = ROOT.TCanvas(cname, cname)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.08)
    c.SetRightMargin(0.08)
    c.cd()
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
        colors = [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kMagenta+2, ROOT.kAzure+2, ROOT.kPink+2]
    if not markers:
        markers = [ROOT.kStar, ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond, ROOT.kFullStar, ROOT.kFullCross]
    max = 0;
    leg = ROOT.TLegend(0.62, 0.90-len(histList)*0.055, 0.92, 0.90)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(16)
    for color,marker,eff in zip(colors,markers,histList):
        h = eff.Clone()
        globalList.append(h)
        h.SetMarkerStyle(marker)
        h.SetMarkerSize(0.9)
        h.SetMarkerColor(color)
        h.SetLineColor(color)
        leg.AddEntry(h, h.GetTitle(), "pe")
        h.Draw("same")
        for i in range(1, h.GetNbinsX()+1):
            y = h.GetBinContent(i)
            if y > max:
                max = y

    if len(histList) > 1:
        leg.Draw()
    globalList.append(leg)
    blank.SetMaximum(max*1.8)

    paveALICE = ROOT.TPaveText(0.13, 0.70, 0.52, 0.90, "NB NDC")
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
    return blank

def CompareUncertainties(file, fileW, config):
    LSlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign", file)
    SBlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand", file)
    IMlist = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum", file)
    IMlist_eff = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum", fileW)

    unc = [IMlist["D0_D_Tagged_Jet_PtD_20_Spectrum_Unc"], SBlist["D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand_Unc"], LSlist["D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign_Unc"]]

    cname = "HQ16_Simulation_UncertaintyComparison"
    c = ROOT.TCanvas(cname, cname, 650, 500)
    canvases.append(c)
    globalList.append(c)
    
    h = IMlist_eff["D0_D_Tagged_Jet_PtD_20_Spectrum_Unc"].DrawCopy("hist")
    globalList.append(h)
    h.SetFillColorAlpha(ROOT.kYellow+4,1.00)
    h.SetLineColorAlpha(ROOT.kYellow+4,1.00)
    #h.SetFillStyle(3002)
    h.SetFillStyle(0)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(17)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(16)
    h.GetYaxis().SetTitleOffset(1.3)
    h.GetYaxis().SetRangeUser(0, 0.45)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetXaxis().SetTitle("#it{p}_{T,jet}^{ch,reco} GeV/#it{c}")
    
    leg1 = ROOT.TLegend(0.11, 0.45, 0.55, 0.66, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(17)
    leg1.SetTextAlign(13)
    leg1.AddEntry(h, "Inv.Mass Fit w/ eff.", "l")

    colors = [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2]
    fillStyles = [0]*3#[3205, 3245, 3254]
    gradients = [1]*3
    titles = ["Inv.Mass Fit", "Side-band Method", "Like-sign Method"]

    for s,color,fillStyle,title,grad in zip(unc,colors,fillStyles,titles,gradients):
        h = s.DrawCopy("same hist")
        globalList.append(h)
        h.SetFillColorAlpha(color,grad)
        h.SetLineColorAlpha(color,1.00)
        h.SetFillStyle(fillStyle)
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

    cname = "HQ16_Simulation_MethodComparison"
    c = ROOT.TCanvas(cname, cname, 700, 700)
    canvases.append(c)
    c.Divide(1, 2)
    globalList.append(c)
    padMain = c.cd(1)
    padMain.SetPad(0, 0.35, 1, 1)
    padMain.SetBottomMargin(0)
    hS = Tspectrum.DrawCopy("hist")
    globalList.append(hS)
    hS.SetLineColorAlpha(ROOT.kBlack,0.25)
    hS.SetFillColorAlpha(ROOT.kBlack,0.25)
    hS.SetFillStyle(1001)
    hS.GetYaxis().SetRangeUser(0.001,1100)
    hS.GetYaxis().SetTitleFont(43)
    hS.GetYaxis().SetTitleSize(17)
    hS.GetYaxis().SetLabelFont(43)
    hS.GetYaxis().SetLabelSize(16)
    hS.GetYaxis().SetTitleOffset(1.6)

    padRatio = c.cd(2)
    padRatio.SetPad(0, 0., 1, 0.35)
    padRatio.SetTopMargin(0)
    padRatio.SetBottomMargin(0.25)
    hSratio = Tspectrum.DrawCopy("hist")
    globalList.append(hSratio)
    hSratio.Divide(Tspectrum)
    hSratio.SetFillStyle(0)
    hSratio.SetLineColor(ROOT.kBlack)
    hSratio.SetLineStyle(2)
    hSratio.SetLineWidth(2)
    hSratio.GetXaxis().SetTitleFont(43)
    hSratio.GetXaxis().SetTitleSize(17)
    hSratio.GetXaxis().SetLabelFont(43)
    hSratio.GetXaxis().SetLabelSize(16)
    hSratio.GetYaxis().SetTitleFont(43)
    hSratio.GetYaxis().SetTitleSize(17)
    hSratio.GetYaxis().SetLabelFont(43)
    hSratio.GetYaxis().SetLabelSize(16)
    hSratio.GetYaxis().SetTitleOffset(1.6)
    hSratio.GetXaxis().SetTitleOffset(3.4)
    hSratio.GetXaxis().SetTitle("#it{p}_{T,jet}^{ch,reco} GeV/#it{c}")
    hSratio.GetYaxis().SetTitle("ratio")
    hSratio.GetYaxis().SetNdivisions(504)

    colors = [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2]
    markers = [ROOT.kFullSquare, ROOT.kFullCircle, ROOT.kFullTriangleUp]
    sizes = [1.2, 1.1, 1.2]
    fillStyles = [3002, 3245, 3254]
    titles = ["Inv.Mass Fit", "Side-band Method", "Like-sign Method"]
    
    leg1 = ROOT.TLegend(0.52, 0.28, 0.96, 0.49, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(18)
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
        
        leg1.AddEntry(h, title, "pe")
        
        padRatio.cd()
        h = s.DrawCopy("same E2")
        globalList.append(h)
        h.Divide(Tspectrum)
        h.SetMarkerColor(color)
        h.SetMarkerStyle(marker)
        h.SetMarkerSize(size)
        h.SetLineColor(color)
        h.SetFillColorAlpha(color,0.55)
        h.SetFillStyle(fillStyle)

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

def LikeSignAnalysisPlots(file, config):
    spectrum = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_LikeSign", file)
    invMass = LoadHistograms("D0_D_Tagged_Jets_PtD_20", file)
    invMassLS = LoadHistograms("2ProngLikeSign_D_Tagged_Jets_PtD_20", file)
    PlotLSInvMass(invMass, invMassLS, spectrum, config)

def PlotLSInvMass(invMassPlotList, invMassLSPlotList, spectrumPlotList, config):
    bins = ["500_600", "1400_1800"]
    binTitles = ["5 < #it{p}_{T,jet}^{ch,reco} < 6 GeV/#it{c}", "14 < #it{p}_{T,jet}^{ch,reco} < 18 GeV/#it{c}"]
    cname = "HQ16_Simulation_InvMassLS"
    c = ROOT.TCanvas(cname, cname, len(bins)*400, 400)
    canvases.append(c)
    c.Divide(len(bins), 1)
    globalList.append(c)
    for i,(bin,binTitle) in enumerate(zip(bins, binTitles)):
        invMassHisto = invMassPlotList["InvMass_D0_JetPt_{0}".format(bin)]
        invMassLSHisto = invMassLSPlotList["InvMass_2ProngLikeSign_JetPt_{0}".format(bin)]
        invMassFitter = invMassPlotList["InvMass_D0_JetPt_{0}_fitter".format(bin)]

        pad = c.cd(i+1)
        pad.SetLeftMargin(0.17)
        pad.SetRightMargin(0.05)
        pad.SetTopMargin(0.09)
        pad.SetBottomMargin(0.13)

        hLS = invMassLSHisto.DrawCopy("hist")
        hLS.SetFillStyle(1001)
        hLS.SetFillColorAlpha(ROOT.kGreen+2, 0.4)
        hLS.SetLineColorAlpha(ROOT.kGreen+2, 0.4)
        globalList.append(hLS)

        h = invMassHisto.DrawCopy("same")
        globalList.append(h)

        hLS.SetMaximum(h.GetMaximum()*1.8)
        hLS.GetXaxis().SetTitleFont(43)
        hLS.GetXaxis().SetTitleOffset(1.1)
        hLS.GetXaxis().SetTitleSize(19)
        hLS.GetXaxis().SetLabelFont(43)
        hLS.GetXaxis().SetLabelOffset(0.009)
        hLS.GetXaxis().SetLabelSize(18)
        hLS.GetYaxis().SetTitleFont(43)
        hLS.GetYaxis().SetTitleOffset(1.5)
        hLS.GetYaxis().SetTitleSize(19)
        hLS.GetYaxis().SetLabelFont(43)
        hLS.GetYaxis().SetLabelOffset(0.009)
        hLS.GetYaxis().SetLabelSize(18)

        htitle = ROOT.TPaveText(0.12, 0.99, 0.95, 0.93, "NB NDC")
        htitle.SetBorderSize(0)
        htitle.SetFillStyle(0)
        htitle.SetTextFont(43)
        htitle.SetTextSize(18)
        htitle.AddText(binTitle)
        htitle.Draw()
        globalList.append(htitle)
        DrawFitResults(invMassFitter)

    c.cd(1)    
    paveALICE = ROOT.TPaveText(0.17, 0.81, 0.96, 0.92, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(14)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Simulation")
    paveALICE.AddText("PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.Draw()

    leg1 = ROOT.TLegend(0.19, 0.42, 0.57, 0.62, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(14)
    leg1.SetTextAlign(13)
    leg1.AddEntry(h, "U-S candidates", "pe")
    leg1.AddEntry(hLS, "L-S candidates", "f")
    leg1.AddEntry(invMassFitter.GetFitFunction(), "Fit Sig+Bkg", "l")
    leg1.AddEntry(invMassFitter.GetBkgFunction(), "Fit Bkg-only", "l")
    leg1.Draw()
    
    c.cd(2)    
    paveALICE = ROOT.TPaveText(0.17, 0.81, 0.96, 0.92, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(14)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4, |#eta_{jet}| < 0.5")
    paveALICE.AddText("2 < #it{p}_{T,D} < 24 GeV/#it{c}")
    paveALICE.Draw()

def SideBandAnalysisPlots(file, config):
    spectrum = LoadHistograms("D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand", file)
    invMass = LoadHistograms("D0_D_SideBandAnalysis_PtD_20", file)
    PlotSBInvMass(invMass, spectrum, config)
    PlotSBSpectra(spectrum, config)

def PlotSBSpectra(spectrumList, config):
    cname = "HQ16_Simulation_SpectraSB"
    c = ROOT.TCanvas(cname, cname, 650, 500)
    globalList.append(c)
    canvases.append(c)
    c.SetLogy()
    c.cd()
    SBWinSpectrum = spectrumList["SideBandAnalysis"]["D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand_SideBandWindowTotal"]
    SigWinSpectrum = spectrumList["SideBandAnalysis"]["D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand_SignalWindowTotal"]
    spectrum = spectrumList["D0_D_Tagged_Jet_PtD_20_Spectrum_SideBand"]

    hSigWin = SigWinSpectrum.DrawCopy()
    globalList.append(hSigWin)
    hSigWin.SetMarkerStyle(ROOT.kFullTriangleUp)
    hSigWin.SetMarkerSize(0.9)
    hSigWin.SetMarkerColor(ROOT.kRed+2)
    hSigWin.SetLineColor(ROOT.kRed+2)
    hSigWin.GetYaxis().SetRangeUser(5e1, 1e4)
    hSigWin.GetXaxis().SetTitleOffset(1.2)
    hSigWin.GetXaxis().SetTitle("#it{p}_{T,jet}^{ch,reco} GeV/#it{c}")

    hSBWin = SBWinSpectrum.DrawCopy("same")
    globalList.append(hSBWin)
    hSBWin.SetMarkerStyle(ROOT.kFullTriangleDown)
    hSBWin.SetMarkerSize(0.9)
    hSBWin.SetMarkerColor(ROOT.kGreen+2)
    hSBWin.SetLineColor(ROOT.kGreen+2)

    h = spectrum.DrawCopy("same")
    globalList.append(h)
    h.SetMarkerStyle(ROOT.kFullCircle)
    h.SetMarkerSize(0.9)
    h.SetMarkerColor(ROOT.kBlue+2)
    h.SetLineColor(ROOT.kBlue+2)
    
    paveALICE = ROOT.TPaveText(0.09, 0.73, 0.88, 0.90, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(17)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Simulation")
    paveALICE.AddText("PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4, |#eta_{jet}| < 0.5 with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
    paveALICE.AddText("2 < #it{p}_{T,D} < 24 GeV/#it{c}")
    paveALICE.Draw()
    
    leg = ROOT.TLegend(0.09, 0.12, 0.46, 0.38, "", "NB NDC")
    globalList.append(leg)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(16)
    leg.SetTextAlign(13)
    leg.AddEntry(hSigWin, "Signal candidates", "pe")
    leg.AddEntry(hSigWin, "#||{#it{m}(#piK) - #it{m}_{fit}} < 2#sigma_{fit}", "")
    leg.AddEntry(hSBWin, "SB candidates", "pe")
    leg.AddEntry(hSBWin, "4#sigma_{fit} < #||{#it{m}(#piK) - #it{m}_{fit}} < 8#sigma_{fit}", "")
    leg.AddEntry(h, "Signal - SB candidantes", "pe")
    leg.Draw()

def PlotSBInvMass(invMassPlotList, spectrumPlotList, config):
    bins = ["200_300", "600_700", "1200_1600"]
    binTitles = ["2 < #it{p}_{T,D} < 3 GeV/#it{c}", "6 < #it{p}_{T,D} < 7 GeV/#it{c}", "12 < #it{p}_{T,D} < 16 GeV/#it{c}"]
    cname = "HQ16_Simulation_InvMassSB"
    c = ROOT.TCanvas(cname, cname, len(bins)*400, 400)
    canvases.append(c)
    c.Divide(len(bins), 1)
    globalList.append(c)
    for i,(bin,binTitle) in enumerate(zip(bins, binTitles)):
        invMassHisto = invMassPlotList["InvMass_D0_DPt_{0}".format(bin)]
        invMassFitter = invMassPlotList["InvMass_D0_DPt_{0}_fitter".format(bin)]
        invMassHistoSB = spectrumPlotList["SideBandAnalysis"]["InvMassSBWindow_D0_DPt_{0}".format(bin)]
        invMassHistoSig = spectrumPlotList["SideBandAnalysis"]["InvMassSigWindow_D0_DPt_{0}".format(bin)]

        pad = c.cd(i+1)
        pad.SetLeftMargin(0.17)
        pad.SetRightMargin(0.05)
        pad.SetTopMargin(0.09)
        pad.SetBottomMargin(0.13)

        (h, hsig) = PlotInvMassSideBands(invMassHisto, invMassFitter, invMassHistoSB, invMassHistoSig)
        globalList.append(h)
        invMassHisto_copy = invMassHisto.DrawCopy("same")
        globalList.append(invMassHisto_copy)

        h.SetMaximum(invMassHisto_copy.GetMaximum()*1.8)
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
        htitle = ROOT.TPaveText(0.12, 0.99, 0.95, 0.93, "NB NDC")
        htitle.SetBorderSize(0)
        htitle.SetFillStyle(0)
        htitle.SetTextFont(43)
        htitle.SetTextSize(18)
        htitle.AddText(binTitle)
        htitle.Draw()
        globalList.append(htitle)
        DrawFitResults(invMassFitter)

    c.cd(1)    
    paveALICE = ROOT.TPaveText(0.17, 0.81, 0.96, 0.92, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(14)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Simulation")
    paveALICE.AddText("PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.Draw()

    c.cd(2)
    leg1 = ROOT.TLegend(0.19, 0.39, 0.57, 0.55, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(14)
    leg1.SetTextAlign(13)
    leg1.AddEntry(invMassHisto_copy, "Inv.Mass", "pe")
    leg1.AddEntry(h, "SB Window", "f")
    leg1.AddEntry(hsig, "Signal Window", "f")
    leg1.Draw()

    c.cd(3)
    leg2 = ROOT.TLegend(0.19, 0.45, 0.53, 0.58, "", "NB NDC")
    globalList.append(leg2)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextFont(43)
    leg2.SetTextSize(14)
    leg2.SetTextAlign(13)
    leg2.AddEntry(invMassFitter.GetFitFunction(), "Fit Sig+Bkg", "l")
    leg2.AddEntry(invMassFitter.GetBkgFunction(), "Fit Bkg-only", "l")
    leg2.Draw()

def PlotInvMassSideBands(invMassHisto, invMassFitter, sideBandWindowHisto, signalWindowHisto):
    hsb = sideBandWindowHisto.DrawCopy("hist")
    hsb.SetFillColorAlpha(ROOT.kGreen+2, 0.4)
    hsb.SetFillStyle(1001)
    hsb.SetLineColorAlpha(ROOT.kGreen+2, 0.4)
    hsig = signalWindowHisto.DrawCopy("hist same")
    hsig.SetFillColorAlpha(ROOT.kRed+2, 0.4)
    hsig.SetFillStyle(1001)
    hsig.SetLineColorAlpha(ROOT.kRed+2, 0.4)
    globalList.append(hsb)
    globalList.append(hsig)

    return hsb,hsig

def DrawFitResults(invMassFitter):
    invMassFitter.Draw("same");

    fitStatus = int(invMassFitter.GetFitStatus())
    if fitStatus == 0:
        chi2Text = invMassFitter.GetChisquareString().Data()
    else:
        chi2Text = "Fit failed"

    paveSig = ROOT.TPaveText(0.18, 0.67, 0.54, 0.80, "NB NDC")
    globalList.append(paveSig)
    paveSig.SetBorderSize(0)
    paveSig.SetFillStyle(0)
    paveSig.SetTextFont(43)
    paveSig.SetTextSize(14)
    paveSig.SetTextAlign(13)
    paveSig.AddText("{0}, {1}".format(invMassFitter.GetSignalString().Data(), 
                                      invMassFitter.GetBackgroundString().Data()))
    paveSig.AddText("{0}, {1}".format(invMassFitter.GetSignalOverSqrtSignalBackgroundString().Data(),
                                      chi2Text))
    paveSig.Draw()

    paveFit = ROOT.TPaveText(0.47, 0.50, 0.96, 0.65, "NB NDC")
    globalList.append(paveFit)
    paveFit.SetBorderSize(0)
    paveFit.SetFillStyle(0)
    paveFit.SetTextFont(43)
    paveFit.SetTextSize(14)
    paveFit.SetTextAlign(33)

    paveFit.AddText(invMassFitter.GetSignalMeanString().Data())
    paveFit.AddText(invMassFitter.GetSignalWidthString().Data())
    paveFit.AddText("Entries={0}".format(int(invMassFitter.GetTotalEntries())))
    paveFit.Draw()

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
                        help='Format (pdf, eps, png,...)', default='pdf')
    args = parser.parse_args()

    main(args.actions, args.o, args.f)
    
    IPython.embed()