#!/usr/bin/env python
#python script to generate plots requested for approval to be presented at HQ16 conference

import argparse
import yaml
import IPython
import ROOT
import subprocess

globalList = []
canvases = []

def main(config, actions):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    file = OpenFile(config)

    if "all" in actions or "SB" in actions:
        SideBandAnalysisPlots(file, config)

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
    spectrumConfig = None
    for s in config["analysis"][0]["spectra"]:
        if s["name"] == "D_Tagged_Jet_PtD_20_Spectrum_SideBand":
            spectrumConfig = s
            break
    minSBsigma = spectrumConfig["side_band"]["min_sigmas"]
    maxSBsigma = spectrumConfig["side_band"]["max_sigmas"]
    sigSigma = spectrumConfig["side_band"]["max_signal_sigmas"]

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

        #h.GetXaxis().SetTitle("#it{m}(K#pi) GeV/#it{c}^{2}")
        #h.GetYaxis().SetTitle("counts")
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

    paveFit = ROOT.TPaveText(0.50, 0.53, 0.99, 0.68, "NB NDC")
    globalList.append(paveFit)
    paveFit.SetBorderSize(0)
    paveFit.SetFillStyle(0)
    paveFit.SetTextFont(43)
    paveFit.SetTextSize(14)
    paveFit.SetTextAlign(23)

    paveFit.AddText(invMassFitter.GetSignalMeanString().Data())
    paveFit.AddText(invMassFitter.GetSignalWidthString().Data())
    #paveFit.AddText(invMassFitter.GetBkgPar1String().Data())
    paveFit.AddText("N={0}".format(int(invMassFitter.GetTotalEntries())))
    paveFit.Draw()

def OpenFile(config):
    ROOT.TH1.AddDirectory(0)
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
    parser.add_argument('--yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument('actions', metavar='action',
                        help='Actions to be taken', nargs='*')
    args = parser.parse_args()
    
    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.actions)
    
    IPython.embed()