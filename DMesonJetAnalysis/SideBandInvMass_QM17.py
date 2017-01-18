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

def PlotSBInvMass(pad, ptmin, ptmax, sbList, dptbinList, plotleg=False):
    pad.SetTicks(1, 1)
    pad.SetLeftMargin(0.16)
    pad.SetRightMargin(0.02)
    pad.SetTopMargin(0.09)
    pad.SetBottomMargin(0.13)
    sbHist = sbList.FindObject("InvMassSBWindow_AnyINT_D0_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
    h = sbHist.DrawCopy("axis")

    minbin_l = 0
    maxbin_l = 0
    minbin_r = 0
    maxbin_r = 0
    for ibin in range(1, sbHist.GetNbinsX() + 1):
        if sbHist.GetBinContent(ibin) > 0 and minbin_l == 0:
            minbin_l = ibin
            continue
        if sbHist.GetBinContent(ibin) == 0 and maxbin_l == 0 and not minbin_l == 0:
            maxbin_l = ibin - 1
            continue
        if sbHist.GetBinContent(ibin) > 0 and minbin_r == 0 and not maxbin_l == 0:
            minbin_r = ibin
            continue
        if sbHist.GetBinContent(ibin) == 0 and maxbin_r == 0 and not minbin_r == 0:
            maxbin_r = ibin - 1
            break

    sbHist_copy_l = sbHist.DrawCopy("hist same")
    globalList.append(sbHist_copy_l)
    sbHist_copy_l.SetFillColor(ROOT.kGreen - 6)
    sbHist_copy_l.SetLineColor(ROOT.kGreen - 6)
    sbHist_copy_l.GetXaxis().SetRange(minbin_l, maxbin_l)

    sbHist_copy_r = sbHist.DrawCopy("hist same")
    globalList.append(sbHist_copy_r)
    sbHist_copy_r.SetFillColor(ROOT.kGreen - 6)
    sbHist_copy_r.SetLineColor(ROOT.kGreen - 6)
    sbHist_copy_r.GetXaxis().SetRange(minbin_r, maxbin_r)

    sigHist = sbList.FindObject("InvMassSigWindow_AnyINT_D0_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
    minbin_s = 0
    maxbin_s = 0
    for ibin in range(1, sigHist.GetNbinsX() + 1):
        if sigHist.GetBinContent(ibin) > 0 and minbin_s == 0:
            minbin_s = ibin
            continue
        if sigHist.GetBinContent(ibin) == 0 and maxbin_s == 0 and not minbin_s == 0:
            maxbin_s = ibin - 1
            break
    sigHist_copy = sigHist.DrawCopy("hist same")
    globalList.append(sigHist_copy)
    sigHist_copy.SetFillColor(ROOT.kRed - 6)
    sigHist_copy.SetLineColor(ROOT.kRed - 6)
    sigHist_copy.GetXaxis().SetRange(minbin_s, maxbin_s)

    invMassHist = dptbinList.FindObject("InvMass_AnyINT_D0_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
    invMassHist_copy = invMassHist.DrawCopy("p0 x0 same")
    globalList.append(invMassHist_copy)
    invMassHist_copy.SetLineColor(ROOT.kBlue + 2)
    invMassHist_copy.SetMarkerColor(ROOT.kBlue + 2)
    invMassHist_copy.SetMarkerStyle(ROOT.kFullCircle)
    invMassHist_copy.SetMarkerSize(1.0)

    fitter = dptbinList.FindObject("InvMass_AnyINT_D0_DPt_{0:.0f}_{1:.0f}_fitter".format(ptmin * 100, ptmax * 100))
    globalList.append(fitter)
    fitter.Draw("same")

    diff = invMassHist_copy.GetMaximum() - invMassHist_copy.GetMinimum()
    miny = invMassHist_copy.GetMinimum() - 0.2 * diff
    if miny < 0: miny = 0
    maxy = invMassHist_copy.GetMaximum() + 0.8 * diff
    h.SetMaximum(maxy)
    h.SetMinimum(miny)

    binTitle = "{0:.0f} < #it{{p}}_{{T,D}} < {1:.0f} GeV/#it{{c}}".format(ptmin, ptmax)

    h.GetYaxis().SetTitle("counts / efficiency")
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleOffset(2.2)
    h.GetXaxis().SetTitleSize(19)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelOffset(0.009)
    h.GetXaxis().SetLabelSize(18)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleOffset(3.2)
    h.GetYaxis().SetTitleSize(19)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelOffset(0.009)
    h.GetYaxis().SetLabelSize(23)
    htitle = ROOT.TPaveText(0.15, 0.90, 0.95, 0.99, "NB NDC")
    htitle.SetBorderSize(0)
    htitle.SetFillStyle(0)
    htitle.SetTextFont(43)
    htitle.SetTextSize(20)
    htitle.SetTextAlign(21)
    htitle.AddText(binTitle)
    htitle.Draw()
    globalList.append(htitle)

    if plotleg:
        leg = ROOT.TLegend(0.18, 0.58, 0.57, 0.87, "", "NB NDC")
        globalList.append(leg)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(20)
        leg.SetTextAlign(13)
        leg.AddEntry(invMassHist_copy, "D^{0}-jet Candidates", "pe")
        leg.AddEntry(fitter.GetFitFunction(), "Fit Sig+Bkg", "l")
        leg.AddEntry(fitter.GetBkgFunction(), "Fit Bkg-only", "l")
        leg.AddEntry(sigHist_copy, "Signal Window", "f")
        leg.AddEntry(sbHist_copy_l, "S-B Window", "f")
        leg.Draw()

def PlotSBSpectra(pad, ptmin, ptmax, sbList):
    pad.SetTicks(1, 1)
    pad.SetLogy()
    pad.SetLeftMargin(0.16)
    pad.SetRightMargin(0.02)
    pad.SetTopMargin(0.04)
    pad.SetBottomMargin(0.13)

    sbHist = sbList.FindObject("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_SideBandWindow_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))
    h = sbHist.DrawCopy("axis")
    globalList.append(h)

    sbHist_copy = sbHist.DrawCopy("p0 same")
    globalList.append(sbHist_copy)
    sbHist_copy.SetMarkerColor(ROOT.kGreen + 2)
    sbHist_copy.SetLineColor(ROOT.kGreen + 2)
    sbHist_copy.SetMarkerStyle(ROOT.kOpenCircle)
    sbHist_copy.SetMarkerSize(1.3)

    sigHist = sbList.FindObject("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_SignalWindow_DPt_{0:.0f}_{1:.0f}".format(ptmin * 100, ptmax * 100))

    sigHist_copy = sigHist.DrawCopy("p0 same")
    globalList.append(sigHist_copy)
    sigHist_copy.SetMarkerColor(ROOT.kRed + 2)
    sigHist_copy.SetLineColor(ROOT.kRed + 2)
    sigHist_copy.SetMarkerStyle(ROOT.kOpenSquare)
    sigHist_copy.SetMarkerSize(1.3)

    subHist_copy = sigHist.DrawCopy("p0 same")
    subHist_copy.Add(sbHist_copy, -1)
    globalList.append(subHist_copy)
    subHist_copy.SetMarkerColor(ROOT.kBlue + 2)
    subHist_copy.SetLineColor(ROOT.kBlue + 2)
    subHist_copy.SetMarkerStyle(ROOT.kOpenDiamond)
    subHist_copy.SetMarkerSize(1.7)

    h.GetYaxis().SetTitle("counts / efficiency")
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleOffset(2.2)
    h.GetXaxis().SetTitleSize(19)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelOffset(0.009)
    h.GetXaxis().SetLabelSize(18)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleOffset(3.2)
    h.GetYaxis().SetTitleSize(19)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelOffset(0.009)
    h.GetYaxis().SetLabelSize(23)

    miny = min([DMesonJetUtils.FindMinimum(sbHist_copy), DMesonJetUtils.FindMinimum(sigHist_copy), DMesonJetUtils.FindMinimum(subHist_copy)])
    maxy = max([DMesonJetUtils.FindMaximum(sbHist_copy), DMesonJetUtils.FindMaximum(sigHist_copy), DMesonJetUtils.FindMaximum(subHist_copy)])
    miny /= 2
    maxy *= 3
    h.SetMaximum(maxy)
    h.SetMinimum(miny)

def SideBandPlot():
    loader = RawYieldSpectrumLoader.RawYieldSpectrumLoader(input_path, "Jets_EMC_pp_823_824_825_826", "LHC10_Train823_efficiency")
    loader.fDMeson = "D0"
    loader.fJetType = "Charged"
    loader.fJetRadius = "R040"
    loader.fSpectrumName = "JetPtSpectrum_DPt_30_SideBand"
    loader.LoadDataListFromDMesonJetAnalysis()
    dptbinList = loader.fDataJetList.FindObject("D0_Charged_R040_DPtBins_JetPt_5_30")
    spectrumList = loader.fDataSpectrumList
    sbList = spectrumList.FindObject("SideBandAnalysis")

    cname = "SideBandInvMass_QM17"
    canvas = ROOT.TCanvas(cname, cname, 1200, 800)
    globalList.append(canvas)
    canvas.Divide(3, 2)
    PlotSBInvMass(canvas.cd(1), 4, 5, sbList, dptbinList)
    PlotSBInvMass(canvas.cd(2), 6, 7, sbList, dptbinList, True)
    PlotSBInvMass(canvas.cd(3), 10, 12, sbList, dptbinList)
    PlotSBSpectra(canvas.cd(4), 4, 5, sbList)
    PlotSBSpectra(canvas.cd(5), 6, 7, sbList)
    PlotSBSpectra(canvas.cd(6), 10, 12, sbList)

    canvas.cd(3)
    htitle = ROOT.TPaveText(0.22, 0.71, 0.57, 0.86, "NB NDC")
    globalList.append(htitle)
    htitle.SetBorderSize(0)
    htitle.SetFillStyle(0)
    htitle.SetTextFont(43)
    htitle.SetTextSize(20)
    htitle.SetTextAlign(11)
    htitle.AddText("|#eta_{jet}| < 0.5")
    htitle.AddText("5 < #it{p}_{T,ch jet} < 30 GeV/#it{c}")
    htitle.Draw()

    canvas.cd(1)
    paveALICE = ROOT.TPaveText(0.18, 0.72, 0.66, 0.89, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(20)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Preliminary")
    paveALICE.AddText("pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.Draw()

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    SideBandPlot()

    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{0}/{1}.pdf".format(input_path, obj.GetName()))

if __name__ == '__main__':
    main()

    IPython.embed()
