#!/usr/bin/env python

import subprocess
import IPython
import ROOT
import RawYieldSpectrumLoader

globalList = []

input_path = "../../workdir"

jet_def = "Charged_R040_pt_scheme"

def PlotSBInvMass(pad, dmeson, ptmin, ptmax, jetptmin, jetptmax, sbList, dptbinList, refl, plotleg1=False, plotleg2=False):
    pad.SetTicks(1, 1)
    pad.SetLeftMargin(0.22)
    pad.SetRightMargin(0.02)
    pad.SetTopMargin(0.03)
    pad.SetBottomMargin(0.14)
    objname = "InvMassSBWindow_AnyINT_{}_DPt_{:.0f}_{:.0f}".format(dmeson, ptmin * 100, ptmax * 100)
    sbHist = sbList.FindObject(objname)
    if not sbHist:
        print("Could not find object '{}' in list '{}'".format(objname, sbList.GetName()))
        sbList.ls()
        exit(1)
    h = sbHist.DrawCopy("axis")
    h.GetXaxis().SetRangeUser(1.709, 2.07)

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

    objname = "InvMassSigWindow_AnyINT_{}_DPt_{:.0f}_{:.0f}".format(dmeson, ptmin * 100, ptmax * 100)
    sigHist = sbList.FindObject(objname)
    if not sigHist:
        print("Could not find object '{}' in list '{}'".format(objname, sbList.GetName()))
        sbList.ls()
        exit(1)
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

    objname = "InvMass_AnyINT_{}_DPt_{:.0f}_{:.0f}".format(dmeson, ptmin * 100, ptmax * 100)
    invMassHist = dptbinList.FindObject(objname)
    if not invMassHist:
        print("Could not find object '{}' in list '{}'".format(objname, sbList.GetName()))
        sbList.ls()
        exit(1)
    invMassHist_copy = invMassHist.DrawCopy("p0 x0 same")
    globalList.append(invMassHist_copy)
    invMassHist_copy.SetLineColor(ROOT.kBlue + 3)
    invMassHist_copy.SetMarkerColor(ROOT.kBlue + 3)
    invMassHist_copy.SetMarkerStyle(ROOT.kFullCircle)
    invMassHist_copy.SetMarkerSize(1.0)

    objname = "InvMass_AnyINT_{}_DPt_{:.0f}_{:.0f}_fitter{}".format(dmeson, ptmin * 100, ptmax * 100, refl)
    fitter = dptbinList.FindObject(objname)
    if not fitter:
        print("Could not find object '{}' in list '{}'".format(objname, sbList.GetName()))
        sbList.ls()
        exit(1)
    globalList.append(fitter)
    fitter.Draw("same")

    diff = invMassHist_copy.GetMaximum() - invMassHist_copy.GetMinimum()
    miny = invMassHist_copy.GetMinimum() - 0.2 * diff
    if miny < 0: miny = 0
    maxy = invMassHist_copy.GetMaximum() + 1.1 * diff
    h.SetMaximum(maxy)
    h.SetMinimum(miny)

    h.GetYaxis().SetTitle("Counts / (6 MeV/#it{c}^{2})")
    h.GetXaxis().SetTitle("#it{M}(K#pi) (GeV/#it{c}^{2})")
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleOffset(1.8)
    h.GetXaxis().SetTitleSize(23)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelOffset(0.009)
    h.GetXaxis().SetLabelSize(18)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleOffset(3.4)
    h.GetYaxis().SetTitleSize(23)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelOffset(0.009)
    h.GetYaxis().SetLabelSize(23)

    binTitle = "{0:.0f} < #it{{p}}_{{T,D}} < {1:.0f} GeV/#it{{c}}".format(ptmin, ptmax)
    binTitle2 = "{0:.0f} < #it{{p}}_{{T,jet}}^{{ch}} < {1:.0f} GeV/#it{{c}}".format(jetptmin, jetptmax)
    htitle = ROOT.TPaveText(0.22, 0.795, 0.98, 0.945, "NB NDC")
    htitle.SetBorderSize(0)
    htitle.SetFillStyle(0)
    htitle.SetTextFont(43)
    htitle.SetTextSize(19)
    htitle.SetTextAlign(22)
    htitle.AddText(binTitle)
    htitle.AddText(binTitle2)
    htitle.Draw()
    globalList.append(htitle)

    fitInfo = ROOT.TPaveText(0.25, 0.705, 0.58, 0.795, "NB NDC")
    fitInfo.SetBorderSize(0)
    fitInfo.SetFillStyle(0)
    fitInfo.SetTextFont(43)
    fitInfo.SetTextSize(18)
    fitInfo.SetTextAlign(12)
    fitInfo.AddText(fitter.GetSignalMeanString().Data())
    fitInfo.AddText(fitter.GetSignalWidthString().Data())
    fitInfo.Draw()
    globalList.append(fitInfo)

    if plotleg1:
        leg = ROOT.TLegend(0.25, 0.48, 0.58, 0.67, "", "NB NDC")
        globalList.append(leg)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(16)
        leg.SetTextAlign(12)
        leg.AddEntry(invMassHist_copy, "D^{0}-Jet Candidates", "p")
        leg.AddEntry(fitter.GetFitFunction(), "Fit Sig+Bkg+Refl", "l")
        leg.AddEntry(fitter.GetBkgFunction(), "Fit Bkg+Refl", "l")
        leg.AddEntry(fitter.GetBkgFunctionNoRefl(), "Fit Bkg-only", "l")
        leg.Draw()

    if plotleg2:
        leg = ROOT.TLegend(0.25, 0.55, 0.58, 0.67, "", "NB NDC")
        globalList.append(leg)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(18)
        leg.SetTextAlign(12)
        leg.AddEntry(sigHist_copy, "Peak Region", "f")
        leg.AddEntry(sbHist_copy_l, "Side Bands", "f")
        leg.Draw()
    
    pad.RedrawAxis()

    return [invMassHist_copy]

def PlotSBSpectra(pad, dmeson, kincuts, ptmin, ptmax, sbList, refl, plotleg=False):
    pad.SetTicks(1, 1)
    # pad.SetLogy()
    pad.SetLeftMargin(0.22)
    pad.SetRightMargin(0.02)
    pad.SetTopMargin(0.01)
    pad.SetBottomMargin(0.16)

    objname = "{}_{}_JetZSpectrum_{}_SideBand{}_SideBandWindow_DPt_{:.0f}_{:.0f}".format(dmeson, jet_def, kincuts, refl, ptmin * 100, ptmax * 100)
    print(objname)
    sbHist = sbList.FindObject(objname)
    h = sbHist.DrawCopy("axis")
    h.GetXaxis().SetTitle("#it{z}_{||}^{ch}")
    globalList.append(h)

    sbHist_copy = sbHist.DrawCopy("p0 same")
    globalList.append(sbHist_copy)
    sbHist_copy.SetMarkerColor(ROOT.kGreen + 2)
    sbHist_copy.SetLineColor(ROOT.kGreen + 2)
    sbHist_copy.SetMarkerStyle(ROOT.kOpenCircle)
    sbHist_copy.SetMarkerSize(1.3)

    sigHist = sbList.FindObject("{}_{}_JetZSpectrum_{}_SideBand{}_SignalWindow_DPt_{:.0f}_{:.0f}".format(dmeson, jet_def, kincuts, refl, ptmin * 100, ptmax * 100))

    sigHist_copy = sigHist.DrawCopy("p0 same")
    globalList.append(sigHist_copy)
    sigHist_copy.SetMarkerColor(ROOT.kRed + 2)
    sigHist_copy.SetLineColor(ROOT.kRed + 2)
    sigHist_copy.SetMarkerStyle(ROOT.kOpenSquare)
    sigHist_copy.SetMarkerSize(1.3)

    subHist_copy = sigHist.DrawCopy("p0 same")
    subHist_copy.Add(sbHist_copy, -1)
    subHist_copy.SetName(subHist_copy.GetName().replace("SignalWindow","Signal"))
    globalList.append(subHist_copy)
    subHist_copy.SetMarkerColor(ROOT.kBlue + 2)
    subHist_copy.SetLineColor(ROOT.kBlue + 2)
    subHist_copy.SetMarkerStyle(ROOT.kOpenDiamond)
    subHist_copy.SetMarkerSize(1.7)

    h.GetYaxis().SetTitle("Raw Yield")
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleOffset(1.9)
    h.GetXaxis().SetTitleSize(23)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelOffset(0.009)
    h.GetXaxis().SetLabelSize(18)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleOffset(3.8)
    h.GetYaxis().SetTitleSize(23)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelOffset(0.009)
    h.GetYaxis().SetLabelSize(23)

    # miny = min([DMesonJetUtils.FindMinimum(sbHist_copy), DMesonJetUtils.FindMinimum(sigHist_copy), DMesonJetUtils.FindMinimum(subHist_copy)])
    # maxy = max([DMesonJetUtils.FindMaximum(sbHist_copy), DMesonJetUtils.FindMaximum(sigHist_copy), DMesonJetUtils.FindMaximum(subHist_copy)])
    # miny /= 2
    # maxy *= 3
    diff = sigHist_copy.GetMaximum() - sigHist_copy.GetMinimum()
    miny = sigHist_copy.GetMinimum() - 0.1 * diff
    maxy = sigHist_copy.GetMaximum() + 0.5 * diff
    h.SetMaximum(maxy)
    h.SetMinimum(miny)

    if plotleg:
        leg = ROOT.TLegend(0.475, 0.77, 0.855, 0.95, "", "NB NDC")
        globalList.append(leg)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(19)
        leg.SetTextAlign(13)
        leg.AddEntry(sigHist_copy, "Peak Region (PR)", "p")
        leg.AddEntry(sbHist_copy, "Side Bands (SB)", "p")
        leg.AddEntry(subHist_copy, "PR - SB", "p")
        leg.Draw()
    
    return [sbHist_copy, sigHist_copy, subHist_copy]

def SideBandPlot():
    loader = RawYieldSpectrumLoader.RawYieldSpectrumLoader(input_path, "Jets_EMC_pp_1116_1117_1118_1119", "LHC10_Train1116")
    loader.fDMeson = "D0_D0toKpiCuts"
    loader.fJetType = "Charged"
    loader.fJetRadius = "R040"
    loader.fJetRecoScheme = "pt_scheme"
    loader.fVariableName = "JetZ"
    loader.fRawYieldMethod = "SideBand"
    loader.fUseReflections = True

    loader.fKinematicCuts = "DPt_20_JetPt_5_15"
    loader.LoadDataListFromDMesonJetAnalysis()
    dptbinList_low = loader.fDataJetList.FindObject("{}_{}_DPtBins_JetPt_5_15".format(loader.fDMeson, jet_def))
    spectrumList_low = loader.fDataSpectrumList
    sbList_low = spectrumList_low.FindObject("SideBandAnalysis")

    loader.fKinematicCuts = "DPt_60_JetPt_15_30"
    loader.LoadDataListFromDMesonJetAnalysis()
    dptbinList_high = loader.fDataJetList.FindObject("{}_{}_DPtBins_JetPt_15_30".format(loader.fDMeson, jet_def))
    spectrumList_high = loader.fDataSpectrumList
    sbList_high = spectrumList_high.FindObject("SideBandAnalysis")

    cname = "SideBandInvMass_Z_Paper"
    canvas = ROOT.TCanvas(cname, cname, 1200, 800)
    globalList.append(canvas)
    canvas.Divide(3, 2)
    PlotSBInvMass(canvas.cd(1), loader.fDMeson, 4, 5, 5, 15, sbList_low, dptbinList_low, "_DoubleGaus", True, False)
    PlotSBInvMass(canvas.cd(2), loader.fDMeson, 6, 7, 5, 15, sbList_low, dptbinList_low, "_DoubleGaus", False, True)
    PlotSBInvMass(canvas.cd(3), loader.fDMeson, 6, 12, 15, 30, sbList_high, dptbinList_high, "_DoubleGaus", False, False)
    PlotSBSpectra(canvas.cd(4), loader.fDMeson, "DPt_20_JetPt_5_15", 4, 5, sbList_low, "_DoubleGaus")
    PlotSBSpectra(canvas.cd(5), loader.fDMeson, "DPt_20_JetPt_5_15", 6, 7, sbList_low, "_DoubleGaus")
    PlotSBSpectra(canvas.cd(6), loader.fDMeson, "DPt_60_JetPt_15_30", 6, 12, sbList_high, "_DoubleGaus", True)

    canvas.cd(5)
    htitle = ROOT.TPaveText(0.22, 0.78, 0.92, 0.92, "NB NDC")
    globalList.append(htitle)
    htitle.SetBorderSize(0)
    htitle.SetFillStyle(0)
    htitle.SetTextFont(43)
    htitle.SetTextSize(17)
    htitle.SetTextAlign(11)
    htitle.AddText("Charged Jets Anti-#it{k}_{T}, #it{R} = 0.4, |#it{#eta}_{jet}| < 0.5")
    htitle.AddText("D^{0} #rightarrow K^{#font[122]{-}}#pi^{+} and charge conj.")
    htitle.Draw()

    canvas.cd(4)
    paveALICE = ROOT.TPaveText(0.32, 0.75, 0.66, 0.92, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(20)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("#it{This Thesis}")
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
            obj.SaveAs("{0}/{1}.C".format(input_path, obj.GetName()))


if __name__ == '__main__':
    main()

    IPython.embed()
