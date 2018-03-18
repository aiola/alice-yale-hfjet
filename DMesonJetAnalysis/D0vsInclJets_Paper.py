#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import math

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"


def GetD0JetCrossSection():
    fname = "{0}/JetPtSpectrum_DPt_30_Systematics.root".format(input_path)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "FinalSpectrum/CentralPointsStatisticalUncertainty")
    if not hStat:
        print("Cannot get measured cross section with statistical uncertainty!")
        exit(1)
    hSyst = DMesonJetUtils.GetObject(file, "FinalSpectrum/CentralPointsSystematicUncertainty")
    if not hSyst:
        print("Cannot get measured cross section with systematic uncertainty!")
        exit(1)
    return hStat, hSyst


def GetInclJetCrossSection():
    sigmapp = 73.2  # published x-section in mbarn
    triggerEff = 62.2 / 73.2  # eff of MB OR trigger to be applied to data
    bin0CorrData = 0.91  # ratio event after vertexNcontributors cut to evts after phys selection

    fname = "../obusch/outData_spec_Bayes_combPtH_rebinned.root"
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {}".format(fname))
        exit(1)
    hStat = file.Get("hSpecComb")
    if not hStat:
        print("Cannot get inclusive cross section with statistical uncertainty!")
        exit(1)
    file.Close()

    hStat.Scale(sigmapp * triggerEff)
    hStat.Scale(bin0CorrData)

    fname = "../obusch/outData_spec_Bayes_combPtH.root"
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {}".format(fname))
        exit(1)
    hStat_orig = file.Get("hSpecComb")
    if not hStat_orig:
        print("Cannot get inclusive cross section with statistical uncertainty (not rebinned)!")
        exit(1)
    file.Close()

    hStat_orig.Scale(sigmapp * triggerEff)
    hStat_orig.Scale(bin0CorrData)

    fname = "../obusch/outSys_tot_spec.root"
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {}".format(fname))
        exit(1)
    hSyst_orig = file.Get("grSysErrTotNoUE")
    if not hSyst_orig:
        print("Cannot get inclusive cross section with systematic uncertainty!")
        exit(1)
    file.Close()

    hSyst = ROOT.TGraphErrors(hStat.GetNbinsX())

    for i in range(1, hStat.GetNbinsX() + 1):
        x = hStat.GetXaxis().GetBinCenter(i)
        y = hStat.GetBinContent(i)
        xerr = hStat.GetXaxis().GetBinWidth(i) / 2
        print("hStat Bin {}, x in [{}, {}], y = {}".format(i, x - xerr, x + xerr, y))
        yerr_rel = 0.
        y_check = 0.
        for iorig in range(0, hSyst_orig.GetN()):
            if hSyst_orig.GetX()[iorig] < hStat.GetXaxis().GetBinLowEdge(i): continue
            if hSyst_orig.GetX()[iorig] > hStat.GetXaxis().GetBinUpEdge(i): break
            yerr_rel += hSyst_orig.GetErrorY(iorig) * hStat_orig.GetBinContent(iorig + 1) * hStat_orig.GetXaxis().GetBinWidth(iorig + 1)
            y_check += hStat_orig.GetBinContent(iorig + 1) * hStat_orig.GetXaxis().GetBinWidth(iorig + 1)
            print("hSyst_orig Bin {}, x = {}, err_y_rel = {}".format(iorig, hSyst_orig.GetX()[iorig], hSyst_orig.GetErrorY(iorig)))
            print("hStat_orig Bin {}, x = {}, y = {}, err_y_rel = {}".format(iorig, hStat_orig.GetXaxis().GetBinCenter(iorig + 1), hStat_orig.GetBinContent(iorig + 1), hStat_orig.GetBinError(iorig + 1)))
        yerr_rel /= y_check
        y_check /= hStat.GetXaxis().GetBinWidth(i)
        yerr = yerr_rel * y
        hSyst.SetPoint(i - 1, x, y)
        hSyst.SetPointError(i - 1, xerr, yerr)
        print("Bin {} Tot Error rel {}".format(i, yerr_rel))
        print("Bin {}, tot y {}, check {}\n".format(i, y, y_check))

    return hStat, hSyst


def PlotCrossSections(d0Stat, d0Syst, inclStat, inclSyst):
    cname = "D0VsInclJets_Paper"
    canvas = ROOT.TCanvas(cname, cname, 700, 700)
    globalList.append(canvas)
    canvas.Divide(1, 2)
    padMain = canvas.cd(1)
    padMain.SetPad(0, 0.35, 1, 1)
    padMain.SetBottomMargin(0)
    padMain.SetLeftMargin(0.15)
    padMain.SetTicks(1, 1)
    padRatio = canvas.cd(2)
    padRatio.SetPad(0, 0., 1, 0.35)
    padRatio.SetTopMargin(0)
    padRatio.SetBottomMargin(0.27)
    padRatio.SetLeftMargin(0.15)
    padRatio.SetGridy()
    padRatio.SetTicks(1, 1)

    padMain.cd()
    padMain.SetLogy()
    h = d0Stat.DrawCopy("axis")
    h.GetYaxis().SetRangeUser(5e-6, 4e1)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(26)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleOffset(1.6)
    # h.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]")
    # h.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")

    d0Syst_copy = d0Syst.Clone("{0}_copy".format(d0Syst.GetName()))
    globalList.append(d0Syst_copy)
    d0Syst_copy.Draw("2")
    d0Syst_copy.SetFillStyle(1001)
    d0Syst_copy.SetLineWidth(0)
    d0Syst_copy.SetFillColor(ROOT.kRed - 10)

    d0Stat_copy = d0Stat.DrawCopy("same p e0 x0")
    globalList.append(d0Stat_copy)
    d0Stat_copy.SetLineColor(ROOT.kRed + 2)
    d0Stat_copy.SetMarkerColor(ROOT.kRed + 2)
    d0Stat_copy.SetMarkerStyle(ROOT.kFullCircle)
    d0Stat_copy.SetMarkerSize(1.2)

    inclSyst_copy = inclSyst.Clone("inclSyst_copy")
    globalList.append(inclSyst_copy)
    inclSyst_copy.Draw("2")
    inclSyst_copy.SetFillStyle(1001)
    inclSyst_copy.SetLineWidth(0)
    inclSyst_copy.SetFillColor(ROOT.kCyan - 10)

    inclStat_copy = inclStat.DrawCopy("same p e0 x0")
    globalList.append(inclStat_copy)
    inclStat_copy.SetLineColor(ROOT.kBlue + 2)
    inclStat_copy.SetMarkerColor(ROOT.kBlue + 2)
    inclStat_copy.SetMarkerStyle(ROOT.kFullSquare)
    inclStat_copy.SetMarkerSize(1.2)

    padRatio.cd()

    hRatio = d0Stat_copy.DrawCopy("axis")
    hRatio.GetYaxis().SetTitle("ratio")
    hRatio.GetXaxis().SetTitleFont(43)
    hRatio.GetXaxis().SetTitleSize(26)
    hRatio.GetXaxis().SetLabelFont(43)
    hRatio.GetXaxis().SetLabelSize(22)
    hRatio.GetYaxis().SetTitleFont(43)
    hRatio.GetYaxis().SetTitleSize(26)
    hRatio.GetYaxis().SetLabelFont(43)
    hRatio.GetYaxis().SetLabelSize(22)
    hRatio.GetYaxis().SetTitleOffset(1.4)
    hRatio.GetXaxis().SetTitleOffset(2.9)
    hRatio.GetYaxis().SetRangeUser(0, 0.079)
    hRatio.GetYaxis().SetNdivisions(509)

    ratioSyst = d0Syst_copy.Clone("ratioSyst")
    globalList.append(ratioSyst)
    ratioSyst.Draw("2")

    ratioStat = d0Stat_copy.DrawCopy("same p e0 x0")
    globalList.append(ratioStat)

    xsec_incl_tot = 0.0
    stat_xsec_incl_tot2 = 0.0
    syst_xsec_incl_tot = 0.0

    xsec_d0_tot = 0.0
    stat_xsec_d0_tot2 = 0.0
    syst_xsec_d0_tot = 0.0

    avg_tot = 0.0
    stat_avg_tot = 0.0
    syst_avg_tot = 0.0

    avg_flat = 0.0
    sum_of_w_avg_flat = 0.0
    stat_avg_flat2 = 0.0
    syst_avg_flat2 = 0.0

    for ibin in range(0, ratioSyst.GetN()):
        ratioSyst.SetPoint(ibin, ratioSyst.GetX()[ibin], d0Syst_copy.GetY()[ibin] / inclSyst_copy.GetY()[ibin + 5])
        syst_erry2 = ((d0Syst_copy.GetErrorY(ibin) / d0Syst_copy.GetY()[ibin]) ** 2 + (inclSyst_copy.GetErrorY(ibin + 5) / inclSyst_copy.GetY()[ibin + 5]) ** 2 - 2 * 0.035 ** 2) * ratioSyst.GetY()[ibin] ** 2
        syst_erry = math.sqrt(syst_erry2)
        ratioSyst.SetPointError(ibin, ratioSyst.GetErrorX(ibin), ratioSyst.GetErrorX(ibin), syst_erry, syst_erry)

        xsec_incl_tot += inclStat_copy.GetBinContent(ibin + 6)
        stat_xsec_incl_tot2 += inclStat_copy.GetBinError(ibin + 6) ** 2
        syst_xsec_incl_tot += inclSyst_copy.GetErrorY(ibin + 5)  # take the weighted average of the rel unc

        ratioStat.SetBinContent(ibin + 1, d0Stat_copy.GetBinContent(ibin + 1) / inclStat_copy.GetBinContent(ibin + 6))
        stat_err_y2 = ((d0Stat_copy.GetBinError(ibin + 1) / d0Stat_copy.GetBinContent(ibin + 1)) ** 2 + (inclStat_copy.GetBinError(ibin + 6) / inclStat_copy.GetBinContent(ibin + 6)) ** 2) * ratioStat.GetBinContent(ibin + 1) ** 2
        stat_err_y = math.sqrt(stat_err_y2)
        ratioStat.SetBinError(ibin + 1, stat_err_y)

        tot_err_y2 = stat_err_y2 + syst_erry2
        tot_err_y = math.sqrt(tot_err_y2)

        xsec_d0_tot += d0Stat_copy.GetBinContent(ibin + 1)
        stat_xsec_d0_tot2 += d0Stat_copy.GetBinError(ibin + 1) ** 2
        syst_xsec_d0_tot += d0Syst_copy.GetErrorY(ibin)  # take the weithed average of the rel unc

        if ratioSyst.GetX()[ibin] > 8:
            avg_flat += ratioStat.GetBinContent(ibin + 1) / tot_err_y
            sum_of_w_avg_flat += 1.0 / tot_err_y
            stat_avg_flat2 += stat_err_y2 / tot_err_y2
            syst_avg_flat2 += syst_erry2 / tot_err_y2

        print("Bin {}, x = {}, ratio = {} +/- {} (stat) +/- {} (syst)".format(ibin, ratioSyst.GetX()[ibin], ratioStat.GetBinContent(ibin + 1), stat_err_y, syst_erry))

    stat_xsec_incl_tot = math.sqrt(stat_xsec_incl_tot2)
    stat_xsec_d0_tot = math.sqrt(stat_xsec_d0_tot2)

    avg_tot = xsec_d0_tot / xsec_incl_tot
    stat_avg_tot = math.sqrt(stat_xsec_incl_tot2 / xsec_incl_tot ** 2 + stat_xsec_d0_tot2 / xsec_d0_tot ** 2) * avg_tot

    syst_xsec_incl_tot2 = syst_xsec_incl_tot ** 2
    syst_xsec_d0_tot2 = syst_xsec_d0_tot ** 2

    syst_avg_tot = math.sqrt(syst_xsec_incl_tot2 / xsec_incl_tot ** 2 + syst_xsec_d0_tot2 / xsec_d0_tot ** 2) * avg_tot

    avg_flat /= sum_of_w_avg_flat
    stat_avg_flat = math.sqrt(stat_avg_flat2) / sum_of_w_avg_flat
    syst_avg_flat = math.sqrt(syst_avg_flat2) / sum_of_w_avg_flat

    print("Average ratio in full range: {} +/- {} (stat) +/- {} (syst)".format(avg_tot, stat_avg_tot, syst_avg_tot))
    print("Average ratio in flat region (pt > 8 GeV/c): {} +/- {} (stat) +/- {} (syst)".format(avg_flat, stat_avg_flat, syst_avg_flat))

    padMain.cd()
    leg1 = ROOT.TLegend(0.50, 0.56, 0.80, 0.70, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(23)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    entry = leg1.AddEntry(None, "D^{0} Jets, #it{p}_{T,D} > 3 GeV/#it{c}", "pf")
    entry.SetFillStyle(d0Syst_copy.GetFillStyle())
    entry.SetFillColor(d0Syst_copy.GetFillColor())
    entry.SetLineColor(d0Syst_copy.GetFillColor())
    entry.SetMarkerColor(d0Stat_copy.GetMarkerColor())
    entry.SetMarkerStyle(d0Stat_copy.GetMarkerStyle())
    entry = leg1.AddEntry(None, "Inclusive Jets", "pf")
    entry.SetFillStyle(inclSyst_copy.GetFillStyle())
    entry.SetFillColor(inclSyst_copy.GetFillColor())
    entry.SetLineColor(inclSyst_copy.GetFillColor())
    entry.SetMarkerColor(inclStat_copy.GetMarkerColor())
    entry.SetMarkerStyle(inclStat_copy.GetMarkerStyle())
    leg1.Draw()

    padMain.cd()
    paveALICE = ROOT.TPaveText(0.16, 0.76, 0.55, 0.90, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(22)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4, |#eta_{jet}| < 0.5")
    paveALICE.Draw()

    padRatio.RedrawAxis("g")
    padRatio.RedrawAxis()
    padMain.RedrawAxis()

    return canvas


def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    d0Stat, d0Syst = GetD0JetCrossSection()
    inclStat, inclSyst = GetInclJetCrossSection()
    canvas = PlotCrossSections(d0Stat, d0Syst, inclStat, inclSyst)
    canvas.SaveAs("{}/{}.pdf".format(input_path, canvas.GetName()))
    canvas.SaveAs("{}/{}.C".format(input_path, canvas.GetName()))


if __name__ == '__main__':
    main()

    IPython.embed()
