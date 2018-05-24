#!/usr/bin/env python
# python script to do extract B feed down correction factors

import math
import IPython
import ROOT
import DMesonJetUtils
import LoadInclusiveJetSpectrum

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"

def GetD0JetCrossSectionTheory(kin_cuts):
    ts = 1483386026
    fname = "{input_path}/PromptDJetsPrediction_{ts}.root".format(input_path=input_path, ts=ts)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    hStat = DMesonJetUtils.GetObject(file, "default/JetPtSpectrum{kin_cuts}/GeneratorLevel_JetPtSpectrum".format(kin_cuts=kin_cuts))
    if not hStat:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)
    hSyst = DMesonJetUtils.GetObject(file, "SystematicUncertainty/JetPtSpectrum{kin_cuts}//GeneratorLevel_JetPtSpectrum/GeneratorLevel_JetPtSpectrum_CentralAsymmSyst".format(kin_cuts=kin_cuts))
    if not hSyst:
        print("Cannot get theory cross section lower systematic uncertainty!")
        exit(1)

    # scale for the bin width and the antiparticle factor
    hStat.Scale(0.5, "width")

    # scale for antiparticle factor
    for i in range(0, hSyst.GetN()):
        hSyst.SetPoint(i, hSyst.GetX()[i], hSyst.GetY()[i] / 2)
        hSyst.SetPointError(i, hSyst.GetErrorXlow(i), hSyst.GetErrorXhigh(i), hSyst.GetErrorYlow(i) / 2, hSyst.GetErrorYhigh(i) / 2)

    return hStat, hSyst

def GetInclJetCrossSectionTheory():
    gen = "powheg+pythia6"
    proc = "dijet"
    ts = 1524742846
    file_name = "FastSimAnalysis_inclusive_jets.root"
    fname = "{input_path}/FastSim_{gen}_{proc}_{ts}/{file_name}".format(input_path=input_path, gen=gen, proc=proc, ts=ts, file_name=file_name)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    h_orig = DMesonJetUtils.GetObject(file, "JetPt")
    if not h_orig:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)

    h = h_orig.Clone()
    return h

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
    return LoadInclusiveJetSpectrum.GetCrossSection()


def PlotCrossSections(d0Stat, d0Syst, inclStat, inclSyst, theory_d0Stat, theory_d0Syst, theory_noptcut_d0Stat, theory_inclStat):
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
    d0Syst_copy.SetFillColor(ROOT.kGray)

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
    inclSyst_copy.SetFillColor(ROOT.kGray)

    inclStat_copy = inclStat.DrawCopy("same p e0 x0")
    globalList.append(inclStat_copy)
    inclStat_copy.SetLineColor(ROOT.kMagenta + 2)
    inclStat_copy.SetMarkerColor(ROOT.kMagenta + 2)
    inclStat_copy.SetMarkerStyle(ROOT.kFullSquare)
    inclStat_copy.SetMarkerSize(1.2)

    theory_noptcut_d0Stat_copy = theory_noptcut_d0Stat.DrawCopy("same p e0 x0")
    globalList.append(theory_noptcut_d0Stat_copy)
    theory_noptcut_d0Stat_copy.SetLineColor(ROOT.kGray + 3)
    theory_noptcut_d0Stat_copy.SetMarkerColor(ROOT.kGray + 3)
    theory_noptcut_d0Stat_copy.SetMarkerStyle(ROOT.kOpenCircle)
    theory_noptcut_d0Stat_copy.SetMarkerSize(1.3)

    theory_d0Stat_copy = theory_d0Stat.DrawCopy("same p e0 x0")
    globalList.append(theory_d0Stat_copy)
    theory_d0Stat_copy.SetLineColor(ROOT.kBlue + 2)
    theory_d0Stat_copy.SetMarkerColor(ROOT.kBlue + 2)
    theory_d0Stat_copy.SetMarkerStyle(ROOT.kOpenCircle)
    theory_d0Stat_copy.SetMarkerSize(1.2)

    theory_d0Syst_copy = theory_d0Syst.Clone("theory_d0_Syst_copy")
    globalList.append(theory_d0Syst_copy)
    theory_d0Syst_copy.Draw("2")
    theory_d0Syst_copy.SetFillStyle(0)
    theory_d0Syst_copy.SetLineWidth(2)
    theory_d0Syst_copy.SetLineColor(ROOT.kBlue + 2)

    theory_inclStat_copy = theory_inclStat.DrawCopy("same p e0 x0")
    globalList.append(theory_inclStat_copy)
    theory_inclStat_copy.SetLineColor(ROOT.kGreen + 2)
    theory_inclStat_copy.SetMarkerColor(ROOT.kGreen + 2)
    theory_inclStat_copy.SetMarkerStyle(ROOT.kOpenSquare)
    theory_inclStat_copy.SetMarkerSize(1.2)

    padRatio.cd()

    hRatio = d0Stat_copy.DrawCopy("axis")
    hRatio.GetYaxis().SetTitle("D^{0} jets / inclusive")
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
        ratioSyst.SetPoint(ibin, ratioSyst.GetX()[ibin], d0Syst_copy.GetY()[ibin] / inclSyst_copy.GetY()[ibin])
        syst_erry2 = ((d0Syst_copy.GetErrorY(ibin) / d0Syst_copy.GetY()[ibin]) ** 2 + (inclSyst_copy.GetErrorY(ibin) / inclSyst_copy.GetY()[ibin]) ** 2 - 2 * 0.035 ** 2) * ratioSyst.GetY()[ibin] ** 2
        syst_erry = math.sqrt(syst_erry2)
        ratioSyst.SetPointError(ibin, ratioSyst.GetErrorX(ibin), ratioSyst.GetErrorX(ibin), syst_erry, syst_erry)

        xsec_incl_tot += inclStat_copy.GetBinContent(ibin + 1)
        stat_xsec_incl_tot2 += inclStat_copy.GetBinError(ibin + 1) ** 2
        syst_xsec_incl_tot += inclSyst_copy.GetErrorY(ibin)  # take the weighted average of the rel unc

        ratioStat.SetBinContent(ibin + 1, d0Stat_copy.GetBinContent(ibin + 1) / inclStat_copy.GetBinContent(ibin + 1))
        stat_err_y2 = ((d0Stat_copy.GetBinError(ibin + 1) / d0Stat_copy.GetBinContent(ibin + 1)) ** 2 + (inclStat_copy.GetBinError(ibin + 1) / inclStat_copy.GetBinContent(ibin + 1)) ** 2) * ratioStat.GetBinContent(ibin + 1) ** 2
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
    theory_d0Stat, theory_d0Syst = GetD0JetCrossSectionTheory("_DPt_30")
    theory_noptcut_d0Stat, theory_noptcut_d0Syst = GetD0JetCrossSectionTheory("")
    theory_inclStat = GetInclJetCrossSectionTheory()
    canvas = PlotCrossSections(d0Stat, d0Syst,
                               inclStat, inclSyst,
                               theory_d0Stat, theory_d0Syst,
                               theory_noptcut_d0Stat,
                               theory_inclStat)
    canvas.SaveAs("{}/{}.pdf".format(input_path, canvas.GetName()))
    canvas.SaveAs("{}/{}.C".format(input_path, canvas.GetName()))


if __name__ == '__main__':
    main()

    IPython.embed()
