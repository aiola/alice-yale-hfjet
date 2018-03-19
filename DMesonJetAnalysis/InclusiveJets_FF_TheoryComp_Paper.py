#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import argparse
import math

globalList = []


def GetD0JetCrossSection(input_path, file_name):
    fname = "{}/{}.root".format(input_path, file_name)
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


def GetD0JetTheoryCrossSectionAll(config, axis):
    for t in config["theory"]:
        h = GetD0JetTheoryCrossSection(config["input_path"], t["gen"], t["proc"], t["ts"], config["theory_spectrum"], axis, config["normalize"])
        t["histogram"] = h


def GetD0JetTheoryCrossSection(input_path, gen, proc, ts, spectrum, axis, normalize):
    fname = "{input_path}/FastSim_{gen}_{proc}_{ts}/FastSimAnalysis_ccbar_{gen}_{proc}_{ts}.root".format(input_path=input_path, gen=gen, proc=proc, ts=ts)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    h_orig = DMesonJetUtils.GetObject(file, "D0_MCTruth/Charged_R040/D0_MCTruth_Charged_R040_{spectrum}/D0_MCTruth_Charged_R040_{spectrum}".format(spectrum=spectrum))
    if not h_orig:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)

    h = DMesonJetUtils.Rebin1D(h_orig, axis)
    if normalize:
        h.Scale(1.0 / h.Integral(1, h.GetNbinsX()), "width")
    else:
        h.Scale(0.5, "width")  # particle/antiparticle

    return h


def GetInclusiveJetTheoryCrossSectionAll(config):
    for t in config["theory"]:
        if not t["inclusive"]: continue
        h = GetInclusiveJetTheoryCrossSection(config["input_path"], t["inclusive"]["gen"], t["inclusive"]["proc"], t["inclusive"]["ts"], config["normalize"])
        t["inclusive_histogram"] = h


def GetInclusiveJetTheoryCrossSection(input_path, gen, proc, ts, normalize):
    fname = "{input_path}/FastSim_{gen}_{proc}_{ts}/Jets.root".format(input_path=input_path, gen=gen, proc=proc, ts=ts)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)
    h_orig = DMesonJetUtils.GetObject(file, "JetPt")
    if not h_orig:
        print("Cannot get theory cross section with statistical uncertainty!")
        exit(1)

    h = h_orig.Clone()
    if normalize: h.Scale(1.0 / h.Integral(1, h.GetNbinsX()), "")

    return h


def GetTotalCrossSection(incl_stat_copy, incl_syst_copy):
    for ibin in range(0, ratioSyst.GetN()):
        xsec_incl_tot += incl_stat_copy.GetBinContent(ibin + 6)
        stat_xsec_incl_tot2 += incl_stat_copy.GetBinError(ibin + 6) ** 2
        syst_xsec_incl_tot += incl_syst_copy.GetErrorY(ibin + 5)  # take the weighted average of the rel unc

    stat_xsec_incl_tot = math.sqrt(stat_xsec_incl_tot2)


def PlotCrossSections(d0jet_stat, d0jet_syst, incl_stat, incl_syst, theory, title, logy, miny, maxy, minr, maxr, legx):
    cname = "D0JetVsInclusiveCrossSectionTheoryComp_Paper"
    canvas = ROOT.TCanvas(cname, cname, 700, 900)
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
    if logy: padMain.SetLogy()
    h = d0jet_stat.DrawCopy("axis")
    h.GetYaxis().SetRangeUser(miny, maxy)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(26)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleOffset(1.6)
    # h.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]")
    # h.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")

    d0jet_syst_copy = d0jet_syst.Clone("{0}_copy".format(d0jet_syst.GetName()))
    d0jet_syst_copy.Draw("2")
    d0jet_syst_copy.SetFillStyle(1001)
    d0jet_syst_copy.SetLineWidth(0)
    d0jet_syst_copy.SetFillColor(ROOT.kRed - 10)
    globalList.append(d0jet_syst_copy)

    d0jet_stat_copy = d0jet_stat.DrawCopy("same p e0 x0")
    d0jet_stat_copy.SetLineColor(ROOT.kRed + 2)
    d0jet_stat_copy.SetMarkerColor(ROOT.kRed + 2)
    d0jet_stat_copy.SetMarkerStyle(ROOT.kFullCircle)
    d0jet_stat_copy.SetMarkerSize(1.2)
    globalList.append(d0jet_stat_copy)

    incl_syst_copy = incl_syst.Clone("incl_syst_copy")
    globalList.append(incl_syst_copy)
    incl_syst_copy.Draw("2")
    incl_syst_copy.SetFillStyle(1001)
    incl_syst_copy.SetLineWidth(0)
    incl_syst_copy.SetFillColor(ROOT.kCyan - 10)

    incl_stat_copy = incl_stat.DrawCopy("same p e0 x0")
    globalList.append(incl_stat_copy)
    incl_stat_copy.SetLineColor(ROOT.kBlue + 2)
    incl_stat_copy.SetMarkerColor(ROOT.kBlue + 2)
    incl_stat_copy.SetMarkerStyle(ROOT.kFullSquare)
    incl_stat_copy.SetMarkerSize(1.2)

    padRatio.cd()

    hRatio = d0jet_stat_copy.DrawCopy("axis")
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

    ratioSyst = d0jet_syst_copy.Clone("ratioSyst")
    globalList.append(ratioSyst)
    ratioSyst.Draw("2")

    ratioStat = d0jet_stat_copy.DrawCopy("same p e0 x0")
    globalList.append(ratioStat)

    xsec_incl_tot, stat_xsec_incl_tot2, syst_xsec_incl_tot = GetTotalCrossSection(incl_stat_copy, incl_syst_copy)

    for t in theory:
        if not t["inclusive"]: continue
        padMain.cd()
        d0jet_h = t["histogram"].Clone("_".join([t["gen"], t["proc"]]))
        d0jet_h.Draw("same e0")
        globalList.append(d0jet_h)
        d0jet_h.SetLineColor(t["color"])
        d0jet_h.SetLineStyle(1)
        d0jet_h.SetLineWidth(2)
        d0jet_h.SetMarkerStyle(1)
        d0jet_h.SetMarkerColor(t["color"])
        t["histogram_plot"] = d0jet_h

        incl_h = t["inclusive_histogram"].Clone("_".join([t["inclusive"]["gen"], t["inclusive"]["proc"]]))
        incl_h.Draw("same e0")
        globalList.append(incl_h)
        incl_h.SetLineColor(t["color"])
        incl_h.SetLineStyle(2)
        incl_h.SetLineWidth(2)
        incl_h.SetMarkerStyle(1)
        incl_h.SetMarkerColor(t["color"])
        t["inclusive_histogram_plot"] = incl_h

        padRatio.cd()
        hratio = d0jet_h.Clone("_".join([t["gen"], t["proc"], "over", t["inclusive"]["gen"], t["inclusive"]["proc"]]))
        hratio.Divide(incl_h)
        hratio.Draw("same e0")
        t["ratio_histogram"] = hratio

    padMain.cd()

    y2 = 0.66
    if len(title) > 3: y2 = 0.66 - 0.06 * (len(title) - 3)
    y1 = y2 - 0.036 * len(t) / 2
    if y1 < 0.1: y1 = 0.1
    leg1 = ROOT.TLegend(0.18, y1, 0.85, y2, "", "NB NDC")
    leg1.SetNColumns(2)
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(19)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    for t in theory:
        if "histogram_plot" in t: leg1.AddEntry(t["histogram_plot"], t["title"], "l")
    leg1.Draw()

    leg1 = ROOT.TLegend(legx, y1 - 0.12, legx + 0.30, y1, "", "NB NDC")
    globalList.append(leg1)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(23)
    leg1.SetTextAlign(12)
    leg1.SetMargin(0.2)
    entry = leg1.AddEntry(None, "D^{0} Jets, #it{p}_{T,D} > 3 GeV/#it{c}", "pf")
    entry.SetFillStyle(d0jet_syst_copy.GetFillStyle())
    entry.SetFillColor(d0jet_syst_copy.GetFillColor())
    entry.SetLineColor(d0jet_syst_copy.GetFillColor())
    entry.SetMarkerColor(d0jet_stat_copy.GetMarkerColor())
    entry.SetMarkerStyle(d0jet_stat_copy.GetMarkerStyle())
    entry = leg1.AddEntry(None, "Inclusive Jets", "pf")
    entry.SetFillStyle(incl_syst_copy.GetFillStyle())
    entry.SetFillColor(incl_syst_copy.GetFillColor())
    entry.SetLineColor(incl_syst_copy.GetFillColor())
    entry.SetMarkerColor(incl_stat_copy.GetMarkerColor())
    entry.SetMarkerStyle(incl_stat_copy.GetMarkerStyle())
    leg1.Draw()

    paveALICE = ROOT.TPaveText(0.16, 0.90 - 0.07 * len(title), 0.55, 0.90, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(22)
    paveALICE.SetTextAlign(13)
    for line in title: paveALICE.AddText(line)
    paveALICE.Draw()

    padRatio.RedrawAxis("g")
    padRatio.RedrawAxis()
    padMain.RedrawAxis()

    return canvas


def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    d0jet_stat, d0jet_syst = GetD0JetCrossSection(config["input_path"], config["data"])
    incl_stat, incl_syst = GetInclJetCrossSection()
    GetD0JetTheoryCrossSectionAll(config, d0jet_stat.GetXaxis())
    GetInclusiveJetTheoryCrossSectionAll(config)
    canvas = PlotCrossSections(d0jet_stat, d0jet_syst, incl_stat, incl_syst, config["theory"], config["title_inclusive"], config["logy"],
                               config["miny_inclusive"], config["maxy_inclusive"], config["minr_inclusive"], config["maxr_inclusive"], config["legx_inclusive"])
    canvas.SaveAs("{}/{}.pdf".format(config["input_path"], config["name"]))
    canvas.SaveAs("{}/{}.C".format(config["input_path"], config["name"]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Jet pt spectrum theory comparison.')
    parser.add_argument('yaml', metavar='conf.yaml')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config)

    IPython.embed()
