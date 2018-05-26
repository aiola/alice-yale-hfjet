#!/usr/bin/env python
# python script to compare measured  D0 and inclusive jet cross sections with theory
# us it with
# JetPtSpectrum_FullChargeComp_Paper.yaml (compare full/charged)
# JetPtSpectrum_TheoryComp_Paper.yaml

import argparse
import math
import yaml
import IPython
import ROOT
import DMesonJetUtils
import LoadInclusiveJetSpectrum
import LoadTheoryCrossSections

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
    return LoadInclusiveJetSpectrum.GetCrossSection()

def GetD0JetTheoryCrossSectionAll(config, axis):
    return LoadTheoryCrossSections.GetD0JetTheoryCrossSectionAll(config, axis, True)

def GetInclusiveJetTheoryCrossSectionAll(config):
    return LoadTheoryCrossSections.GetInclusiveJetTheoryCrossSectionAll(config)

def PlotCrossSections(d0jet_stat, d0jet_syst, incl_stat, incl_syst, theory, title, logy, miny, maxy, minr, maxr, legx, no_pt_cut):
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
    h.GetYaxis().SetTitleSize(23)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleOffset(1.6)
    # h.GetYaxis().SetTitle("#frac{d^{2}#sigma}{d#it{p}_{T}d#eta} [mb (GeV/#it{c})^{-1}]")
    # h.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")

    d0jet_syst_copy = d0jet_syst.Clone("{0}_copy".format(d0jet_syst.GetName()))
    d0jet_syst_copy.Draw("2")
    d0jet_syst_copy.SetFillStyle(1001)
    d0jet_syst_copy.SetLineWidth(0)
    d0jet_syst_copy.SetFillColor(ROOT.kGray)
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
    incl_syst_copy.SetFillColor(ROOT.kGray)

    incl_stat_copy = incl_stat.DrawCopy("same p e0 x0")
    globalList.append(incl_stat_copy)
    incl_stat_copy.SetLineColor(ROOT.kMagenta + 2)
    incl_stat_copy.SetMarkerColor(ROOT.kMagenta + 2)
    incl_stat_copy.SetMarkerStyle(ROOT.kFullSquare)
    incl_stat_copy.SetMarkerSize(1.2)

    padRatio.cd()

    hRatio = d0jet_stat_copy.DrawCopy("axis")
    hRatio.GetYaxis().SetTitle("R(#it{p}_{T}) / #Delta#it{p}_{T} (GeV/#it{c})^{-1}")
    hRatio.GetXaxis().SetTitleFont(43)
    hRatio.GetXaxis().SetTitleSize(23)
    hRatio.GetXaxis().SetLabelFont(43)
    hRatio.GetXaxis().SetLabelSize(22)
    hRatio.GetYaxis().SetTitleFont(43)
    hRatio.GetYaxis().SetTitleSize(23)
    hRatio.GetYaxis().SetLabelFont(43)
    hRatio.GetYaxis().SetLabelSize(22)
    hRatio.GetYaxis().SetTitleOffset(1.9)
    hRatio.GetXaxis().SetTitleOffset(2.9)
    hRatio.GetYaxis().SetRangeUser(minr, maxr)
    hRatio.GetYaxis().SetNdivisions(509)

    ratioSyst = d0jet_syst_copy.Clone("ratioSyst")
    globalList.append(ratioSyst)
    ratioSyst.Draw("2")

    ratioStat = d0jet_stat_copy.DrawCopy("same p e0 x0")
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
        ratioSyst.SetPoint(ibin, ratioSyst.GetX()[ibin], d0jet_syst_copy.GetY()[ibin] / incl_syst_copy.GetY()[ibin])
        syst_erry2 = ((d0jet_syst_copy.GetErrorY(ibin) / d0jet_syst_copy.GetY()[ibin]) ** 2 + (incl_syst_copy.GetErrorY(ibin) / incl_syst_copy.GetY()[ibin]) ** 2 - 2 * 0.035 ** 2) * ratioSyst.GetY()[ibin] ** 2
        syst_erry = math.sqrt(syst_erry2)
        ratioSyst.SetPointError(ibin, ratioSyst.GetErrorX(ibin), ratioSyst.GetErrorX(ibin), syst_erry, syst_erry)

        xsec_incl_tot += incl_stat_copy.GetBinContent(ibin + 1)
        stat_xsec_incl_tot2 += incl_stat_copy.GetBinError(ibin + 1) ** 2
        syst_xsec_incl_tot += incl_syst_copy.GetErrorY(ibin)  # take the weighted average of the rel unc

        ratioStat.SetBinContent(ibin + 1, d0jet_stat_copy.GetBinContent(ibin + 1) / incl_stat_copy.GetBinContent(ibin + 1))
        stat_err_y2 = ((d0jet_stat_copy.GetBinError(ibin + 1) / d0jet_stat_copy.GetBinContent(ibin + 1)) ** 2 + (incl_stat_copy.GetBinError(ibin + 1) / incl_stat_copy.GetBinContent(ibin + 1)) ** 2) * ratioStat.GetBinContent(ibin + 1) ** 2
        stat_err_y = math.sqrt(stat_err_y2)
        ratioStat.SetBinError(ibin + 1, stat_err_y)

        tot_err_y2 = stat_err_y2 + syst_erry2
        tot_err_y = math.sqrt(tot_err_y2)

        xsec_d0_tot += d0jet_stat_copy.GetBinContent(ibin + 1)
        stat_xsec_d0_tot2 += d0jet_stat_copy.GetBinError(ibin + 1) ** 2
        syst_xsec_d0_tot += d0jet_syst_copy.GetErrorY(ibin)  # take the weithed average of the rel unc

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

    for t in theory:
        if not t["active"]: continue
        if not t["inclusive"]: continue
        padMain.cd()
        d0jet_h = t["histogram"].Clone("_".join([t["gen"], t["proc"]]))
        d0jet_h.Draw("same e0 p x0")
        globalList.append(d0jet_h)
        d0jet_h.SetLineColor(t["color"])
        d0jet_h.SetLineStyle(1)
        d0jet_h.SetLineWidth(1)
        d0jet_h.SetMarkerStyle(ROOT.kOpenCircle)
        d0jet_h.SetMarkerColor(t["color"])
        t["histogram_plot"] = d0jet_h

        incl_h = t["inclusive_histogram"].Clone("_".join([t["inclusive"]["gen"], t["inclusive"]["proc"]]))
        incl_h.Draw("same e0 p x0")
        globalList.append(incl_h)
        incl_h.SetLineColor(t["color"])
        incl_h.SetLineStyle(1)
        incl_h.SetLineWidth(1)
        incl_h.SetMarkerStyle(ROOT.kOpenSquare)
        incl_h.SetMarkerColor(t["color"])
        t["inclusive_histogram_plot"] = incl_h

        padRatio.cd()
        hratio = d0jet_h.Clone("_".join([t["gen"], t["proc"], "over", t["inclusive"]["gen"], t["inclusive"]["proc"]]))
        hratio.Divide(incl_h)
        hratio.Draw("same e0")
        t["ratio_histogram"] = hratio


        if no_pt_cut:
            padMain.cd()
            d0jet_no_pt_cut_h = t["histogram_no_pt_cut"].Clone("_".join([t["gen"], t["proc"]]))
            d0jet_no_pt_cut_h.Draw("same e0 p x0")
            globalList.append(d0jet_no_pt_cut_h)
            d0jet_no_pt_cut_h.SetLineColor(t["color"])
            d0jet_no_pt_cut_h.SetLineStyle(1)
            d0jet_no_pt_cut_h.SetLineWidth(1)
            d0jet_no_pt_cut_h.SetMarkerStyle(ROOT.kOpenTriangleUp)
            d0jet_no_pt_cut_h.SetMarkerColor(t["color"])
            t["histogram_no_pt_cut_plot"] = d0jet_no_pt_cut_h

            padRatio.cd()
            hratio = d0jet_no_pt_cut_h.Clone("_".join([t["gen"], t["proc"], "over", t["inclusive"]["gen"], t["inclusive"]["proc"], "no_pt_cut"]))
            hratio.Divide(incl_h)
            hratio.Draw("same e0")
            t["ratio_no_pt_cut_histogram"] = hratio

    padMain.cd()

    y2 = 0.90 - 0.07 * len(title)
    paveALICE = ROOT.TPaveText(0.16, y2, 0.55, 0.90, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(22)
    paveALICE.SetTextAlign(13)
    for line in title: paveALICE.AddText(line)
    paveALICE.Draw()

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
        if "inclusive_histogram_plot" in t: leg1.AddEntry(t["inclusive_histogram_plot"], "Inclusive, {}".format(t["title"]), "p")
        if "histogram_plot" in t: leg1.AddEntry(t["histogram_plot"], t["title"], "p")
        if "histogram_no_pt_cut_plot" in t: leg1.AddEntry(t["histogram_no_pt_cut_plot"], "#it{{p}}_{{T,D}} > 0, {}".format(t["title"]), "p")
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

    padRatio.RedrawAxis("g")
    padRatio.RedrawAxis()
    padMain.RedrawAxis()

    return canvas


def main(config, no_pt_cut):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    d0jet_stat, d0jet_syst = GetD0JetCrossSection(config["input_path"], config["data"])
    incl_stat, incl_syst = GetInclJetCrossSection()
    GetD0JetTheoryCrossSectionAll(config, d0jet_stat.GetXaxis())
    GetInclusiveJetTheoryCrossSectionAll(config)
    canvas = PlotCrossSections(d0jet_stat, d0jet_syst, incl_stat, incl_syst, config["theory"], config["title_inclusive"], config["logy"],
                               config["miny_inclusive"], config["maxy_inclusive"], config["minr_inclusive"], config["maxr_inclusive"], config["legx_inclusive"],
                               no_pt_cut)
    canvas.SaveAs("{}/{}.pdf".format(config["input_path"], config["name_inclusive"]))
    canvas.SaveAs("{}/{}.C".format(config["input_path"], config["name_inclusive"]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Jet pt spectrum theory comparison.')
    parser.add_argument('yaml', metavar='conf.yaml')
    parser.add_argument("--no-pt-cut", action='store_const',
                        default=False, const=True,
                        help='Include theory without D pt cut.')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.no_pt_cut)

    IPython.embed()
