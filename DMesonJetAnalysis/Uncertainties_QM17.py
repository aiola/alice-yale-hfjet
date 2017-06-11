#!/usr/local/bin/python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import DataSystematics
import DMesonJetCompare

globalList = []

def PlotUncertainties(results):
    asymm = False  # whether there is any asymmetric uncertainty

    sourcesUp = []
    sourcesLow = []
    colorsUp = []
    colorsLow = []
    h = results["Uncertainties"]["tot_unc_up"]
    tot_unc_up = h.Clone("{0}_copy".format(h.GetName()))
    sourcesUp.append(tot_unc_up)
    colorsUp.append(ROOT.kBlue + 2)

    print("Source & \\multicolumn{{{0}}}{{c}}{{Uncertainty (\\%)}} \\\\ \\hline".format(tot_unc_up.GetNbinsX()))
    print(" & ".join(["\\ptchjet\\ (\\GeVc)"] + ["{0:.0f} - {1:.0f}".format(tot_unc_up.GetXaxis().GetBinLowEdge(ibin), tot_unc_up.GetXaxis().GetBinUpEdge(ibin)) for ibin in range(1, tot_unc_up.GetNbinsX() + 1)]) + "\\\\ \hline")

    h = results["Uncertainties"]["tot_unc_low"]
    tot_unc_low = h.Clone("{0}_copy".format(h.GetName()))
    sourcesLow.append(tot_unc_low)
    colorsLow.append(ROOT.kBlue + 2)

    h = results["Uncertainties"]["stat_unc"]
    stat_unc = h.Clone("{0}_copy".format(h.GetName()))
    sourcesUp.append(stat_unc)
    colorsUp.append(ROOT.kRed + 2)

    h = results["Uncertainties"]["tot_rel_syst_unc_up"]
    tot_rel_syst_unc_up = h.Clone("{0}_copy".format(h.GetName()))
    sourcesUp.append(tot_rel_syst_unc_up)
    colorsUp.append(ROOT.kGreen + 2)

    h = results["Uncertainties"]["tot_rel_syst_unc_low"]
    tot_rel_syst_unc_low = h.Clone("{0}_copy".format(h.GetName()))
    sourcesLow.append(tot_rel_syst_unc_low)
    colorsLow.append(ROOT.kGreen + 2)

    colorsPart = [ROOT.kOrange + 2, ROOT.kAzure + 2, ROOT.kMagenta + 2, ROOT.kYellow + 2, ROOT.kCyan + 2, ROOT.kPink + 1, ROOT.kTeal + 2]
    cont = 0
    for hUp, hLow in zip(results["Uncertainties"]["PartialSystematicUncertaintiesUp"], results["Uncertainties"]["PartialSystematicUncertaintiesLow"]):
        if hLow:
            print(" & ".join(["\multirow{{2}}{{*}}{{{}}}".format(hUp.GetTitle())] + ["+{0:.0f}".format(hUp.GetBinContent(ibin) * 100) for ibin in range(1, h.GetNbinsX() + 1)]) + "\\\\")
            print(" & ".join([" "] + ["-{0:.0f}".format(hLow.GetBinContent(ibin) * 100) for ibin in range(1, hLow.GetNbinsX() + 1)]) + "\\\\")
            asymm = True
        else:
            print(" & ".join([hUp.GetTitle()] + ["{0:.0f}".format(hUp.GetBinContent(ibin) * 100) for ibin in range(1, h.GetNbinsX() + 1)]) + "\\\\")
        h_copy = hUp.Clone("{0}_copy".format(hUp.GetName()))
        sourcesUp.append(h_copy)
        colorsUp.append(colorsPart[cont])

        if hLow:
            h_copy = hLow.Clone("{0}_copy".format(hLow.GetName()))
            sourcesLow.append(h_copy)
            colorsLow.append(colorsPart[cont])

        cont += 1

    print("\\hline")
    print("\\pt-independent Uncertainty & \\multicolumn{{{0}}}{{c}}{{{1:.1f}}} \\\\".format(tot_unc_up.GetNbinsX(), results["Uncertainties"]["fixed_syst_unc"] * 100))
    print("\\hline")
    if asymm:
        print(" & ".join(["\multirow{2}{*}{Total Systematic Uncertainty}"] + ["+{0:.1f}".format(tot_rel_syst_unc_up.GetBinContent(ibin) * 100) for ibin in range(1, tot_rel_syst_unc_up.GetNbinsX() + 1)]) + "\\\\")
        print(" & ".join([" "] + ["-{0:.1f}".format(tot_rel_syst_unc_low.GetBinContent(ibin) * 100) for ibin in range(1, tot_rel_syst_unc_low.GetNbinsX() + 1)]) + "\\\\")
    else:
        print(" & ".join(["Total Systematic Uncertainty"] + ["{0:.1f}".format(tot_rel_syst_unc_up.GetBinContent(ibin) * 100) for ibin in range(1, tot_rel_syst_unc_up.GetNbinsX() + 1)]) + "\\\\")
    print("\\hline")
    print(" & ".join([stat_unc.GetTitle()] + ["{0:.1f}".format(stat_unc.GetBinContent(ibin) * 100) for ibin in range(1, stat_unc.GetNbinsX() + 1)]) + "\\\\")
    print("\\hline")
    if asymm:
        print(" & ".join(["\multirow{2}{*}{Total Uncertainty}"] + ["+{0:.1f}".format(tot_unc_up.GetBinContent(ibin) * 100) for ibin in range(1, tot_unc_up.GetNbinsX() + 1)]) + "\\\\")
        print(" & ".join([" "] + ["-{0:.1f}".format(tot_unc_low.GetBinContent(ibin) * 100) for ibin in range(1, tot_unc_low.GetNbinsX() + 1)]) + "\\\\")
    else:
        print(" & ".join(["Total Uncertainty"] + ["{0:.1f}".format(tot_unc_up.GetBinContent(ibin) * 100) for ibin in range(1, tot_unc_up.GetNbinsX() + 1)]) + "\\\\")
    globalList.extend(sourcesUp)
    globalList.extend(sourcesLow)
    comp = DMesonJetCompare.DMesonJetCompare("CompareUncertainties_{0}".format(config["name"]))
    comp.fOptSpectrum = "hist"
    comp.fOptSpectrumBaseline = "hist"
    comp.fX1LegSpectrum = 0.16
    comp.fX2LegSpectrum = 0.96
    comp.fY1LegSpectrum = 0.79
    comp.fLegLineHeight = 0.05
    comp.fLegTextSize = 19
    comp.fNColsLegSpectrum = 2
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = False
    comp.fLineWidths = [3, 3, 3] + [2] * (len(colorsUp) - 3)
    comp.fColors = colorsUp
    comp.fLinUpperSpace = 1.0
    if asymm:
        comp.fLines = [1] * len(colorsUp)
    else:
        comp.fLines = range(1, len(colorsUp))
    for h in sourcesUp: h.Scale(100)
    r = comp.CompareSpectra(sourcesUp[0], sourcesUp[1:])
    for obj in r:
        globalList.append(obj)
    if asymm:
        comp.fOptSpectrumBaseline = "hist same"
        comp.fColors = colorsLow
        comp.fLines = [2] * len(colorsLow)
        comp.fDoSpectrumLegend = False
        for h in sourcesLow: h.Scale(100)
        r = comp.CompareSpectra(sourcesLow[0], sourcesLow[1:])
        for obj in r:
            globalList.append(obj)

    comp.fLegendSpectra.SetMargin(0.1)
    canvas = comp.fCanvasSpectra
    canvas.SetTicks(1, 1)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    canvas.SetBottomMargin(0.15)
    canvas.cd()

    h = comp.fMainHistogram

    h.GetYaxis().SetTitle("Relative Uncertainty (%)")
    h.GetXaxis().SetTitle("#it{p}_{T,ch jet} (GeV/#it{c})")
    h.GetXaxis().SetTitleFont(43)
    h.GetXaxis().SetTitleSize(26)
    h.GetXaxis().SetLabelFont(43)
    h.GetXaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleFont(43)
    h.GetYaxis().SetTitleSize(26)
    h.GetYaxis().SetLabelFont(43)
    h.GetYaxis().SetLabelSize(22)
    h.GetYaxis().SetTitleOffset(0.9)
    # h.GetYaxis().SetRangeUser(0, 0.59)

    paveALICE = ROOT.TPaveText(0.135, 0.81, 0.53, 0.93, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(20)
    paveALICE.SetTextAlign(12)
    paveALICE.AddText("ALICE Preliminary, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4, |#eta_{jet}| < 0.5 with D^{0}, #it{p}_{T,D} > 3 GeV/#it{c}")
    paveALICE.Draw()

    if asymm:
        legUpLow = ROOT.TLegend(0.55, 0.70, 0.93, 0.60)
        globalList.append(legUpLow)
        legUpLow.SetFillStyle(0)
        legUpLow.SetBorderSize(0)
        legUpLow.SetTextFont(43)
        legUpLow.SetTextSize(20)
        legUpLow.SetMargin(0.1)
        entry = legUpLow.AddEntry(ROOT.nullptr, "Symmetric/Upper Uncertainty", "l")
        entry.SetLineColor(ROOT.kBlack)
        entry.SetLineWidth(2)
        entry.SetLineStyle(1)
        entry = legUpLow.AddEntry(ROOT.nullptr, "Lower Uncertainty", "l")
        entry.SetLineColor(ROOT.kBlack)
        entry.SetLineWidth(2)
        entry.SetLineStyle(2)
        legUpLow.Draw()

    pavePtIndep = ROOT.TPaveText(0.15, 0.49, 0.42, 0.57, "NB NDC")
    globalList.append(pavePtIndep)
    pavePtIndep.SetBorderSize(0)
    pavePtIndep.SetFillStyle(0)
    pavePtIndep.SetTextFont(43)
    pavePtIndep.SetTextSize(18)
    pavePtIndep.SetTextAlign(13)
    pavePtIndep.AddText("Normalization: 3.5% lumi, 1% BR")
    pavePtIndep.Draw()

    return canvas

def main(config):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    histograms = DataSystematics.LoadHistograms(config)
    results = dict()
    results["Variations"] = DataSystematics.CompareVariations(config, histograms)
    results["Uncertainties"] = DataSystematics.GenerateUncertainties(config, histograms)
    canvas = PlotUncertainties(results)
    canvas.SaveAs("{0}/Uncertainties_QM17.pdf".format(config["input_path"]))

if __name__ == '__main__':

    f = open("DataSystematics_LHC10.yaml", 'r')
    config = yaml.load(f)
    f.close()

    main(config)

    IPython.embed()
