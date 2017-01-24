#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import DataSystematics
import DMesonJetCompare

globalList = []

def PlotUncertainties(results):
    sources = []
    h = results["Uncertainties"]["tot_unc"]
    baseline = h.Clone("{0}_copy".format(h.GetName()))

    print("Source & \\multicolumn{{{0}}}{{c}}{{Uncertainty (\\%)}} \\\\ \\hline".format(baseline.GetNbinsX()))
    print(" & ".join(["\\ptchjet (\\GeVc)"] + ["{0:.0f} - {1:.0f}".format(baseline.GetXaxis().GetBinLowEdge(ibin), baseline.GetXaxis().GetBinUpEdge(ibin)) for ibin in range(1, baseline.GetNbinsX() + 1)]) + "\\\\ \hline")

    h = results["Uncertainties"]["stat_unc"]
    stat_unc = h.Clone("{0}_copy".format(h.GetName()))
    sources.append(stat_unc)

    h = results["Uncertainties"]["tot_rel_syst_unc"]
    tot_rel_syst_unc = h.Clone("{0}_copy".format(h.GetName()))
    sources.append(tot_rel_syst_unc)
    for h in results["Uncertainties"]["PartialSystematicUncertainties"]:
        print(" & ".join([h.GetTitle()] + ["{0:.1f}".format(h.GetBinContent(ibin) * 100) for ibin in range(1, h.GetNbinsX() + 1)]) + "\\\\")
        h_copy = h.Clone("{0}_copy".format(h.GetName()))
        sources.append(h_copy)

    print("\\hline")
    print("Correlated Uncertainty & \\multicolumn{{{0}}}{{c}}{{{1:.1f}}} \\\\".format(baseline.GetNbinsX(), results["Uncertainties"]["correlated_uncertainty"] * 100))
    print("\\hline")
    print(" & ".join([tot_rel_syst_unc.GetTitle()] + ["{0:.1f}".format(tot_rel_syst_unc.GetBinContent(ibin) * 100) for ibin in range(1, tot_rel_syst_unc.GetNbinsX() + 1)]) + "\\\\")
    print("\\hline")
    print(" & ".join([stat_unc.GetTitle()] + ["{0:.1f}".format(stat_unc.GetBinContent(ibin) * 100) for ibin in range(1, stat_unc.GetNbinsX() + 1)]) + "\\\\")
    print("\\hline")
    print(" & ".join([baseline.GetTitle()] + ["{0:.1f}".format(baseline.GetBinContent(ibin) * 100) for ibin in range(1, baseline.GetNbinsX() + 1)]) + "\\\\")

    globalList.extend(sources)
    globalList.append(baseline)
    baseline.Scale(100)
    for s in sources: s.Scale(100)
    comp = DMesonJetCompare.DMesonJetCompare("Uncertainties_QM17")
    comp.fOptSpectrum = "hist"
    comp.fOptSpectrumBaseline = "hist"
    comp.fX1LegSpectrum = 0.15
    comp.fY1LegSpectrum = 0.77
    comp.fLegLineHeight = 0.05
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = False
    comp.fNColsLegSpectrum = 2
    r = comp.CompareSpectra(baseline, sources)
    for obj in r:
        globalList.append(obj)

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

    paveALICE = ROOT.TPaveText(0.14, 0.84, 0.53, 0.95, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(20)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Preliminary, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4, |#eta_{jet}| < 0.5 with D^{0}, #it{p}_{T,D} > 3 GeV/#it{c}")
    paveALICE.Draw()

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
