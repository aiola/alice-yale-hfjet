#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import DMesonJetCompare
import array

globalList = []


def EfficiencyComparison(config):
    fname = "{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"])
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)

    spectrumName = "JetPtDPtSpectrum"
    jetName = "Jet_AKTChargedR040_pt_scheme"
    dmesonName = "D0_D0toKpiCuts"
    prefix = "Prompt_{0}_{1}_{2}".format(dmesonName, jetName, spectrumName)

    DPtBins = [2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 30]

    pt_lim = [(5, 15), (15, 30)]
    histos = []
    for (minJetPt, maxJetPt) in pt_lim:
        recoTruthName = "{0}/{0}_ReconstructedTruth_JetPt_{1}_{2}".format(prefix, minJetPt * 100, maxJetPt * 100)
        truthName = "{0}/{0}_Truth_JetPt_{1}_{2}".format(prefix, minJetPt * 100, maxJetPt * 100)
        h = DMesonJetUtils.GetObject(file, recoTruthName)
        recoTruth = h.Rebin(len(DPtBins) - 1, "{0}_rebin".format(recoTruthName), array.array('d', DPtBins))
        h = DMesonJetUtils.GetObject(file, truthName)
        truth = h.Rebin(len(DPtBins) - 1, "{0}_rebin".format(truthName), array.array('d', DPtBins))
        hist_Ratio = recoTruth.Clone(recoTruthName.replace("RecontructedTruth", "Efficiency"))
        hist_Ratio.Divide(truth)
        hist_Ratio.SetTitle("{} < #it{{p}}_{{T,ch jet}} < {} GeV/#it{{c}}".format(minJetPt, maxJetPt))
        globalList.append(hist_Ratio)
        histos.append(hist_Ratio)

    cname = "EfficiencyVsJetPt_Paper"
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    # comp.fOptSpectrum = "hist"
    # comp.fOptSpectrumBaseline = "hist"
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = None
    comp.fX1LegSpectrum = 0.47
    comp.fX2LegSpectrum = 0.93
    comp.fY1LegSpectrum = 0.38
    comp.fLinUpperSpace = 0.50
    comp.fLegLineHeight = 0.065
    comp.fColors = [ROOT.kOrange + 2, ROOT.kGreen + 2]
    comp.fMarkers = [ROOT.kFullCircle, ROOT.kOpenCircle]
    r = comp.CompareSpectra(histos[0], histos[1:])
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

    h.GetYaxis().SetTitle("D^{0} Efficiency #times Acceptance")
    h.GetXaxis().SetTitle("#it{p}_{T,D} (GeV/#it{c})")
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

    paveALICE = ROOT.TPaveText(0.14, 0.63, 0.53, 0.95, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(21)
    paveALICE.SetTextAlign(13)
    # paveALICE.AddText("ALICE Preliminary")
    paveALICE.AddText("PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Prompt D^{0} #rightarrow K^{-}#pi^{+} and charge conj.")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4")
    paveALICE.AddText("|#eta_{jet}| < 0.5")
    paveALICE.Draw()

    return canvas


def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    f = open("LHC15i2response_Train1399.yaml", 'r')
    config = yaml.load(f)
    f.close()

    canvas = EfficiencyComparison(config)
    canvas.SaveAs("{0}/EfficiencyVsJetPt_Paper.pdf".format(config["input_path"]))
    canvas.SaveAs("{0}/EfficiencyVsJetPt_Paper.C".format(config["input_path"]))


if __name__ == '__main__':

    main()

    IPython.embed()
