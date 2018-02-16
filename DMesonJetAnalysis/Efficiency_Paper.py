#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import DMesonJetCompare
import array

globalList = []


def EfficiencyComparison(config_c, config_b):
    fname = "{0}/{1}/{2}.root".format(config_c["input_path"], config_c["train"], config_c["name"])
    file_c = ROOT.TFile(fname)
    if not file_c or file_c.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)

    fname = "{0}/{1}/{2}.root".format(config_b["input_path"], config_b["train"], config_b["name"])
    file_b = ROOT.TFile(fname)
    if not file_b or file_b.IsZombie():
        print("Could not open file {0}".format(fname))
        exit(1)

    spectrumName = "JetPtDPtSpectrum"
    jetName = "Jet_AKTChargedR040_pt_scheme"
    dmesonName = "D0_D0toKpiCuts"
    prefix_c = "Prompt_{0}_{1}_{2}".format(dmesonName, jetName, spectrumName)
    prefix_b = "NonPrompt_{0}_{1}_{2}".format(dmesonName, jetName, spectrumName)

    DPtBins = [2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 30]

    minJetPt = 5
    maxJetPt = 30
    recoTruthName_c = "{0}/{0}_ReconstructedTruth_JetPt_{1}_{2}".format(prefix_c, minJetPt * 100, maxJetPt * 100)
    truthName_c = "{0}/{0}_Truth_JetPt_{1}_{2}".format(prefix_c, minJetPt * 100, maxJetPt * 100)
    recoTruthName_b = "{0}/{0}_ReconstructedTruth_JetPt_{1}_{2}".format(prefix_b, minJetPt * 100, maxJetPt * 100)
    truthName_b = "{0}/{0}_Truth_JetPt_{1}_{2}".format(prefix_b, minJetPt * 100, maxJetPt * 100)
    h = DMesonJetUtils.GetObject(file_c, recoTruthName_c)
    cRecoTruth = h.Rebin(len(DPtBins) - 1, "{0}_c_rebin".format(recoTruthName_c), array.array('d', DPtBins))
    h = DMesonJetUtils.GetObject(file_c, truthName_c)
    cTruth = h.Rebin(len(DPtBins) - 1, "{0}_c_rebin".format(truthName_c), array.array('d', DPtBins))
    h = DMesonJetUtils.GetObject(file_b, recoTruthName_b)
    bRecoTruth = h.Rebin(len(DPtBins) - 1, "{0}_b_rebin".format(recoTruthName_b), array.array('d', DPtBins))
    h = DMesonJetUtils.GetObject(file_b, truthName_b)
    bTruth = h.Rebin(len(DPtBins) - 1, "{0}_b_rebin".format(truthName_b), array.array('d', DPtBins))
    c_hist_Ratio = cRecoTruth.Clone(recoTruthName_c.replace("RecontructedTruth", "Efficiency"))
    c_hist_Ratio.Divide(cTruth)
    c_hist_Ratio.SetTitle("D^{0}, Prompt")
    b_hist_Ratio = bRecoTruth.Clone(recoTruthName_b.replace("RecontructedTruth", "Efficiency"))
    b_hist_Ratio.Divide(bTruth)
    b_hist_Ratio.SetTitle("D^{0}, Non-Prompt")
    globalList.append(c_hist_Ratio)
    globalList.append(b_hist_Ratio)

    cname = "Efficiency_Paper"
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    # comp.fOptSpectrum = "hist"
    # comp.fOptSpectrumBaseline = "hist"
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = None
    comp.fX1LegSpectrum = 0.57
    comp.fX2LegSpectrum = 0.93
    comp.fY1LegSpectrum = 0.38
    comp.fLinUpperSpace = 0.50
    comp.fLegLineHeight = 0.05
    comp.fColors = [ROOT.kRed + 2, ROOT.kBlue + 2]
    comp.fMarkers = [ROOT.kFullCircle, ROOT.kFullSquare]
    r = comp.CompareSpectra(c_hist_Ratio, [b_hist_Ratio])
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

    paveALICE = ROOT.TPaveText(0.14, 0.60, 0.53, 0.95, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(21)
    paveALICE.SetTextAlign(13)
    # paveALICE.AddText("ALICE Preliminary")
    paveALICE.AddText("PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("D^{0} #rightarrow K^{-}#pi^{+} and charge conj.")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4")
    paveALICE.AddText("5 < #it{p}_{T,ch jet} < 30 GeV/#it{c}")
    paveALICE.AddText("|#eta_{jet}| < 0.5")
    paveALICE.Draw()

    return canvas


def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    f = open("LHC15i2response_Train1399.yaml", 'r')
    config_b = yaml.load(f)
    f.close()

    f = open("LHC15i2response_Train1399.yaml", 'r')
    config_c = yaml.load(f)
    f.close()

    canvas = EfficiencyComparison(config_c, config_b)
    canvas.SaveAs("{0}/Efficiency_Paper.pdf".format(config_c["input_path"]))
    canvas.SaveAs("{0}/Efficiency_Paper.C".format(config_c["input_path"]))


if __name__ == '__main__':

    main()

    IPython.embed()
