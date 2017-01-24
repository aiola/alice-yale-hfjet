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
    dmesonName = "D0"
    prefix = "{0}_{1}_{2}".format(dmesonName, jetName, spectrumName)

    DPtBins = [3, 4, 5, 6, 7, 8, 10, 12, 16, 30]

    minJetPt = 5
    maxJetPt = 30
    recoTruthName = "{0}/{0}_ReconstructedTruth_JetPt_{1}_{2}".format(prefix, minJetPt * 100, maxJetPt * 100)
    truthName = "{0}/{0}_Truth_JetPt_{1}_{2}".format(prefix, minJetPt * 100, maxJetPt * 100)
    h = DMesonJetUtils.GetObject(file_c, recoTruthName)
    cRecoTruth = h.Rebin(len(DPtBins) - 1, "{0}_c_rebin".format(recoTruthName), array.array('d', DPtBins))
    h = DMesonJetUtils.GetObject(file_c, truthName)
    cTruth = h.Rebin(len(DPtBins) - 1, "{0}_c_rebin".format(truthName), array.array('d', DPtBins))
    h = DMesonJetUtils.GetObject(file_b, recoTruthName)
    bRecoTruth = h.Rebin(len(DPtBins) - 1, "{0}_b_rebin".format(recoTruthName), array.array('d', DPtBins))
    h = DMesonJetUtils.GetObject(file_b, truthName)
    bTruth = h.Rebin(len(DPtBins) - 1, "{0}_b_rebin".format(truthName), array.array('d', DPtBins))
    c_hist_Ratio = cRecoTruth.Clone(recoTruthName.replace("RecontructedTruth", "Efficiency"))
    c_hist_Ratio.Divide(cTruth)
    c_hist_Ratio.SetTitle("D^{0}, prompt")
    b_hist_Ratio = bRecoTruth.Clone(recoTruthName.replace("RecontructedTruth", "Efficiency"))
    b_hist_Ratio.Divide(bTruth)
    b_hist_Ratio.SetTitle("D^{0}, non-prompt")
    globalList.append(c_hist_Ratio)
    globalList.append(b_hist_Ratio)

    cname = "Efficiency_QM17"
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fOptSpectrum = "hist"
    comp.fOptSpectrumBaseline = "hist"
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = None
    comp.fX1LegSpectrum = 0.16
    comp.fY1LegSpectrum = 0.73
    comp.fColors = [ROOT.kRed + 2, ROOT.kBlue + 2]
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

    h.GetYaxis().SetTitle("Reconstruction Efficiency #times Acceptance")
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

    paveALICE = ROOT.TPaveText(0.14, 0.79, 0.53, 0.95, "NB NDC")
    globalList.append(paveALICE)
    paveALICE.SetBorderSize(0)
    paveALICE.SetFillStyle(0)
    paveALICE.SetTextFont(43)
    paveALICE.SetTextSize(20)
    paveALICE.SetTextAlign(13)
    paveALICE.AddText("ALICE Simulation, PYTHIA6, pp, #sqrt{#it{s}} = 7 TeV")
    paveALICE.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4, |#eta_{jet}| < 0.5")
    paveALICE.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c., #it{p}_{T,D} > 3 GeV/#it{c}")
    paveALICE.Draw()

    return canvas

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    f = open("LHC15i2analysis_Train973.yaml", 'r')
    config_b = yaml.load(f)
    f.close()

    f = open("LHC15i2analysis_Train961.yaml", 'r')
    config_c = yaml.load(f)
    f.close()

    canvas = EfficiencyComparison(config_c, config_b)
    canvas.SaveAs("{0}/Efficiency_QM17.pdf".format(config_c["input_path"]))

if __name__ == '__main__':

    main()

    IPython.embed()
