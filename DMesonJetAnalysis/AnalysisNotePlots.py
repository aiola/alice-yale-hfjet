#!/usr/bin/env python
# python script to generate plots for the analysis note

import argparse
import yaml
import IPython
import ROOT
import DMesonJetUtils
import array

globalList = []
canvases = []

def main(actions, output_path, output_type):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    configs = dict()

    f = open("LHC14j4analysis_Train969.yaml", 'r')
    configs["LHC14j4_cb"] = yaml.load(f)
    f.close()

    f = open("LHC14j4analysis_Train969_efficiency.yaml", 'r')
    configs["LHC14j4_cb_eff"] = yaml.load(f)
    f.close()

    f = open("LHC14j4analysis_Train982.yaml", 'r')
    configs["LHC14j4_c_eff"] = yaml.load(f)
    f.close()

    f = open("LHC14j4analysis_Train982_efficiency.yaml", 'r')
    configs["LHC14j4_c_eff"] = yaml.load(f)
    f.close()

    f = open("LHC15i2analysis_Train961.yaml", 'r')
    configs["LHC15i2_c"] = yaml.load(f)
    f.close()

    f = open("LHC15i2analysis_Train961_efficiency.yaml", 'r')
    configs["LHC15i2_c_eff"] = yaml.load(f)
    f.close()

    f = open("LHC15i2analysis_Train973.yaml", 'r')
    configs["LHC15i2_b"] = yaml.load(f)
    f.close()

    f = open("LHC15i2analysis_Train973_efficiency.yaml", 'r')
    configs["LHC15i2_b_eff"] = yaml.load(f)
    f.close()

    f = open("LHC10analysis_Train823.yaml", 'r')
    configs["LHC10"] = yaml.load(f)
    f.close()

    f = open("LHC10analysis_Train823_efficiency.yaml", 'r')
    configs["LHC10_eff"] = yaml.load(f)
    f.close()

    f = open("LHC10_Train823_LHC15i2_Train961_efficiency.yaml", 'r')
    configs["data_unfolding"] = yaml.load(f)
    f.close()

    f = open("LHC14j4_Train982_LHC15i2_Train961_efficiency.yaml", 'r')
    configs["mc_unfolding"] = yaml.load(f)
    f.close()

    histograms = dict()

    for name, c in configs.iteritems():
        file = OpenFile(c)
        histograms[name] = ExtractRootList(file)

    file = ROOT.TFile("/Volumes/DATA/ALICE/JetResults/FastSim_powheg_beauty_1478869008/stage_1/output/BFeedDown_powheg_beauty_1478869008.root")
    histograms["BFeedDown"] = ExtractRootList(file)

    if "all" in actions or "jetptres_comp" in actions:
        JetPtResolutionComparison(histograms["LHC15i2_c_eff"], histograms["LHC15i2_b_eff"])

    if "all" in actions or "eff_comp" in actions:
        EfficiencyComparison(histograms["LHC15i2_c"], histograms["LHC15i2_b"])

    if "all" in actions or "fd_fold_unfold" in actions:
        FD_FoldUnfold_Comparison(histograms["BFeedDown"])

    for c in canvases:
        c.SaveAs("{0}/{1}.{2}".format(output_path, c.GetName(), output_type))

    if "data":
        CopyDataFilesWithoutEff(configs["LHC10"])
        CopyDataFilesWithEff(configs["LHC10_eff"])

    if "mc":
        CopyMCFilesWithoutEff(configs["LHC14j4_cb"])
        CopyMCFilesWithEff(configs["LHC14j4_cb_eff"])
        CopyMCFilesWithoutEffCharmOnly(configs["LHC14j4_c"])
        CopyMCFilesWithEffCharmOnly(configs["LHC14j4_c_eff"])

    if "data_unfolding":
        CopyDataUnfoldingFiles(configs["data_unfolding"])

    if "mc_unfolding":
        CopyMCUnfoldingFiles(configs["mc_unfolding"])


def CopyDataFilesWithoutEff(config):
    pass

def CopyDataFilesWithEff(config):
    pass

def CopyMCFilesWithoutEff(config):
    pass

def CopyMCFilesWithEff(config):
    pass

def CopyMCFilesWithoutEffCharmOnly(config):
    pass

def CopyMCFilesWithEffCharmOnly(config):
    pass

def CopyDataUnfoldingFiles(config):
    pass

def CopyMCUnfoldingFiles(config):
    pass

def FD_FoldUnfold_Comparison(histograms):
    baseline = histograms["D0_MCTruth_Charged_R040_Jet_Pt_D_Pt_Spectrum_bEff_efficiency_jetpt_DPt_20"].Clone()
    detector = histograms["D0_MCTruth_Charged_R040_Jet_Pt_D_Pt_Spectrum_bEff_efficiency_jetpt_DPt_20_detector"].Clone()
    unfolded = histograms["D0_MCTruth_Charged_R040_Jet_Pt_D_Pt_Spectrum_bEff_efficiency_jetpt_DPt_20_unfolded"].Clone()
    unfolded_b = histograms["D0_MCTruth_Charged_R040_Jet_Pt_D_Pt_Spectrum_bEff_efficiency_jetpt_DPt_20_unfolded_b"].Clone()
    baseline.SetTitle("POWHEG+PYTHIA FD spectrum")
    detector.SetTitle("Smeared w/ b #rightarrow D^{0} detector response")
    unfolded.SetTitle("Unfolded w/ c #rightarrow D^{0} detector response")
    unfolded_b.SetTitle("Unfolded w/ b #rightarrow D^{0} detector response")
    hist = [detector, unfolded, unfolded_b]
    globalList.append(baseline)
    globalList.extend(hist)
    r = DMesonJetUtils.CompareSpectra(baseline, hist, "FD_FoldUnfold_Comparison")
    for obj in r:
        globalList.append(obj)
        if isinstance(obj, ROOT.TCanvas):
            canvases.append(obj)

def EfficiencyComparison(hist_c, hist_b):
    spectrumName = "JetPtDPtSpectrum"
    jetName = "Jet_AKTChargedR040_pt_scheme"
    dmesonName = "D0"
    prefix = "{0}_{1}_{2}".format(dmesonName, jetName, spectrumName)
    jetPtLimits = [5, 30]
    cname = "ReconstructionEfficiencyComparison"
    opt = "hist"
    c = None
    leg = None
    detResp_c = hist_c[prefix]
    detResp_b = hist_b[prefix]
    DPtBins = [2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 30]
    colors = [ROOT.kBlue + 2, ROOT.kGreen + 2, ROOT.kRed + 2, ROOT.kMagenta + 2, ROOT.kCyan + 2, ROOT.kOrange + 2]
    markers = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenDiamond, ROOT.kOpenStar]
    lines = [1, 2, 9, 5, 7, 10]
    styles = dict()
    for i, (minPt, maxPt) in enumerate(zip(jetPtLimits[:-1], jetPtLimits[1:])):
        styles["colors"] = colors[i * 2:(i + 1) * 2]
        styles["lines"] = lines[i * 2:(i + 1) * 2]
        styles["markers"] = markers[i * 2:(i + 1) * 2]
        recoTruthName = "{0}_RecontructedTruth_JetPt_{1}_{2}".format(prefix, minPt * 100, maxPt * 100)
        truthName = "{0}_Truth_JetPt_{1}_{2}".format(prefix, minPt * 100, maxPt * 100)
        cRecoTruth = detResp_c[recoTruthName].Rebin(len(DPtBins) - 1, "{0}_c_rebin".format(recoTruthName), array.array('d', DPtBins))
        cTruth = detResp_c[truthName].Rebin(len(DPtBins) - 1, "{0}_c_rebin".format(truthName), array.array('d', DPtBins))
        bRecoTruth = detResp_b[recoTruthName].Rebin(len(DPtBins) - 1, "{0}_b_rebin".format(recoTruthName), array.array('d', DPtBins))
        bTruth = detResp_b[truthName].Rebin(len(DPtBins) - 1, "{0}_b_rebin".format(truthName), array.array('d', DPtBins))
        c_hist_Ratio = cRecoTruth.Clone(recoTruthName.replace("RecontructedTruth", "Efficiency"))
        c_hist_Ratio.Divide(cTruth)
        c_hist_Ratio.SetTitle("c #rightarrow D^{{0}}, {0}".format(c_hist_Ratio.GetTitle()))
        c_hist_Ratio.GetYaxis().SetTitle("Reconstruction Efficiency")
        c_hist_Ratio.GetXaxis().SetTitle("#it{p}_{T,D} (GeV/#it{c})")
        b_hist_Ratio = bRecoTruth.Clone(recoTruthName.replace("RecontructedTruth", "Efficiency"))
        b_hist_Ratio.Divide(bTruth)
        b_hist_Ratio.SetTitle("b #rightarrow D^{{0}}, {0}".format(b_hist_Ratio.GetTitle()))
        globalList.append(c_hist_Ratio)
        globalList.append(b_hist_Ratio)
        r = DMesonJetUtils.CompareSpectra(c_hist_Ratio, [b_hist_Ratio], cname, opt, "", "", "lineary", None, c, None, leg, None, styles)
        opt += "same"
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TLegend):
                leg = obj
            elif isinstance(obj, ROOT.TCanvas):
                c = obj
    globalList.append(c)
    canvases.append(c)

def JetPtResolutionComparison(hist_c, hist_b):
    spectrumName = "JetPtSpectrum_DPt_20"
    jetName = "Jet_AKTChargedR040_pt_scheme"
    dmesonName = "D0"
    listName = "{0}_{1}_{2}".format(dmesonName, jetName, spectrumName)
    detResp_c = hist_c[listName]["DetectorResponse"]
    detResp_b = hist_b[listName]["DetectorResponse"]

    canvas = DMesonJetUtils.GenerateMultiCanvas("DetectorJetPtResolutionComparison", len(detResp_c))
    globalList.append(canvas)
    canvases.append(canvas)
    for i, (c, b) in enumerate(zip(detResp_c.itervalues(), detResp_b.itervalues())):
        pad = canvas.cd(i + 1)
        # pad.SetLogy()
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.05)
        pad.SetTopMargin(0.08)
        pad.SetBottomMargin(0.18)
        h = c.DrawCopy()
        h.GetYaxis().SetTitle("Probability density")
        h.GetXaxis().SetRangeUser(-1, 0.5)
        h.SetMarkerColor(ROOT.kRed + 2)
        h.SetMarkerStyle(ROOT.kFullCircle)
        h.SetMarkerSize(0.8)
        h.SetLineColor(ROOT.kRed + 2)
        globalList.append(h)
        h.GetXaxis().SetTitleFont(43)
        h.GetXaxis().SetTitleOffset(2.8)
        h.GetXaxis().SetTitleSize(19)
        h.GetXaxis().SetLabelFont(43)
        h.GetXaxis().SetLabelOffset(0.009)
        h.GetXaxis().SetLabelSize(18)
        h.GetYaxis().SetTitleFont(43)
        h.GetYaxis().SetTitleOffset(2.6)
        h.GetYaxis().SetTitleSize(19)
        h.GetYaxis().SetLabelFont(43)
        h.GetYaxis().SetLabelOffset(0.009)
        h.GetYaxis().SetLabelSize(18)
        # h.SetMaximum(h.GetMaximum()*1.3)
        h.SetMinimum(0)
        htitle = ROOT.TPaveText(0.12, 0.99, 0.95, 0.93, "NB NDC")
        htitle.SetBorderSize(0)
        htitle.SetFillStyle(0)
        htitle.SetTextFont(43)
        htitle.SetTextSize(18)
        htitle.AddText(h.GetTitle())
        htitle.Draw()
        globalList.append(htitle)
        hbis = b.DrawCopy("same")
        hbis.SetMarkerColor(ROOT.kBlue + 2)
        hbis.SetMarkerStyle(ROOT.kOpenCircle)
        hbis.SetMarkerSize(0.9)
        hbis.SetLineColor(ROOT.kBlue + 2)
    canvas.cd(1)
    leg = ROOT.TLegend(0.19, 0.77, 0.58, 0.95, "NB NDC")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(24)
    leg.SetTextAlign(13)
    leg.AddEntry(h, "c #rightarrow D^{0}", "pe")
    leg.AddEntry(hbis, "b #rightarrow D^{0}", "pe")
    leg.Draw()
    globalList.append(leg)

def OpenFile(config):
    file = ROOT.TFile.Open("{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"]))
    return file

def ExtractRootList(input):
    result = dict()
    if isinstance(input, ROOT.TDirectory):
        keys = input.GetListOfKeys()
        for k in keys:
            obj = input.Get(k.GetName())
            if isinstance(obj, ROOT.TH1) or isinstance(obj, ROOT.THnBase):
                result[k.GetName()] = obj
            elif isinstance(obj, ROOT.TDirectory) or isinstance(obj, ROOT.TCollection):
                result[k.GetName()] = ExtractRootList(obj)
            else:
                print("Object {0} not recognized".format(obj.GetName()))
                print(obj)
    else:
        for h in input:
            if isinstance(h, ROOT.TH1) or isinstance(h, ROOT.THnBase):
                result[h.GetName()] = h
            elif isinstance(h, ROOT.TDirectory) or isinstance(h, ROOT.TCollection):
                result[h.GetName()] = ExtractRootList(h)
            else:
                print("Object {0} not recognized".format(h.GetName()))
                print(h)
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Side Band analysis.')
    parser.add_argument('actions', metavar='action',
                        help='Actions to be taken', nargs='*')
    parser.add_argument('-o', metavar='path',
                        help='Output path', default='../../DJetAnalysisNote/pp_plots')
    parser.add_argument('-f', metavar='format',
                        help='Format (pdf, eps, png,...)', default='pdf')
    args = parser.parse_args()

    main(args.actions, args.o, args.f)

    IPython.embed()
