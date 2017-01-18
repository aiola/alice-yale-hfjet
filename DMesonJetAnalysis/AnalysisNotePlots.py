#!/usr/bin/env python
# python script to generate plots for the analysis note

import argparse
import yaml
import IPython
import ROOT
import DMesonJetUtils
import DMesonJetCompare
import array
import shutil
import os

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
    configs["LHC14j4_c"] = yaml.load(f)
    f.close()

    f = open("LHC14j4analysis_Train982_efficiency.yaml", 'r')
    configs["LHC14j4_c_eff"] = yaml.load(f)
    f.close()

    f = open("LHC15i2analysis_Train961.yaml", 'r')
    configs["LHC15i2_c"] = yaml.load(f)
    f.close()

    f = open("LHC15i2analysis_Train961_mcshape.yaml", 'r')
    configs["LHC15i2_c_mcshape"] = yaml.load(f)
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
        if "LHC" in name:
            file = OpenFile(c)
            histograms[name] = ExtractRootList(file)

    file = ROOT.TFile("/Volumes/DATA/ALICE/JetResults/BFeedDown.root")
    histograms["BFeedDown"] = ExtractRootList(file)

    if "all" in actions or "jetptres_comp" in actions:
        JetPtResolutionComparison(histograms["LHC15i2_c_eff"], histograms["LHC15i2_b_eff"])

    if "all" in actions or "eff_comp" in actions:
        cname = "BFeedDown/ReconstructionEfficiencyPromptNonPromptComparison"
        title_c = "c #rightarrow D^{0}"
        title_b = "b #rightarrow D^{0}"
        EfficiencyComparison(cname, title_c, histograms["LHC15i2_c"], title_b, histograms["LHC15i2_b"])

    if "all" in actions or "mcshape" in actions:
        cname = "DataSystematics/ReconstructionEfficiencyMCShape"
        title_c = "PYTHIA6"
        title_mcshape = "POWHEG+PYTHIA6"
        EfficiencyComparison(cname, title_c, histograms["LHC15i2_c"], title_mcshape, histograms["LHC15i2_c_mcshape"])
        CopyMCShapeSystematics("/Volumes/DATA/ALICE/JetResults", output_path, output_type)

    if "all" in actions or "fd_fold_unfold" in actions:
        FD_FoldUnfold_Comparison(histograms["BFeedDown"]["default"]["JetPtSpectrum_DPt_30"])

    if "all" in actions or "fit_params" in actions:
        FitParameterComparison(histograms["LHC14j4_cb_eff"], histograms["LHC10_eff"])

    BFeedDownPath = "{0}/BFeedDown".format(output_path)
    if not os.path.isdir(BFeedDownPath):
        os.makedirs(BFeedDownPath)

    for c in canvases:
        c.SaveAs("{0}/{1}.{2}".format(output_path, c.GetName(), output_type))

    if "all" in actions or "data" in actions:
        CopyDataFilesWithoutEff(configs["LHC10"], output_path, output_type)
        CopyDataFilesWithEff(configs["LHC10_eff"], output_path, output_type)

    if "all" in actions or "mc" in actions:
        CopyMCFilesWithoutEff(configs["LHC14j4_cb"], output_path, output_type)
        CopyMCFilesWithEff(configs["LHC14j4_cb_eff"], output_path, output_type)
        CopyMCFilesWithoutEffCharmOnly(configs["LHC14j4_c"], output_path, output_type)
        CopyMCFilesWithEffCharmOnly(configs["LHC14j4_c_eff"], output_path, output_type)

    if "all" in actions or "data_unfolding" in actions:
        CopyDataUnfoldingFiles(configs["data_unfolding"], output_path, output_type)

    if "all" in actions or "mc_unfolding" in actions:
        CopyMCUnfoldingFiles(configs["mc_unfolding"], output_path, output_type)

    if "all" in actions or "efficiency_c" in actions:
        CopyEfficiencyFiles(configs["LHC15i2_c"], output_path, output_type)

    if "all" in actions or "response_c" in actions:
        CopyResponseFiles(configs["LHC15i2_c_eff"], output_path, output_type)

    if "all" in actions or "bfeed_down" in actions:
        CopyBFeedDown("/Volumes/DATA/ALICE/JetResults", output_path, output_type)

    if "all" in actions or "promptd_sim" in actions:
        CopyPromptDSimulation("/Volumes/DATA/ALICE/JetResults", output_path, output_type)

    if "all" in actions or "data_systematics" in actions:
        CopDataSystematics("/Volumes/DATA/ALICE/JetResults", output_path, output_type)

    if "all" in actions or "theory_comparison" in actions:
        CopyTheoryComparisonFiles(configs["data_unfolding"], output_path, output_type)

    if "all" in actions or "ppb_comparison" in actions:
        CopypPbComparisonFiles(configs["data_unfolding"], output_path, output_type)

def CopyFiles(input_path, output_path, file_list, output_type):
    if not os.path.isdir(output_path):
        os.makedirs(output_path)
    for file_name in file_list:
        print("Copying {0}...".format(file_name))
        shutil.copy("{0}/{1}.{2}".format(input_path, file_name, output_type), output_path)

def CopyTheoryComparisonFiles(config, output_path, output_type):
    full_output_path = "{0}/TheoryComparison".format(output_path)
    file_list = []
    file_list.append("D0JetCrossSection_pp7TeV")
    CopyFiles(config["input_path"], full_output_path, file_list, output_type)

def CopyMCShapeSystematics(input_path, output_path, output_type):
    full_output_path = "{0}/DataSystematics".format(output_path)
    file_list = []
    file_list.append("PYTHIA_POWHEG_DPtSpectrumComparison_Charged_R040_LHC15i2analysis_Train961_charm_1483386026")
    file_list.append("PYTHIA_POWHEG_DPtSpectrumComparison_Charged_R040_LHC15i2analysis_Train961_charm_1483386026_Ratio")
    file_list.append("PYTHIA_POWHEG_JetPtSpectrumComparison_Charged_R040_LHC15i2analysis_Train961_charm_1483386026")
    file_list.append("PYTHIA_POWHEG_JetPtSpectrumComparison_Charged_R040_LHC15i2analysis_Train961_charm_1483386026_Ratio")
    CopyFiles(input_path, full_output_path, file_list, output_type)

def CopypPbComparisonFiles(config, output_path, output_type):
    full_output_path = "{0}/pPbComparison".format(output_path)
    file_list = []
    file_list.append("pPbComparison_Charged_R040")
    file_list.append("pPbComparison_Charged_R040_Ratio")
    CopyFiles(config["input_path"], full_output_path, file_list, output_type)

def CopyBFeedDown(input_path, output_path, output_type):
    full_output_path = "{0}/BFeedDown".format(output_path)
    file_list = []
    file_list.append("BFeedDownVsPtD_powheg_Charged_R040")
    file_list.append("BFeedDownVsPtD_powheg_Charged_R040_Ratio")
    file_list.append("BFeedDownVsPtJet_powheg_Charged_R040")
    file_list.append("BFeedDownVsPtJet_powheg_Charged_R040_Ratio")
    file_list.append("BFeedDownVsZ_powheg_Charged_R040")
    file_list.append("BFeedDownVsZ_powheg_Charged_R040_Ratio")
    file_list.append("BFeedDown_DPtSpectrum_GeneratorLevel_DPtSpectrum")
    file_list.append("BFeedDown_DPtSpectrum_GeneratorLevel_DPtSpectrum_Ratio")
    file_list.append("BFeedDown_DPtSpectrum_GeneratorLevel_DPtSpectrum_canvas")
    file_list.append("BFeedDown_JetPtSpectrum_DPt_30_GeneratorLevel_JetPtSpectrum")
    file_list.append("BFeedDown_JetPtSpectrum_DPt_30_GeneratorLevel_JetPtSpectrum_Ratio")
    file_list.append("BFeedDown_JetPtSpectrum_DPt_30_GeneratorLevel_JetPtSpectrum_canvas")
    CopyFiles(input_path, full_output_path, file_list, output_type)

def CopDataSystematics(input_path, output_path, output_type):
    full_output_path = "{0}/DataSystematics".format(output_path)
    file_list = []
    file_list.append("DataSystematics_LHC10")
    file_list.append("CompareVariations_DataSystematics_LHC10")
    file_list.append("CompareVariations_DataSystematics_LHC10_Ratio")
    file_list.append("CompareUncertainties_DataSystematics_LHC10")
    CopyFiles(input_path, full_output_path, file_list, output_type)

def CopyPromptDSimulation(input_path, output_path, output_type):
    full_output_path = "{0}/PromptDSim".format(output_path)
    file_list = []
    file_list.append("PromptDJetsPrediction_DPtSpectrum_GeneratorLevel_DPtSpectrum")
    file_list.append("PromptDJetsPrediction_DPtSpectrum_GeneratorLevel_DPtSpectrum_Ratio")
    file_list.append("PromptDJetsPrediction_DPtSpectrum_GeneratorLevel_DPtSpectrum_canvas")
    file_list.append("PromptDJetsPrediction_JetPtSpectrum_DPt_30_GeneratorLevel_JetPtSpectrum")
    file_list.append("PromptDJetsPrediction_JetPtSpectrum_DPt_30_GeneratorLevel_JetPtSpectrum_Ratio")
    file_list.append("PromptDJetsPrediction_JetPtSpectrum_DPt_30_GeneratorLevel_JetPtSpectrum_canvas")
    CopyFiles(input_path, full_output_path, file_list, output_type)

def CopyDataFilesWithoutEff(config, output_path, output_type):
    full_input_path = "{0}/{1}/{2}/{3}".format(config["input_path"], config["train"], config["name"], output_type)
    full_output_path = "{0}/RawYieldExtractionWithoutEff".format(output_path)
    file_list = []
    file_list.append("D0_Charged_R040_JetPtBins_DPt_30")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_InvMassFit_Normalized_canvas")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_InvMassFit_Unc_canvas")
    file_list.append("D0_Charged_R040_DPtBins_JetPt_5_30_SideBand_D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_Normalized_canvas")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_Unc_canvas")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_BkgVsSig")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_TotalBkgVsSig")
    file_list.append("D0_Charged_R040_jet_pt_SpectraComparison")
    file_list.append("D0_Charged_R040_jet_pt_SpectraComparison_Ratio")
    CopyFiles(full_input_path, full_output_path, file_list, output_type)

def CopyDataFilesWithEff(config, output_path, output_type):
    full_input_path = "{0}/{1}/{2}/{3}".format(config["input_path"], config["train"], config["name"], output_type)
    full_output_path = "{0}/RawYieldExtractionWithEff".format(output_path)
    file_list = []
    file_list.append("D0_Charged_R040_JetPtBins_DPt_30")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_InvMassFit_Normalized_canvas")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_InvMassFit_Unc_canvas")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_Normalized_canvas")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_Unc_canvas")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_TotalBkgVsSig")
    file_list.append("D0_Charged_R040_jet_pt_SpectraComparison")
    file_list.append("D0_Charged_R040_jet_pt_SpectraComparison_Ratio")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_FDCorrection")
    file_list.append("D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand_FDCorrection_Ratio")
    CopyFiles(full_input_path, full_output_path, file_list, output_type)

    full_input_path = "{0}/{1}/{2}/RawYieldUnc/{3}".format(config["input_path"], config["train"], config["name"], output_type)
    full_output_path = "{0}/RawYieldExtractionWithEff".format(output_path)
    file_list = []
    file_list.append("AverageRawYieldVsDefault")
    file_list.append("AverageRawYieldVsDefault_Ratio")
    CopyFiles(full_input_path, full_output_path, file_list, output_type)

def CopyEfficiencyFiles(config, output_path, output_type):
    full_input_path = "{0}/{1}/{2}/{3}".format(config["input_path"], config["train"], config["name"], output_type)
    full_output_path = "{0}/RecoEfficiency".format(output_path)
    file_list = []
    file_list.append("D0_Jet_AKTChargedR040_pt_scheme_JetPtDPtSpectrum_PartialEfficiency")
    file_list.append("D0_Jet_AKTChargedR040_pt_scheme_JetPtDPtSpectrum_PartialEfficiencyRatios")
    CopyFiles(full_input_path, full_output_path, file_list, output_type)

def CopyResponseFiles(config, output_path, output_type):
    full_input_path = "{0}/{1}/{2}/{3}".format(config["input_path"], config["train"], config["name"], output_type)
    full_output_path = "{0}/DetectorResponse".format(output_path)
    file_list = []
    file_list.append("D0_Jet_AKTChargedR040_pt_scheme_JetPtSpectrum_DPt_30_DetectorResponse")
    file_list.append("D0_Jet_AKTChargedR040_pt_scheme_JetPtSpectrum_DPt_30_DetectorResponse_Resolution_canvas")
    file_list.append("D0_Jet_AKTChargedR040_pt_scheme_JetPtSpectrum_DPt_30_DetectorResponse_EnergyScaleShift_canvas")
    file_list.append("D0_Jet_AKTChargedR040_pt_scheme_JetPtSpectrum_DPt_30_DetectorResponse_canvas")
    CopyFiles(full_input_path, full_output_path, file_list, output_type)

def CopyMCFilesWithoutEff(config, output_path, output_type):
    pass

def CopyMCFilesWithEff(config, output_path, output_type):
    full_input_path = "{0}/{1}/{2}/{3}".format(config["input_path"], config["train"], config["name"], output_type)
    full_output_path = "{0}/MCRawYieldExtractionWithEff".format(output_path)
    file_list = []
    file_list.append("D0_Charged_R040_JetPtBins_DPt_30")
    file_list.append("D0_kSignalOnly_Charged_R040_JetPtBins_DPt_30")
    file_list.append("D0_Charged_R040_DPtBins_JetPt_5_30_SideBand_D0_Charged_R040_JetPtSpectrum_DPt_30_SideBand")
    file_list.append("D0_kSignalOnly_Charged_R040_DPtBins_JetPt_5_30")
    file_list.append("Charged_R040_JetPtSpectrum_DPt_30_SpectraComparison")
    file_list.append("Charged_R040_JetPtSpectrum_DPt_30_SpectraComparison_Ratio")
    file_list.append("D0_kSignalOnly_Charged_R040_jet_pt_SpectraComparison")
    file_list.append("D0_kSignalOnly_Charged_R040_jet_pt_SpectraComparison_Ratio")
    file_list.append("D0_MCTruth_Charged_R040_jet_pt_SpectraComparison")
    file_list.append("D0_MCTruth_Charged_R040_jet_pt_SpectraComparison_Ratio")
    CopyFiles(full_input_path, full_output_path, file_list, output_type)

def CopyMCFilesWithoutEffCharmOnly(config, output_path, output_type):
    pass

def CopyMCFilesWithEffCharmOnly(config, output_path, output_type):
    pass

def CopyDataUnfoldingFiles(config, output_path, output_type):
    full_input_path = "{0}/{1}_mt/{2}".format(config["input_path"], config["name"], output_type)
    full_output_path = "{0}/DataUnfolding".format(output_path)
    file_list = []
    file_list.append("SideBand_DPt_30_Response_PriorResponseTruth")
    file_list.append("SideBand_DPt_30_UnfoldingSummary_Bayes")
    file_list.append("SideBand_DPt_30_UnfoldingSummary_Bayes_UnfoldedOverMeasured")
    file_list.append("SideBand_DPt_30_BinByBinCorrectionFactors")
    file_list.append("SideBand_DPt_30_Priors")
    file_list.append("SideBand_DPt_30_Priors_Ratio")
    file_list.append("SideBand_DPt_30_UnfoldingPrior_Bayes")
    file_list.append("SideBand_DPt_30_UnfoldingPrior_Bayes_Ratio")
    file_list.append("SideBand_DPt_30_UnfoldingRegularization_Bayes_PriorResponseTruth")
    file_list.append("SideBand_DPt_30_UnfoldingRegularization_Bayes_PriorResponseTruth_Ratio")
    file_list.append("SideBand_DPt_30_UnfoldingMethod")
    file_list.append("SideBand_DPt_30_UnfoldingMethod_Ratio")
    CopyFiles(full_input_path, full_output_path, file_list, output_type)

def CopyMCUnfoldingFiles(config, output_path, output_type):
    full_input_path = "{0}/{1}/{2}".format(config["input_path"], config["name"], output_type)
    full_output_path = "{0}/MCUnfolding".format(output_path)
    file_list = []
    file_list.append("SideBand_DPt_30_UnfoldingPrior_Bayes_Ratio")
    file_list.append("SideBand_DPt_30_UnfoldingRegularization_Bayes_PriorResponseTruth_Ratio")
    file_list.append("SideBand_DPt_30_UnfoldingMethod_Ratio")
    file_list.append("SignalOnly_DPt_30_UnfoldingPrior_Bayes_Ratio")
    file_list.append("SignalOnly_DPt_30_UnfoldingRegularization_Bayes_PriorResponseTruth_Ratio")
    file_list.append("SignalOnly_DPt_30_UnfoldingMethod_Ratio")
    CopyFiles(full_input_path, full_output_path, file_list, output_type)

def FD_FoldUnfold_Comparison(histograms):
    baseline = histograms["GeneratorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide"].Clone()
    detector = histograms["DetectorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide"].Clone()
    unfolded = histograms["Unfolded_c_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide"].Clone()
    unfolded_b = histograms["Unfolded_b_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide"].Clone()
    baseline.SetTitle("POWHEG+PYTHIA FD spectrum")
    detector.SetTitle("Smeared w/ b #rightarrow D^{0} detector response")
    unfolded.SetTitle("Unfolded w/ c #rightarrow D^{0} detector response")
    unfolded_b.SetTitle("Unfolded w/ b #rightarrow D^{0} detector response")
    hist = [detector, unfolded, unfolded_b]
    globalList.append(baseline)
    globalList.extend(hist)
    comp = DMesonJetCompare.DMesonJetCompare("BFeedDown/FD_FoldUnfold_Comparison")
    comp.fX1LegRatio = 0.25
    comp.fX1LegSpectrum = 0.25
    r = comp.CompareSpectra(baseline, hist)
    for obj in r:
        globalList.append(obj)
        if isinstance(obj, ROOT.TCanvas):
            canvases.append(obj)

def EfficiencyComparison(cname, title_c, hist_c, title_b, hist_b):
    spectrumName = "JetPtDPtSpectrum"
    jetName = "Jet_AKTChargedR040_pt_scheme"
    dmesonName = "D0"
    prefix = "{0}_{1}_{2}".format(dmesonName, jetName, spectrumName)
    jetPtLimits = [5, 30]
    detResp_c = hist_c[prefix]
    detResp_b = hist_b[prefix]
    DPtBins = [3, 4, 5, 6, 7, 8, 10, 12, 16, 30]
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fOptSpectrum = "hist"
    comp.fDoSpectraPlot = "lineary"
    comp.fDoRatioPlot = None
    comp.fX1LegSpectrum = 0.35
    colors = [ROOT.kBlue + 2, ROOT.kGreen + 2, ROOT.kRed + 2, ROOT.kMagenta + 2, ROOT.kCyan + 2, ROOT.kOrange + 2]
    markers = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenDiamond, ROOT.kOpenStar]
    lines = [1, 2, 9, 5, 7, 10]
    for i, (minPt, maxPt) in enumerate(zip(jetPtLimits[:-1], jetPtLimits[1:])):
        comp.fColors = colors[i * 2:(i + 1) * 2]
        comp.fLines = lines[i * 2:(i + 1) * 2]
        comp.fMarkers = markers[i * 2:(i + 1) * 2]
        recoTruthName = "{0}_ReconstructedTruth_JetPt_{1}_{2}".format(prefix, minPt * 100, maxPt * 100)
        truthName = "{0}_Truth_JetPt_{1}_{2}".format(prefix, minPt * 100, maxPt * 100)
        cRecoTruth = detResp_c[recoTruthName].Rebin(len(DPtBins) - 1, "{0}_c_rebin".format(recoTruthName), array.array('d', DPtBins))
        cTruth = detResp_c[truthName].Rebin(len(DPtBins) - 1, "{0}_c_rebin".format(truthName), array.array('d', DPtBins))
        bRecoTruth = detResp_b[recoTruthName].Rebin(len(DPtBins) - 1, "{0}_b_rebin".format(recoTruthName), array.array('d', DPtBins))
        bTruth = detResp_b[truthName].Rebin(len(DPtBins) - 1, "{0}_b_rebin".format(truthName), array.array('d', DPtBins))
        c_hist_Ratio = cRecoTruth.Clone(recoTruthName.replace("RecontructedTruth", "Efficiency"))
        c_hist_Ratio.Divide(cTruth)
        c_hist_Ratio.SetTitle("{0}, {1}".format(title_c, c_hist_Ratio.GetTitle()))
        c_hist_Ratio.GetYaxis().SetTitle("Reconstruction Efficiency")
        c_hist_Ratio.GetXaxis().SetTitle("#it{p}_{T,D} (GeV/#it{c})")
        b_hist_Ratio = bRecoTruth.Clone(recoTruthName.replace("RecontructedTruth", "Efficiency"))
        b_hist_Ratio.Divide(bTruth)
        b_hist_Ratio.SetTitle("{0}, {1}".format(title_b, b_hist_Ratio.GetTitle()))
        globalList.append(c_hist_Ratio)
        globalList.append(b_hist_Ratio)
        r = comp.CompareSpectra(c_hist_Ratio, [b_hist_Ratio])
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TLegend):
                leg = obj
            elif isinstance(obj, ROOT.TCanvas):
                c = obj
    globalList.append(c)
    canvases.append(c)

def JetPtResolutionComparison(hist_c, hist_b):
    spectrumName = "JetPtSpectrum_DPt_30"
    jetName = "Jet_AKTChargedR040_pt_scheme"
    dmesonName = "D0"
    listName = "{0}_{1}_{2}".format(dmesonName, jetName, spectrumName)
    detResp_c = hist_c[listName]["DetectorResponse"]
    detResp_b = hist_b[listName]["DetectorResponse"]

    canvas = DMesonJetUtils.GenerateMultiCanvas("BFeedDown/DetectorJetPtResolutionPromptNonPromptComparison", len(detResp_c))
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

def FitParameterComparison(MC_histos, data_histos):
    spectrumName = ["JetPtSpectrum_DPt_30", "DPtSpectrum_JetPt_5_30"]
    jetName = "Charged_R040"
    variables = ["Mass", "MassWidth"]
    DoFitParameterComparison(MC_histos, data_histos, spectrumName, jetName, variables)

def DoFitParameterComparison(MC_histos, data_histos, spectrumName, jetName, variables):
    for var in variables:
        for spectrum in spectrumName:
            cname = "MCRawYieldExtractionWithEff/{0}_{1}_{2}_Comparison".format(var, jetName, spectrum)
            hname = "D0_kSignalOnly_{0}_{1}_SignalOnly".format(jetName, spectrum)
            baseline = MC_histos["D0_kSignalOnly"][jetName][hname]["{0}_{1}".format(hname, var)]
            baseline.SetTitle("MC w/o background")
            globalList.append(baseline)

            histos = []
            hname = "D0_{0}_{1}_InvMassFit".format(jetName, spectrum)
            h = MC_histos["D0"][jetName][hname]["{0}_{1}".format(hname, var)]
            h.SetTitle("MC w/ background")
            histos.append(h)
            globalList.append(h)
            h = data_histos["D0"][jetName][hname]["{0}_{1}".format(hname, var)]
            h.SetTitle("Data")
            histos.append(h)
            globalList.append(h)
            comp = DMesonJetCompare.DMesonJetCompare(cname)
            comp.fDoSpectraPlot = "lineary"
            comp.fDoRatioPlot = None
            r = comp.CompareSpectra(baseline, histos)

            for obj in r:
                globalList.append(obj)
                if isinstance(obj, ROOT.TCanvas):
                    canvases.append(obj)

def OpenFile(config):
    file = ROOT.TFile.Open("{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"]))
    return file

def ExtractRootList(input):
    result = dict()
    if isinstance(input, ROOT.TDirectory):
        keys = input.GetListOfKeys()
        for k in keys:
            obj = input.Get(k.GetName())
            if isinstance(obj, ROOT.TH1) or isinstance(obj, ROOT.THnBase) or isinstance(obj, ROOT.TGraph):
                result[k.GetName()] = obj
            elif isinstance(obj, ROOT.TDirectory) or isinstance(obj, ROOT.TCollection):
                result[k.GetName()] = ExtractRootList(obj)
            else:
                print("Object {0} not recognized".format(obj.GetName()))
                print(obj)
    else:
        for h in input:
            if isinstance(h, ROOT.TH1) or isinstance(h, ROOT.THnBase) or isinstance(h, ROOT.TGraph):
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
