#!/usr/bin/env python
# python script to prepare B feed-down correction file

import argparse
import IPython
import ROOT
import array
import numpy
import math
import yaml
import DMesonJetUnfolding
import DMesonJetUtils
import DMesonJetCompare
import copy
from collections import OrderedDict

globalList = []

def main(config, unfolding_debug):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gSystem.Load("libRooUnfold")

    if "b_response" in config:
        print("Opening the b->D0 response matrix file")
        bResponseFile = OpenResponseFile(config["input_path"], config["b_response"], False)

        print("Opening the b->D0 response matrix file (w/ efficiency)")
        bResponseFile_efficiency = OpenResponseFile(config["input_path"], config["b_response"], True)
    else:
        print("No b->D0 response provided.")
        bResponseFile = None
        bResponseFile_efficiency = None

    if "c_response" in config:
        if config["c_response"]["train"] == config["b_response"]["train"] and config["c_response"]["name"] == config["b_response"]["name"]:
            cResponseFile = bResponseFile
            cResponseFile_efficiency = bResponseFile_efficiency
        else:
            print("Opening the c->D0 response matrix file")
            cResponseFile = OpenResponseFile(config["input_path"], config["c_response"], False)

            print("Opening the c->D0 response matrix file (w/ efficiency)")
            cResponseFile_efficiency = OpenResponseFile(config["input_path"], config["c_response"], True)
    else:
        print("No c->D0 response provided.")
        cResponseFile = None
        cResponseFile_efficiency = None

    results = OrderedDict()
    robjects = []
    for v in config["variations"]:
        if not v["active"]: continue
        name = v["name"]
        suffix = "_".join([config["generator"], str(v["ts"])])
        input_file_name = "{0}/FastSim_{1}/FastSimAnalysis_Reduced_{1}.root".format(config["input_path"], suffix)
        (FDhistogram_jetpt_orig, FDhistogram_jetz_orig, FDhistogram_dpt_orig) = LoadFDHistogram(input_file_name)
        results[name] = OrderedDict()
        results[name].update(PrepareFDhist_dpt(v["ts"], FDhistogram_dpt_orig, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency))
        results[name].update(PrepareFDhist_jetpt(v["ts"], FDhistogram_jetpt_orig, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug))
        results[name].update(PrepareFDhist_jetz(v["ts"], FDhistogram_jetz_orig, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug))
        robjects.append(GenerateRootList(results[name], name))

    results["SystematicUncertainty"] = CompareVariations(config["variations"], results)
    robjects.append(GenerateRootList(results["SystematicUncertainty"], "SystematicUncertainty"))
    PlotFSspectraAndSyst(results)

    output_file_name = "{0}/{1}.root".format(config["input_path"], config["name"])
    outputFile = ROOT.TFile(output_file_name, "recreate")
    if not outputFile or outputFile.IsZombie():
        print("Could not open output file {0}".format(output_file_name))
        exit(1)
    for obj in robjects:
        obj.Write(obj.GetName(), ROOT.TObject.kSingleKey)
    outputFile.Close()

    SaveCanvases(config["name"], config["input_path"])

    print("Done: {0}.".format(output_file_name))

def PlotFSspectraAndSyst(results):
    spectrumNames = ["DPtSpectrum/GeneratorLevel_DPtSpectrum",
                     "DPtSpectrum/GeneratorLevel_DPtSpectrum_bEfficiencyMultiply",
                     "DPtSpectrum/GeneratorLevel_DPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetPtSpectrum_DPt_30/GeneratorLevel_JetPtSpectrum",
                     "JetPtSpectrum_DPt_30/GeneratorLevel_JetPtSpectrum_bEfficiencyMultiply",
                     "JetPtSpectrum_DPt_30/DetectorLevel_JetPtSpectrum_bEfficiencyMultiply",
                     "JetPtSpectrum_DPt_30/GeneratorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetPtSpectrum_DPt_30/DetectorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetZSpectrum_DPt_30_JetPt_5_30/Unfolded_c_JetZSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetZSpectrum_DPt_30_JetPt_5_30/GeneratorLevel_JetZSpectrum",
                     "JetZSpectrum_DPt_30_JetPt_5_30/GeneratorLevel_JetZSpectrum_bEfficiencyMultiply",
                     "JetZSpectrum_DPt_30_JetPt_5_30/DetectorLevel_JetZSpectrum_bEfficiencyMultiply",
                     "JetZSpectrum_DPt_30_JetPt_5_30/GeneratorLevel_JetZSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetZSpectrum_DPt_30_JetPt_5_30/DetectorLevel_JetZSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetZSpectrum_DPt_30_JetPt_5_30/Unfolded_c_JetZSpectrum_bEfficiencyMultiply_cEfficiencyDivide"
                     ]

    for name in spectrumNames:
        stat = GetSpectrum(results["default"], name)
        if not stat: continue
        hname = name[name.rfind("/") + 1:]
        syst = results["SystematicUncertainty"][name]["{0}_CentralAsymmSyst".format(hname)]
        canvas = ROOT.TCanvas("{0}_canvas".format(name), "{0}_canvas".format(name))
        if not "JetZ" in name: canvas.SetLogy()
        canvas.SetLeftMargin(0.13)
        canvas.cd()
        syst_copy = syst.Clone()
        syst_copy.SetFillColor(ROOT.kCyan + 1)
        syst_copy.SetMarkerColor(ROOT.kCyan + 1)
        syst_copy.GetYaxis().SetTitleOffset(1.5)
        syst_copy.Draw("a e2")
        stat_copy = stat.DrawCopy("same p e0 x0")
        stat_copy.Scale(1., "width")
        stat_copy.SetMarkerColor(ROOT.kBlue + 1)
        stat_copy.SetMarkerStyle(ROOT.kFullCircle)
        stat_copy.SetMarkerSize(0.6)
        stat_copy.SetLineColor(ROOT.kBlue + 1)
        if "JetZ" in name:
            leg = ROOT.TLegend(0.15, 0.71, 0.70, 0.89, "NB")
        else:
            leg = ROOT.TLegend(0.35, 0.71, 0.89, 0.89, "NB")
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(20)
        leg.AddEntry(stat_copy, "Central values with stat. unc.", "pe")
        leg.AddEntry(syst_copy, "Systematic uncertainty", "f")
        leg.Draw()
        globalList.append(syst_copy)
        globalList.append(syst_copy.GetHistogram())
        globalList.append(stat_copy)
        globalList.append(leg)
        globalList.append(canvas)

def GetSpectrum(results, name):
    subdirs = name.split("/")
    if len(subdirs) < 1:
        print("{0} is not a valid spectrum name".format(name))
        exit(1)
    for subdir in subdirs:
        if not (isinstance(results, dict) or isinstance(results, OrderedDict)):
            print("ERROR: {0} is not a dictionary!".format(results))
            exit(1)
        if not subdir in results:
            print("Could not find key '{0}' in dictionary '{1}'".format(subdir, results))
            return None
        results = results[subdir]
    return results

def CompareVariationsForSpectrum(comp_template, variations, results, name):
    h = GetSpectrum(results["default"], name)
    if not h: return None
    comp = copy.deepcopy(comp_template)
    comp.fName = name
    baseline = h.Clone("{0}_copy".format(h.GetName()))
    baseline.SetTitle(variations[0]["title"])
    spectra = []
    for v in variations:
        vname = v["name"]
        if not v["active"] or vname == "default": continue
        h = GetSpectrum(results[vname], name)
        h_copy = h.Clone("{0}_copy".format(h.GetName()))
        h_copy.SetTitle(v["title"])
        spectra.append(h_copy)
    comp.CompareSpectra(baseline, spectra)
    globalList.append(baseline)
    globalList.extend(spectra)
    globalList.extend(comp.fResults)
    return GenerateSystematicUncertainty(baseline, spectra)

def GenerateSystematicUncertainty(baseline, spectra):
    hname = baseline.GetName().replace("_copy", "")
    upperLimitsHist = ROOT.TH1D("{0}_UpperSyst".format(hname), "{0}_UpperSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    lowerLimitsHist = ROOT.TH1D("{0}_LowerSyst".format(hname), "{0}_LowerSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    symmetricLimitsHist = ROOT.TH1D("{0}_SymmSyst".format(hname), "{0}_SymmSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    relativeSystHist = ROOT.TH1D("{0}_RelSyst".format(hname), "{0}_RelSyst".format(hname), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    result = OrderedDict()
    result[upperLimitsHist.GetName()] = upperLimitsHist
    result[lowerLimitsHist.GetName()] = lowerLimitsHist
    result[symmetricLimitsHist.GetName()] = symmetricLimitsHist
    result[relativeSystHist.GetName()] = relativeSystHist
    for ibin in range(1, baseline.GetNbinsX() + 1):
        centralValue = baseline.GetBinContent(ibin)
        for var in spectra:
            diff = var.GetBinContent(ibin) - centralValue
            if diff > upperLimitsHist.GetBinContent(ibin):
                upperLimitsHist.SetBinContent(ibin, diff)
                print("Bin {0}, upper limit {1}".format(ibin, diff))
            if -diff > lowerLimitsHist.GetBinContent(ibin):
                lowerLimitsHist.SetBinContent(ibin, -diff)
                print("Bin {0}, lower limit {1}".format(ibin, -diff))
        symmetricLimitsHist.SetBinContent(ibin, max(upperLimitsHist.GetBinContent(ibin), lowerLimitsHist.GetBinContent(ibin)))
        if baseline.GetBinContent(ibin) != 0:
            relativeSystHist.SetBinContent(ibin, symmetricLimitsHist.GetBinContent(ibin) / baseline.GetBinContent(ibin))
    xArray = numpy.array([baseline.GetXaxis().GetBinCenter(ibin) for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    xArrayErr = numpy.array([baseline.GetXaxis().GetBinWidth(ibin) / 2 for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    yArray = numpy.array([baseline.GetBinContent(ibin) / baseline.GetXaxis().GetBinWidth(ibin) for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrUp = numpy.array([upperLimitsHist.GetBinContent(ibin) / upperLimitsHist.GetXaxis().GetBinWidth(ibin) for ibin in range(1, upperLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrLow = numpy.array([lowerLimitsHist.GetBinContent(ibin) / lowerLimitsHist.GetXaxis().GetBinWidth(ibin) for ibin in range(1, lowerLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrSym = numpy.array([symmetricLimitsHist.GetBinContent(ibin) / symmetricLimitsHist.GetXaxis().GetBinWidth(ibin) for ibin in range(1, symmetricLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    symmetricUncGraph = ROOT.TGraphErrors(baseline.GetNbinsX(), xArray, yArray, xArrayErr, yArrayErrSym)
    symmetricUncGraph.SetName("{0}_CentralSymmSyst".format(hname))
    symmetricUncGraph.GetXaxis().SetTitle(baseline.GetXaxis().GetTitle())
    asymmetricUncGraph = ROOT.TGraphAsymmErrors(baseline.GetNbinsX(), xArray, yArray, xArrayErr, xArrayErr, yArrayErrLow, yArrayErrUp)
    asymmetricUncGraph.SetName("{0}_CentralAsymmSyst".format(hname))
    asymmetricUncGraph.GetXaxis().SetTitle(baseline.GetXaxis().GetTitle())
    if "it{p}" in baseline.GetXaxis().GetTitle():
        symmetricUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [(mb) (GeV/#it{c})^{-1}]")
        asymmetricUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [(mb) (GeV/#it{c})^{-1}]")
    elif "it{z}" in baseline.GetXaxis().GetTitle():
        symmetricUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{z}_{||,D}^{ch jet}} (mb)")
        asymmetricUncGraph.GetYaxis().SetTitle("#frac{d#sigma}{d#it{z}_{||,D}^{ch jet}} (mb)")

    result[symmetricUncGraph.GetName()] = symmetricUncGraph
    result[asymmetricUncGraph.GetName()] = asymmetricUncGraph
    return result

def CompareVariations(variations, results):
    comp_template = DMesonJetCompare.DMesonJetCompare("DMesonJetCompare")
    comp_template.fOptRatio = "hist"
    comp_template.fLogUpperSpace = 3
    comp_template.fNColsLegRatio = 2
    comp_template.fX1LegRatio = 0.2

    result = OrderedDict()

    spectrumNames = ["DPtSpectrum/GeneratorLevel_DPtSpectrum",
                     "DPtSpectrum/GeneratorLevel_DPtSpectrum_bEfficiencyMultiply",
                     "DPtSpectrum/GeneratorLevel_DPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetPtSpectrum_DPt_30/GeneratorLevel_JetPtSpectrum",
                     "JetPtSpectrum_DPt_30/GeneratorLevel_JetPtSpectrum_bEfficiencyMultiply",
                     "JetPtSpectrum_DPt_30/DetectorLevel_JetPtSpectrum_bEfficiencyMultiply",
                     "JetPtSpectrum_DPt_30/GeneratorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetPtSpectrum_DPt_30/DetectorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetZSpectrum_DPt_30_JetPt_5_30/Unfolded_c_JetZSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetZSpectrum_DPt_30_JetPt_5_30/GeneratorLevel_JetZSpectrum",
                     "JetZSpectrum_DPt_30_JetPt_5_30/GeneratorLevel_JetZSpectrum_bEfficiencyMultiply",
                     "JetZSpectrum_DPt_30_JetPt_5_30/DetectorLevel_JetZSpectrum_bEfficiencyMultiply",
                     "JetZSpectrum_DPt_30_JetPt_5_30/GeneratorLevel_JetZSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetZSpectrum_DPt_30_JetPt_5_30/DetectorLevel_JetZSpectrum_bEfficiencyMultiply_cEfficiencyDivide",
                     "JetZSpectrum_DPt_30_JetPt_5_30/Unfolded_c_JetZSpectrum_bEfficiencyMultiply_cEfficiencyDivide"
                     ]

    for name in spectrumNames:
        if "JetZ" in name:
            comp_template.fDoSpectraPlot = "lineary"
            comp_template.fX1LegSpectrum = 0.15
            comp_template.fX2LegSpectrum = 0.70
        else:
            comp_template.fDoSpectraPlot = "logy"
            comp_template.fX1LegSpectrum = 0.55
            comp_template.fX2LegSpectrum = 0.90
        r = CompareVariationsForSpectrum(comp_template, variations, results, name)
        if r: result[name] = r

    return result

def SaveCanvases(name, input_path):
    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            oname = obj.GetName().replace("/", "_")
            obj.SaveAs("{0}/{1}_{2}.pdf".format(input_path, name, oname))

def GenerateRootList(pdict, name):
    rlist = ROOT.TList()
    rlist.SetName(name)
    print("Loop in {0}".format(name))
    for nobj, obj in pdict.iteritems():
        if isinstance(obj, ROOT.TObject):
            rlist.Add(obj)
        elif isinstance(obj, dict) or isinstance(obj, OrderedDict):
            print("Recursion in {0}".format(nobj))
            rlist.Add(GenerateRootList(obj, nobj))
        else:
            print("Error: type of object {0} not recognized!".format(obj))
    return rlist

def PrepareFDhist_dpt(ts, FDhistogram_old, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency):
    print("Preparing D pt FD histograms")
    result = OrderedDict()

    dpt = OrderedDict()
    result["DPtSpectrum"] = dpt

    if bResponseFile or cResponseFile:
        responseList = OrderedDict()
        dpt["DetectorResponse"] = responseList

    dptbins = [2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 30]

    print("Rebinning the FD original histogram")
    FDhistogram_orig = FDhistogram_old.Rebin(len(dptbins) - 1, FDhistogram_old.GetName(), array.array('d', dptbins))
    FDhistogram_orig.SetName("GeneratorLevel_DPtSpectrum")
    dpt[FDhistogram_orig.GetName()] = FDhistogram_orig

    if bResponseFile:
        print("Loading the response matrix b->D0")
        temp, bEfficiency_dpt = LoadResponse(bResponseFile, "JetPtDPtSpectrum_NonPrompt", "_NoJet", "b", FDhistogram_orig.GetNbinsX(), FDhistogram_orig.GetXaxis().GetXbins().GetArray())
        bEfficiency_dpt.SetName("EfficiencyVsDPt_b")
        responseList[bEfficiency_dpt.GetName()] = bEfficiency_dpt
    else:
        bEfficiency_dpt = None

    if cResponseFile:
        print("Loading the response matrix c->D0")
        temp, cEfficiency_dpt = LoadResponse(cResponseFile, "JetPtDPtSpectrum_Prompt", "_NoJet", "c", FDhistogram_orig.GetNbinsX(), FDhistogram_orig.GetXaxis().GetXbins().GetArray())
        cEfficiency_dpt.SetName("EfficiencyVsDPt_c")
        responseList[cEfficiency_dpt.GetName()] = cEfficiency_dpt
    else:
        cEfficiency_dpt = None

    if bEfficiency_dpt:
        print("Applying the b->D0 reconstruction efficiency")
        FDhistogram_dpt = FDhistogram_orig.Clone("{0}_bEff".format(FDhistogram_orig.GetName()))
        FDhistogram_dpt.Multiply(bEfficiency_dpt)
        FDhistogram_dpt.SetName("GeneratorLevel_DPtSpectrum_bEfficiencyMultiply")
        dpt[FDhistogram_dpt.GetName()] = FDhistogram_dpt
    else:
        FDhistogram_dpt = None

    if FDhistogram_dpt and cEfficiency_dpt:
        print("Applying the correction for the c->D0 reconstruction efficiency")
        FDhistogram_efficiency_dpt = FDhistogram_dpt.Clone("{0}_efficiency".format(FDhistogram_dpt.GetName()))
        FDhistogram_efficiency_dpt.Divide(cEfficiency_dpt)
        FDhistogram_efficiency_dpt.SetName("GeneratorLevel_DPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")
        dpt[FDhistogram_efficiency_dpt.GetName()] = FDhistogram_efficiency_dpt
    else:
        FDhistogram_efficiency_dpt = None

    return result

def PrepareFDhist_jetpt(ts, FDhistogram_old, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug):
    jetptbins = [5, 6, 8, 10, 14, 20, 30]
    return PrepareFDhist_jet("JetPt", jetptbins, "", ts, FDhistogram_old, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug)

def PrepareFDhist_jetz(ts, FDhistogram_old, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug):
    jetzbins = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    return PrepareFDhist_jet("JetZ", jetzbins, "_JetPt_5_30", ts, FDhistogram_old, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug)

def PrepareFDhist_jet(jet_var_name, jetxbins, additional_kin_cuts, ts, FDhistogram_old, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug):
    print("Preparing {} FD histograms".format(jet_var_name))
    result = OrderedDict()

    spectrum_2d_name = "{}DPtSpectrum".format(jet_var_name)

    jetxdpt = OrderedDict()
    result[spectrum_2d_name] = jetxdpt

    dptbins = [2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 30]

    if bResponseFile or cResponseFile:
        responseList = OrderedDict()
        jetxdpt["DetectorResponse"] = responseList

    print("Rebinning the FD original histogram")
    FDhistogram_orig = DMesonJetUtils.Rebin2D_fromBins(FDhistogram_old, FDhistogram_old.GetName(), len(jetxbins) - 1, array.array('d', jetxbins), len(dptbins) - 1, array.array('d', dptbins))
    FDhistogram_orig.SetName("GeneratorLevel_{}".format(spectrum_2d_name))
    jetxdpt[FDhistogram_orig.GetName()] = FDhistogram_orig

    jet_var_limits = "_{}_{}_{}".format(jet_var_name, int(jetxbins[0] * 100), int(jetxbins[-1] * 100))

    if bResponseFile:
        print("Loading the response matrix b->D0")
        temp, bEfficiency_dpt = LoadResponse(bResponseFile, "{}_NonPrompt".format(spectrum_2d_name), jet_var_limits, "b", FDhistogram_orig.GetNbinsY(), FDhistogram_orig.GetYaxis().GetXbins().GetArray())
        bEfficiency_dpt.SetName("EfficiencyVsDPt_b")
        responseList[bEfficiency_dpt.GetName()] = bEfficiency_dpt
    else:
        bEfficiency_dpt = None

    if cResponseFile:
        print("Loading the response matrix c->D0")
        temp, cEfficiency_dpt = LoadResponse(cResponseFile, "{}_Prompt".format(spectrum_2d_name), jet_var_limits, "c", FDhistogram_orig.GetNbinsY(), FDhistogram_orig.GetYaxis().GetXbins().GetArray())
        cEfficiency_dpt.SetName("EfficiencyVsDPt_c")
        responseList[cEfficiency_dpt.GetName()] = cEfficiency_dpt
    else:
        cEfficiency_dpt = None

    if bEfficiency_dpt:
        print("Applying the b->D0 reconstruction efficiency")
        FDhistogram = FDhistogram_orig.Clone("{0}_bEff".format(FDhistogram_orig.GetName()))
        ApplyEfficiency(FDhistogram, bEfficiency_dpt, False)
        FDhistogram.SetName("GeneratorLevel_{}_bEfficiencyMultiply".format(spectrum_2d_name))
        jetxdpt[FDhistogram.GetName()] = FDhistogram
    else:
       FDhistogram = None

    if FDhistogram and cEfficiency_dpt:
        print("Applying the correction for the c->D0 reconstruction efficiency")
        FDhistogram_efficiency = FDhistogram.Clone("{0}_efficiency".format(FDhistogram.GetName()))
        ApplyEfficiency(FDhistogram_efficiency, cEfficiency_dpt, True)
        FDhistogram_efficiency.SetName("GeneratorLevel_{}_bEfficiencyMultiply_cEfficiencyDivide".format(spectrum_2d_name))
        jetxdpt[FDhistogram_efficiency.GetName()] = FDhistogram_efficiency
    else:
        FDhistogram_efficiency = None

    if bEfficiency_dpt:
        # Finer bin versions used as priors in the folding/unfolding
        print("Applying the b->D0 reconstruction efficiency (fine bins)")
        FDhistogram_fineBins = DMesonJetUtils.Rebin2D_fromBins(FDhistogram_old, "{0}_fineBins_bEff".format(FDhistogram_old.GetName()), FDhistogram_old.GetNbinsX(), FDhistogram_old.GetXaxis().GetXbins().GetArray(), len(dptbins) - 1, array.array('d', dptbins))
        ApplyEfficiency(FDhistogram_fineBins, bEfficiency_dpt, False)
        FDhistogram_fineBins.SetName("GeneratorLevel_{}_FineBins_bEfficiencyMultiply".format(spectrum_2d_name))
        jetxdpt[FDhistogram_fineBins.GetName()] = FDhistogram_fineBins
    else:
        FDhistogram_fineBins = None

    if FDhistogram_fineBins and cEfficiency_dpt:
        print("Applying the correction for the c->D0 reconstruction efficiency (fine bins)")
        FDhistogram_fineBins_efficiency = FDhistogram_fineBins.Clone("{0}_efficiency".format(FDhistogram_fineBins.GetName()))
        ApplyEfficiency(FDhistogram_fineBins_efficiency, cEfficiency_dpt, True)
        FDhistogram_fineBins_efficiency.SetName("GeneratorLevel_{}_FineBins_bEfficiencyMultiply_cEfficiencyDivide".format(spectrum_2d_name))
        jetxdpt[FDhistogram_fineBins_efficiency.GetName()] = FDhistogram_fineBins_efficiency
    else:
        FDhistogram_fineBins_efficiency = None

    for ptd in range(3, 4):
        spectrumName = "{}Spectrum".format(jet_var_name)
        fullSpectrumName = "{}_DPt_{}{}".format(spectrumName, ptd * 10, additional_kin_cuts)
        nonPromptFullSpectrumName = "{}_NonPrompt".format(fullSpectrumName)
        promptFullSpectrumName = "{}_Prompt".format(fullSpectrumName)

        ptdList = OrderedDict()
        result[fullSpectrumName] = ptdList

        print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/o b or c efficiency)".format(ptd))
        FDhistogram_jetx_orig = FDhistogram_orig.ProjectionX("{}_{}_DPt_{}".format(FDhistogram_orig.GetName(), jet_var_name, ptd * 10), FDhistogram_orig.GetYaxis().FindBin(ptd), FDhistogram_orig.GetNbinsY() + 1)
        FDhistogram_jetx_orig.SetName("GeneratorLevel_{}".format(spectrumName))
        FDhistogram_jetx_orig.GetYaxis().SetTitle(FDhistogram_orig.GetZaxis().GetTitle())
        ptdList[FDhistogram_jetx_orig.GetName()] = FDhistogram_jetx_orig

        if FDhistogram:
            print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/o efficiency)".format(ptd))
            FDhistogram_jetx = FDhistogram.ProjectionX("{}_{}_DPt_{}".format(FDhistogram.GetName(), jet_var_name, ptd * 10), FDhistogram.GetYaxis().FindBin(ptd), FDhistogram.GetNbinsY() + 1)
            FDhistogram_jetx.SetName("GeneratorLevel_{}_bEfficiencyMultiply".format(spectrumName))
            FDhistogram_jetx.GetYaxis().SetTitle(FDhistogram.GetZaxis().GetTitle())
            ptdList[FDhistogram_jetx.GetName()] = FDhistogram_jetx
        else:
            FDhistogram_jetx = None

        if FDhistogram_fineBins:
            print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/o efficiency, fine bins)".format(ptd))
            FDhistogram_fineBins_jetx = FDhistogram_fineBins.ProjectionX("{}_{}_DPt_{}".format(FDhistogram_fineBins.GetName(), jet_var_name, ptd * 10), FDhistogram_fineBins.GetYaxis().FindBin(ptd), FDhistogram_fineBins.GetNbinsY() + 1)
            FDhistogram_fineBins_jetx.SetName("GeneratorLevel_{}_FineBins_bEfficiencyMultiply".format(spectrumName))
            FDhistogram_fineBins_jetx.GetYaxis().SetTitle(FDhistogram_fineBins.GetZaxis().GetTitle())
            ptdList[FDhistogram_fineBins_jetx.GetName()] = FDhistogram_fineBins_jetx
        else:
            FDhistogram_fineBins_jetx = None

        if bResponseFile:
            print("Applying the jet b response matrix")
            FDhistogram_jetx_detector = ApplyResponse(FDhistogram_jetx, bResponseFile, FDhistogram_fineBins_jetx, nonPromptFullSpectrumName)
            FDhistogram_jetx_detector.SetName("DetectorLevel_{}_bEfficiencyMultiply".format(spectrumName))
            ptdList[FDhistogram_jetx_detector.GetName()] = FDhistogram_jetx_detector
        else:
            FDhistogram_jetx_detector = None

        if cResponseFile and FDhistogram_jetx_detector:
            print("Unfolding using the c response matrix")
            FDhistogram_jetx_unfolded = GetUnfoldedSpectrum("{0}_{1}_unfolded_c".format(ts, FDhistogram_jetx_detector.GetName()), FDhistogram_jetx_detector, cResponseFile, FDhistogram_fineBins_jetx, promptFullSpectrumName, unfolding_debug)
            FDhistogram_jetx_unfolded.SetName("Unfolded_c_{}_bEfficiencyMultiply".format(spectrumName))
            ptdList[FDhistogram_jetx_unfolded.GetName()] = FDhistogram_jetx_unfolded
        else:
            FDhistogram_jetx_unfolded = None

        if FDhistogram_jetx_detector and bResponseFile:
            print("Unfolding using the b response matrix")
            FDhistogram_jetx_unfolded = GetUnfoldedSpectrum("{0}_{1}_unfolded_b".format(ts, FDhistogram_jetx_detector.GetName()), FDhistogram_jetx_detector, bResponseFile, FDhistogram_fineBins_jetx, nonPromptFullSpectrumName, unfolding_debug)
            FDhistogram_jetx_unfolded.SetName("Unfolded_b_{}_bEfficiencyMultiply".format(spectrumName))
            ptdList[FDhistogram_jetx_unfolded.GetName()] = FDhistogram_jetx_unfolded
        else:
            FDhistogram_jetx_unfolded = None

        if FDhistogram_efficiency:
            print("Projecting FD histogram into the jet axis for D pt > {0} GeV/c (w/ efficiency)".format(ptd))
            FDhistogram_efficiency_jetx = FDhistogram_efficiency.ProjectionX("{}_{}_DPt_{}".format(FDhistogram_efficiency.GetName(), jet_var_name, ptd * 10), FDhistogram_efficiency.GetYaxis().FindBin(ptd), FDhistogram_efficiency.GetNbinsY() + 1)
            FDhistogram_efficiency_jetx.SetName("GeneratorLevel_{}_bEfficiencyMultiply_cEfficiencyDivide".format(spectrumName))
            FDhistogram_efficiency_jetx.GetYaxis().SetTitle(FDhistogram_efficiency.GetZaxis().GetTitle())
            ptdList[FDhistogram_efficiency_jetx.GetName()] = FDhistogram_efficiency_jetx
        else:
            FDhistogram_efficiency_jetx = None

        if FDhistogram_fineBins_efficiency:
            print("Projecting FD histogram into the jet axis for D pt > {0} GeV/c (w/ efficiency, fine bins)".format(ptd))
            FDhistogram_fineBins_efficiency_jetx = FDhistogram_fineBins_efficiency.ProjectionX("{}_{}_DPt_{}".format(FDhistogram_fineBins_efficiency.GetName(), jet_var_name, ptd * 10), FDhistogram_fineBins_efficiency.GetYaxis().FindBin(ptd), FDhistogram_fineBins_efficiency.GetNbinsY() + 1)
            FDhistogram_fineBins_efficiency_jetx.SetName("GeneratorLevel_{}_FineBins_bEfficiencyMultiply_cEfficiencyDivide".format(spectrumName))
            FDhistogram_fineBins_efficiency_jetx.GetYaxis().SetTitle(FDhistogram_fineBins_efficiency.GetZaxis().GetTitle())
            ptdList[FDhistogram_fineBins_efficiency_jetx.GetName()] = FDhistogram_fineBins_efficiency_jetx
        else:
            FDhistogram_fineBins_efficiency_jetx = None

        if FDhistogram_efficiency_jetx and bResponseFile_efficiency:
            print("Applying the jet b response matrix")
            FDhistogram_efficiency_jetx_detector = ApplyResponse(FDhistogram_efficiency_jetx, bResponseFile_efficiency, FDhistogram_fineBins_efficiency_jetx, nonPromptFullSpectrumName)
            FDhistogram_efficiency_jetx_detector.SetName("DetectorLevel_{}_bEfficiencyMultiply_cEfficiencyDivide".format(spectrumName))
            ptdList[FDhistogram_efficiency_jetx_detector.GetName()] = FDhistogram_efficiency_jetx_detector
        else:
            FDhistogram_efficiency_jetx_detector = None

        if FDhistogram_efficiency_jetx_detector and cResponseFile_efficiency:
            print("Unfolding using the c response matrix")
            FDhistogram_efficiency_jetx_unfolded = GetUnfoldedSpectrum("{0}_{1}_unfolded_c".format(ts, FDhistogram_efficiency_jetx_detector.GetName()), FDhistogram_efficiency_jetx_detector, cResponseFile_efficiency, FDhistogram_fineBins_efficiency_jetx, promptFullSpectrumName, unfolding_debug)
            FDhistogram_efficiency_jetx_unfolded.SetName("Unfolded_c_{}_bEfficiencyMultiply_cEfficiencyDivide".format(spectrumName))
            ptdList[FDhistogram_efficiency_jetx_unfolded.GetName()] = FDhistogram_efficiency_jetx_unfolded
        else:
            FDhistogram_efficiency_jetx_unfolded = None

        if FDhistogram_efficiency_jetx_detector and bResponseFile_efficiency:
            print("Unfolding using the b response matrix")
            FDhistogram_efficiency_jetx_unfolded_b = GetUnfoldedSpectrum("{0}_{1}_unfolded_b".format(ts, FDhistogram_efficiency_jetx_detector.GetName()), FDhistogram_efficiency_jetx_detector, bResponseFile_efficiency, FDhistogram_fineBins_efficiency_jetx, nonPromptFullSpectrumName, unfolding_debug)
            FDhistogram_efficiency_jetx_unfolded_b.SetName("Unfolded_b_{}_bEfficiencyMultiply_cEfficiencyDivide".format(spectrumName))
            ptdList[FDhistogram_efficiency_jetx_unfolded_b.GetName()] = FDhistogram_efficiency_jetx_unfolded_b
        else:
            FDhistogram_efficiency_jetx_unfolded_b = None

    return result

def GenerateUnfoldingEngine(name, responseFile, prior, spectrumName, nbins):
    analysis = dict()
    analysis["name"] = name
    analysis["d_meson"] = "_d_meson_"
    analysis["variable"] = "_variable_"
    analysis["jet_type"] = "Charged"
    analysis["jet_radius"] = "R040"
    analysis["active"] = True
    analysis["kinematic_cuts"] = None
    analysis["raw_yield_method"] = None
    analysis["d_meson_response"] = "D0"
    analysis["spectrum_response"] = spectrumName
    analysis["priors"] = ["GeneratedSpectrum"]
    analysis["default_prior"] = "GeneratedSpectrum"
    analysis["default_method"] = "Svd"
    svd = dict()
    default_reg_svd = {"GeneratedSpectrum" : nbins}
    svd["default_reg"] = default_reg_svd
    bayes = dict()
    default_reg_bayes = {"GeneratedSpectrum" : 6}
    bayes["default_reg"] = default_reg_bayes
    bayes["iter_min"] = 1
    bayes["iter_max"] = 6
    bayes["iter_step"] = 1
    binBybin = dict()
    binBybin["default_reg"] = None
    analysis["methods"] = {"Svd" : svd, "Bayes" : bayes, "BinByBin" : binBybin}
    unf = DMesonJetUnfolding.DMesonJetUnfoldingEngine(None, None, "unfolding", analysis)
    unf.fDoErrorCrossChecks = False
    unf.SetCustomPrior("GeneratedSpectrum", prior)
    unf.fUseOverflow = True
    unf.fResponseFile = responseFile
    unf.LoadResponse(False)
    return unf

def GetUnfoldedSpectrum(name, histogram, responseFile, prior, spectrumName, unfolding_debug):
    unf = GenerateUnfoldingEngine(name, responseFile, prior, spectrumName, histogram.GetXaxis().GetNbins())
    unf.fInputSpectrum = histogram
    unf.fTruthSpectrum = None
    unf.Start(unfolding_debug)
    default_reg = unf.GetDefaultRegularization(unf.fDefaultMethod, unf.fDefaultPrior)
    default_unfolded = unf.fUnfoldedSpectra[unf.fDefaultMethod, default_reg, unf.fDefaultPrior]
    default_unfolded.SetName(histogram.GetName().replace("_detector", "_unfolded"))
    return default_unfolded

def ApplyEfficiency(hist, efficiency, reverse):
    print("Applying efficiency {0} to histogram {1}".format(efficiency.GetName(), hist.GetName()))
    for ybin in range(0, hist.GetNbinsY() + 2):
        if reverse and efficiency.GetBinContent(ybin) > 0:
            print("(bin {0:.2f} ({1:.2f}): dividing by {2}".format(hist.GetYaxis().GetBinUpEdge(ybin), efficiency.GetXaxis().GetBinUpEdge(ybin), efficiency.GetBinContent(ybin)))
            eff = 1. / efficiency.GetBinContent(ybin)
        else:
            print("(bin {0:.2f} ({1:.2f}): multiplying by {2}".format(hist.GetYaxis().GetBinUpEdge(ybin), efficiency.GetXaxis().GetBinUpEdge(ybin), efficiency.GetBinContent(ybin)))
            eff = efficiency.GetBinContent(ybin)
        for xbin in range(0, hist.GetNbinsX() + 2):
            hist.SetBinContent(xbin, ybin, hist.GetBinContent(xbin, ybin) * eff)
            hist.SetBinError(xbin, ybin, hist.GetBinError(xbin, ybin) * eff)

def ApplyResponse(truth, responseFile, prior, spectrumName):
    unf = GenerateUnfoldingEngine(truth.GetName(), responseFile, prior, spectrumName, truth.GetXaxis().GetNbins())
    unf.fInputSpectrum = truth
    unf.fTruthSpectrum = None
    unf.GenerateResponse()
    resp = unf.fResponseMatrices["GeneratedSpectrum"].fNormResponse.Clone("temp")
    for ybin in range(0, resp.GetYaxis().GetNbins() + 2):
        inty = resp.Integral(0, resp.GetXaxis().GetNbins() + 1, ybin, ybin)
        if inty == 0:
            continue
        print("Integral is {0:.3f}".format(inty))
        scaling = truth.GetBinContent(ybin) / inty
        for xbin in range(0, resp.GetXaxis().GetNbins() + 2):
            binValue = resp.GetBinContent(xbin, ybin) * scaling
            binErr = resp.GetBinError(xbin, ybin) * scaling
            resp.SetBinContent(xbin, ybin, binValue)
            resp.SetBinError(xbin, ybin, binErr)
    result = resp.ProjectionX("{0}_detector".format(truth.GetName()), 0, -1)
    result.GetYaxis().SetTitle(truth.GetYaxis().GetTitle())
    return result

def OpenResponseFile(input_path, response, efficiency):
    if efficiency:
        file_name = "{0}/{1}/{2}_efficiency.root".format(input_path, response["train"], response["name"])
    else:
        file_name = "{0}/{1}/{2}.root".format(input_path, response["train"], response["name"])
    file = ROOT.TFile(file_name)
    print("Response file: {0}".format(file_name))
    if not file or file.IsZombie():
        print("Could not open input file {0}".format(file_name))
        exit(1)
    return file

def LoadResponse(responseFile, spectrumName, suffix_in, suffix_out, detAxisNbins, detAxisBin):
    rlistName = "D0_Jet_AKTChargedR040_pt_scheme_{0}".format(spectrumName)
    rlist = responseFile.Get(rlistName)
    if not rlist:
        print("Could not get list {0} from file {1}".format(rlistName, responseFile.GetName()))
        exit(1)
    truthName = "D0_Jet_AKTChargedR040_pt_scheme_{0}_Truth{1}".format(spectrumName, suffix_in)
    truth = rlist.FindObject(truthName)
    if not truth:
        print("Could not get histogram {0}".format(truthName))
        exit(1)
    truth.SetName("{0}_{1}".format(truthName, suffix_out))
    truth_coarse = truth.Rebin(detAxisNbins, "{0}_coarse_{1}".format(truthName, suffix_out), detAxisBin)
    respName = "D0_Jet_AKTChargedR040_pt_scheme_{0}_DetectorResponse{1}".format(spectrumName, suffix_in)
    resp = rlist.FindObject(respName)
    if not resp:
        print("Could not get histogram {0}".format(respName))
        exit(1)
    resp.SetName("{0}_{1}".format(respName, suffix_out))
    resp_coarse = DMesonJetUtils.Rebin2D_fromBins(resp, "{0}_coarse_{1}".format(respName, suffix_out), detAxisNbins, detAxisBin, detAxisNbins, detAxisBin)
    eff_coarse = resp_coarse.ProjectionY(truth_coarse.GetName().replace("Truth", "Efficiency"))
    eff_coarse.Divide(truth_coarse)
    return resp_coarse, eff_coarse

def LoadFDHistogram(file_name):
    result = []
    for jet_var_name, label in zip(["JetPt", "JetZ"], ["#it{p}_{T}", "#it{z}"]):
        spectrum_name = "{}DPtSpectrum".format(jet_var_name)
        file = ROOT.TFile(file_name)
        if not file or file.IsZombie():
            print("Could not open input file {0}".format(file_name))
            exit(1)
        else:
            print("File {0} open".format(file_name))
        dlist = file.Get("D0_MCTruth")
        jlist = dlist.FindObject("Charged_R040")
        slist = jlist.FindObject("D0_MCTruth_Charged_R040_{}".format(spectrum_name))
        if not slist:
            print("Could not find FD histogram ({})!".format(jet_var_name))
            jlist.Print()
            exit(1)
        jet_hist = slist.FindObject("D0_MCTruth_Charged_R040_{}".format(spectrum_name))
        if not jet_hist:
            print("Could not find FD histogram ({})!".format(jet_var_name))
            slist.Print()
            exit(1)
        jet_hist.GetZaxis().SetTitle("#frac{{d#sigma}}{{d{lab}}} #times #Delta{lab} (mb)".format(lab=label))
        print("Histogram {} loaded".format(jet_hist))
        result.append(jet_hist)

    slist = dlist.FindObject("D0_MCTruth_DPtSpectrum")
    if not slist:
        print("Could not find FD histogram (d pt)!")
        dlist.ls()
        exit(1)
    dpt_hist = slist.FindObject("D0_MCTruth_DPtSpectrum")
    if not dpt_hist:
        print("Could not find FD histogram (d pt)!")
        slist.ls()
        exit(1)
    dpt_hist.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} #times #Delta#it{p}_{T} (mb)")
    result.append(dpt_hist)

    return result

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare B feed-down correction file.')
    parser.add_argument('yaml', metavar='conf.yaml')
    parser.add_argument("--unfolding-debug", action='store_const',
                        default=False, const=True,
                        help='Unfolding debug plots.')

    args = parser.parse_args()

    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.unfolding_debug)

    IPython.embed()
