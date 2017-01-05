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

    print("Opening the b->D0 response matrix file")
    bResponseFile = OpenResponseFile(config["input_path"], config["b_response"], False)

    print("Opening the c->D0 response matrix file")
    cResponseFile = OpenResponseFile(config["input_path"], config["c_response"], False)

    print("Opening the b->D0 response matrix file (w/ efficiency)")
    bResponseFile_efficiency = OpenResponseFile(config["input_path"], config["b_response"], True)

    print("Opening the c->D0 response matrix file (w/ efficiency)")
    cResponseFile_efficiency = OpenResponseFile(config["input_path"], config["c_response"], True)

    results = OrderedDict()
    robjects = []
    for v in config["variations"]:
        if not v["active"]: continue
        name = v["name"]
        suffix = "_".join([config["generator"], str(v["ts"])])
        input_file_name = "{0}/FastSim_{1}/stage_{2}/output/FastSimAnalysis_{1}.root".format(config["input_path"], suffix, v["stage"])
        (FDhistogram_jetpt_orig, FDhistogram_dpt_orig) = LoadFDHistogram(input_file_name)
        results[name] = OrderedDict()
        results[name].update(PrepareFDhist_dpt(v["ts"], FDhistogram_dpt_orig, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency))
        results[name].update(PrepareFDhist_jetpt(v["ts"], FDhistogram_jetpt_orig, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug))
        robjects.append(GenerateRootList(results[name], name))

    results["SystematicUncertainty"] = CompareVariations(config["variations"], results)
    robjects.append(GenerateRootList(results["SystematicUncertainty"], "SystematicUncertainty"))

    output_file_name = "{0}/{1}.root".format(config["input_path"], config["name"])
    outputFile = ROOT.TFile(output_file_name, "recreate")
    if not outputFile or outputFile.IsZombie():
        print("Could not open output file {0}".format(output_file_name))
        exit(1)
    for obj in robjects:
        obj.Write(obj.GetName(), ROOT.TObject.kSingleKey)
    outputFile.Close()

    SaveCanvases(config["input_path"])

    print("Done: {0}.".format(output_file_name))

def GetSpectrum(results, name):
    subdirs = name.split("/")
    if len(subdirs) < 1:
        print("{0} is not a valid spectrum name".format(name))
        exit(1)
    for subdir in subdirs:
        if not (isinstance(results, dict) or isinstance(results, OrderedDict)):
            print("ERROR: {0} is not a dictionary!".format(results))
            exit(1)
        results = results[subdir]
    return results

def CompareVariationsForSpectrum(comp_template, variations, results, name, spectrum):
    comp = copy.deepcopy(comp_template)
    comp.fName = name
    h = GetSpectrum(results["default"], spectrum)
    baseline = h.Clone("{0}_copy".format(h.GetName()))
    baseline.SetTitle(variations[0]["title"])
    spectra = []
    for v in variations:
        name = v["name"]
        if not v["active"] or name == "default": continue
        h = GetSpectrum(results[name], spectrum)
        h_copy = h.Clone("{0}_copy".format(h.GetName()))
        h_copy.SetTitle(v["title"])
        spectra.append(h_copy)
    comp.CompareSpectra(baseline, spectra)
    globalList.append(baseline)
    globalList.extend(spectra)
    globalList.extend(comp.fResults)
    return GenerateSystematicUncertainty(baseline, spectra)

def GenerateSystematicUncertainty(baseline, spectra):
    upperLimitsHist = ROOT.TH1D("{0}_UpperSyst".format(baseline.GetName()), "{0}_UpperSyst".format(baseline.GetName()), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    lowerLimitsHist = ROOT.TH1D("{0}_LowerSyst".format(baseline.GetName()), "{0}_LowerSyst".format(baseline.GetName()), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    symmetricLimitsHist = ROOT.TH1D("{0}_SymmSyst".format(baseline.GetName()), "{0}_SymmSyst".format(baseline.GetName()), baseline.GetNbinsX(), baseline.GetXaxis().GetXbins().GetArray())
    result = OrderedDict()
    result[upperLimitsHist.GetName()] = upperLimitsHist
    result[lowerLimitsHist.GetName()] = lowerLimitsHist
    result[symmetricLimitsHist.GetName()] = symmetricLimitsHist
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
    xArray = numpy.array([baseline.GetXaxis().GetBinCenter(ibin) for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    xArrayErr = numpy.array([baseline.GetXaxis().GetBinWidth(ibin) / 2 for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    yArray = numpy.array([baseline.GetBinContent(ibin) / baseline.GetXaxis().GetBinWidth(ibin) for ibin in range(1, baseline.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrUp = numpy.array([upperLimitsHist.GetBinContent(ibin) / upperLimitsHist.GetXaxis().GetBinWidth(ibin) for ibin in range(1, upperLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrLow = numpy.array([lowerLimitsHist.GetBinContent(ibin) / lowerLimitsHist.GetXaxis().GetBinWidth(ibin) for ibin in range(1, lowerLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    yArrayErrSym = numpy.array([symmetricLimitsHist.GetBinContent(ibin) / symmetricLimitsHist.GetXaxis().GetBinWidth(ibin) for ibin in range(1, symmetricLimitsHist.GetNbinsX() + 1)], dtype=numpy.float32)
    symmetricUncGraph = ROOT.TGraphErrors(baseline.GetNbinsX(), xArray, yArray, xArrayErr, yArrayErrSym)
    symmetricUncGraph.SetName("{0}_CentralSymmSyst".format(baseline.GetName()))
    asymmetricUncGraph = ROOT.TGraphAsymmErrors(baseline.GetNbinsX(), xArray, yArray, xArrayErr, xArrayErr, yArrayErrLow, yArrayErrUp)
    asymmetricUncGraph.SetName("{0}_CentralAsymmSyst".format(baseline.GetName()))
    result[symmetricUncGraph.GetName()] = symmetricUncGraph
    result[asymmetricUncGraph.GetName()] = asymmetricUncGraph
    return result

def CompareVariations(variations, results):
    comp_template = DMesonJetCompare.DMesonJetCompare("DMesonJetCompare")
    comp_template.fOptRatio = "hist"
    comp_template.fLogUpperSpace = 3

    result = OrderedDict()

    name = "GeneratorLevel_DPtSpectrum"
    result[name] = CompareVariationsForSpectrum(comp_template, variations, results, name,
                                                "DPtSpectrum/GeneratorLevel_DPtSpectrum")

    name = "GeneratorLevel_JetPtSpectrum_DPt_30_bEfficiencyMultiply_cEfficiencyDivide"
    result[name] = CompareVariationsForSpectrum(comp_template, variations, results, name,
                                                "JetPtSpectrum_DPt_30/GeneratorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")

    name = "DetectorLevel_JetPtSpectrum_DPt_30_bEfficiencyMultiply_cEfficiencyDivide"
    result[name] = CompareVariationsForSpectrum(comp_template, variations, results, name,
                                                "JetPtSpectrum_DPt_30/DetectorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")

    name = "Unfolded_c_JetPtSpectrum_DPt_30_bEfficiencyMultiply_cEfficiencyDivide"
    result[name] = CompareVariationsForSpectrum(comp_template, variations, results, name,
                                                "JetPtSpectrum_DPt_30/Unfolded_c_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")
    return result

def SaveCanvases(input_path):
    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{0}/BFeedDown_{1}.pdf".format(input_path, obj.GetName()))

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

    dptbins = [2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 24, 32]

    responseList = OrderedDict()
    dpt["DetectorResponse"] = responseList

    print("Rebinning the FD original histogram")
    FDhistogram_orig = FDhistogram_old.Rebin(len(dptbins) - 1, FDhistogram_old.GetName(), array.array('d', dptbins))
    FDhistogram_orig.SetName("GeneratorLevel_DPtSpectrum")
    dpt[FDhistogram_orig.GetName()] = FDhistogram_orig

    print("Loading the response matrix b->D0")
    temp, bEfficiency_dpt = LoadResponse(bResponseFile, "JetPtDPtSpectrum_CoarseBins", "_NoJet", "b", FDhistogram_orig.GetNbinsX(), FDhistogram_orig.GetXaxis().GetXbins().GetArray())
    bEfficiency_dpt.SetName("EfficiencyVsDPt_b")
    responseList[bEfficiency_dpt.GetName()] = bEfficiency_dpt

    print("Loading the response matrix c->D0")
    temp, cEfficiency_dpt = LoadResponse(cResponseFile, "JetPtDPtSpectrum_CoarseBins", "_NoJet", "c", FDhistogram_orig.GetNbinsX(), FDhistogram_orig.GetXaxis().GetXbins().GetArray())
    cEfficiency_dpt.SetName("EfficiencyVsDPt_c")
    responseList[cEfficiency_dpt.GetName()] = cEfficiency_dpt

    print("Applying the b->D0 reconstruction efficiency")
    FDhistogram_dpt = FDhistogram_orig.Clone("{0}_bEff".format(FDhistogram_orig.GetName()))
    FDhistogram_dpt.Multiply(bEfficiency_dpt)
    FDhistogram_dpt.SetName("GeneratorLevel_DPtSpectrum_bEfficiencyMultiply")
    dpt[FDhistogram_dpt.GetName()] = FDhistogram_dpt

    print("Applying the correction for the c->D0 reconstruction efficiency")
    FDhistogram_efficiency_dpt = FDhistogram_dpt.Clone("{0}_efficiency".format(FDhistogram_dpt.GetName()))
    FDhistogram_efficiency_dpt.Divide(cEfficiency_dpt)
    FDhistogram_efficiency_dpt.SetName("GeneratorLevel_DPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")
    dpt[FDhistogram_efficiency_dpt.GetName()] = FDhistogram_efficiency_dpt

    return result

def PrepareFDhist_jetpt(ts, FDhistogram_old, bResponseFile, cResponseFile, bResponseFile_efficiency, cResponseFile_efficiency, unfolding_debug):
    print("Preparing jet pt FD histograms")
    result = OrderedDict()

    jetptdpt = OrderedDict()
    result["JetPtDPtSpectrum"] = jetptdpt

    dptbins = [2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 24, 32]
    jetptbins = [5, 6, 8, 10, 14, 20, 32]

    responseList = OrderedDict()
    jetptdpt["DetectorResponse"] = responseList

    print("Rebinning the FD original histogram")
    FDhistogram_orig = DMesonJetUtils.Rebin2D_fromBins(FDhistogram_old, FDhistogram_old.GetName(), len(jetptbins) - 1, array.array('d', jetptbins), len(dptbins) - 1, array.array('d', dptbins))
    FDhistogram_orig.SetName("GeneratorLevel_JetPtDPtSpectrum")
    jetptdpt[FDhistogram_orig.GetName()] = FDhistogram_orig

    print("Loading the response matrix b->D0")
    temp, bEfficiency_dpt = LoadResponse(bResponseFile, "JetPtDPtSpectrum_FineBins", "_JetPt_500_3200", "b", FDhistogram_orig.GetNbinsY(), FDhistogram_orig.GetYaxis().GetXbins().GetArray())
    bEfficiency_dpt.SetName("EfficiencyVsDPt_b")
    responseList[bEfficiency_dpt.GetName()] = bEfficiency_dpt

    print("Loading the response matrix c->D0")
    temp, cEfficiency_dpt = LoadResponse(cResponseFile, "JetPtDPtSpectrum_FineBins", "_JetPt_500_3200", "c", FDhistogram_orig.GetNbinsY(), FDhistogram_orig.GetYaxis().GetXbins().GetArray())
    cEfficiency_dpt.SetName("EfficiencyVsDPt_c")
    responseList[cEfficiency_dpt.GetName()] = cEfficiency_dpt

    print("Applying the b->D0 reconstruction efficiency")
    FDhistogram = FDhistogram_orig.Clone("{0}_bEff".format(FDhistogram_orig.GetName()))
    ApplyEfficiency(FDhistogram, bEfficiency_dpt, False)
    FDhistogram.SetName("GeneratorLevel_JetPtDPtSpectrum_bEfficiencyMultiply")
    jetptdpt[FDhistogram.GetName()] = FDhistogram

    print("Applying the correction for the c->D0 reconstruction efficiency")
    FDhistogram_efficiency = FDhistogram.Clone("{0}_efficiency".format(FDhistogram.GetName()))
    ApplyEfficiency(FDhistogram_efficiency, cEfficiency_dpt, True)
    FDhistogram_efficiency.SetName("GeneratorLevel_JetPtDPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")
    jetptdpt[FDhistogram_efficiency.GetName()] = FDhistogram_efficiency

    # Finer bin versions used as priors in the folding/unfolding
    print("Applying the b->D0 reconstruction efficiency (fine bins)")
    FDhistogram_fineBins = DMesonJetUtils.Rebin2D_fromBins(FDhistogram_old, "{0}_fineBins_bEff".format(FDhistogram_old.GetName()), FDhistogram_old.GetNbinsX(), FDhistogram_old.GetXaxis().GetXbins().GetArray(), len(dptbins) - 1, array.array('d', dptbins))
    ApplyEfficiency(FDhistogram_fineBins, bEfficiency_dpt, False)
    FDhistogram_fineBins.SetName("GeneratorLevel_JetPtDPtSpectrum_FineBins_bEfficiencyMultiply")
    jetptdpt[FDhistogram_fineBins.GetName()] = FDhistogram_fineBins

    print("Applying the correction for the c->D0 reconstruction efficiency (fine bins)")
    FDhistogram_fineBins_efficiency = FDhistogram_fineBins.Clone("{0}_efficiency".format(FDhistogram_fineBins.GetName()))
    ApplyEfficiency(FDhistogram_fineBins_efficiency, cEfficiency_dpt, True)
    FDhistogram_fineBins_efficiency.SetName("GeneratorLevel_JetPtDPtSpectrum_FineBins_bEfficiencyMultiply_cEfficiencyDivide")
    jetptdpt[FDhistogram_fineBins_efficiency.GetName()] = FDhistogram_fineBins_efficiency

    for ptd in range(2, 5):
        spectrumName = "JetPtSpectrum_DPt_{0}".format(ptd * 10)

        ptdList = OrderedDict()
        result[spectrumName] = ptdList

        print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/o efficiency)".format(ptd))
        FDhistogram_jetpt = FDhistogram.ProjectionX("{0}_jetpt_DPt_{1}".format(FDhistogram.GetName(), ptd * 10), FDhistogram.GetYaxis().FindBin(ptd), FDhistogram.GetNbinsY() + 1)
        FDhistogram_jetpt.SetName("GeneratorLevel_JetPtSpectrum_bEfficiencyMultiply")
        FDhistogram_jetpt.GetYaxis().SetTitle(FDhistogram.GetZaxis().GetTitle())
        ptdList[FDhistogram_jetpt.GetName()] = FDhistogram_jetpt

        print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/o efficiency, fine bins)".format(ptd))
        FDhistogram_fineBins_jetpt = FDhistogram_fineBins.ProjectionX("{0}_jetpt_DPt_{1}".format(FDhistogram_fineBins.GetName(), ptd * 10), FDhistogram_fineBins.GetYaxis().FindBin(ptd), FDhistogram_fineBins.GetNbinsY() + 1)
        FDhistogram_fineBins_jetpt.SetName("GeneratorLevel_JetPtSpectrum_FineBins_bEfficiencyMultiply")
        FDhistogram_fineBins_jetpt.GetYaxis().SetTitle(FDhistogram_fineBins.GetZaxis().GetTitle())
        ptdList[FDhistogram_fineBins_jetpt.GetName()] = FDhistogram_fineBins_jetpt

        print("Applying the jet pt b response matrix")
        FDhistogram_jetpt_detector = ApplyResponse(FDhistogram_jetpt, bResponseFile, FDhistogram_fineBins_jetpt, spectrumName)
        FDhistogram_jetpt_detector.SetName("DetectorLevel_JetPtSpectrum_bEfficiencyMultiply")
        ptdList[FDhistogram_jetpt_detector.GetName()] = FDhistogram_jetpt_detector

        print("Unfolding using the c response matrix")
        FDhistogram_jetpt_unfolded = GetUnfoldedSpectrum("{0}_{1}_unfolded_c".format(ts, FDhistogram_jetpt_detector.GetName()), FDhistogram_jetpt_detector, cResponseFile, FDhistogram_fineBins_jetpt, spectrumName, unfolding_debug)
        FDhistogram_jetpt_unfolded.SetName("Unfolded_c_JetPtSpectrum_bEfficiencyMultiply")
        ptdList[FDhistogram_jetpt_unfolded.GetName()] = FDhistogram_jetpt_unfolded

        print("Unfolding using the b response matrix")
        FDhistogram_jetpt_unfolded = GetUnfoldedSpectrum("{0}_{1}_unfolded_b".format(ts, FDhistogram_jetpt_detector.GetName()), FDhistogram_jetpt_detector, bResponseFile, FDhistogram_fineBins_jetpt, spectrumName, unfolding_debug)
        FDhistogram_jetpt_unfolded.SetName("Unfolded_b_JetPtSpectrum_bEfficiencyMultiply")
        ptdList[FDhistogram_jetpt_unfolded.GetName()] = FDhistogram_jetpt_unfolded

        print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/ efficiency)".format(ptd))
        FDhistogram_efficiency_jetpt = FDhistogram_efficiency.ProjectionX("{0}_jetpt_DPt_{1}".format(FDhistogram_efficiency.GetName(), ptd * 10), FDhistogram_efficiency.GetYaxis().FindBin(ptd), FDhistogram_efficiency.GetNbinsY() + 1)
        FDhistogram_efficiency_jetpt.SetName("GeneratorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")
        FDhistogram_efficiency_jetpt.GetYaxis().SetTitle(FDhistogram_efficiency.GetZaxis().GetTitle())
        ptdList[FDhistogram_efficiency_jetpt.GetName()] = FDhistogram_efficiency_jetpt

        print("Projecting FD histogram into the jet pt axis for D pt > {0} GeV/c (w/ efficiency, fine bins)".format(ptd))
        FDhistogram_fineBins_efficiency_jetpt = FDhistogram_fineBins_efficiency.ProjectionX("{0}_jetpt_DPt_{1}".format(FDhistogram_fineBins_efficiency.GetName(), ptd * 10), FDhistogram_fineBins_efficiency.GetYaxis().FindBin(ptd), FDhistogram_fineBins_efficiency.GetNbinsY() + 1)
        FDhistogram_fineBins_efficiency_jetpt.SetName("GeneratorLevel_JetPtSpectrum_FineBins_bEfficiencyMultiply_cEfficiencyDivide")
        FDhistogram_fineBins_efficiency_jetpt.GetYaxis().SetTitle(FDhistogram_fineBins_efficiency.GetZaxis().GetTitle())
        ptdList[FDhistogram_fineBins_efficiency_jetpt.GetName()] = FDhistogram_fineBins_efficiency_jetpt

        print("Applying the jet pt b response matrix")
        FDhistogram_efficiency_jetpt_detector = ApplyResponse(FDhistogram_efficiency_jetpt, bResponseFile_efficiency, FDhistogram_fineBins_efficiency_jetpt, spectrumName)
        FDhistogram_efficiency_jetpt_detector.SetName("DetectorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")
        ptdList[FDhistogram_efficiency_jetpt_detector.GetName()] = FDhistogram_efficiency_jetpt_detector

        print("Unfolding using the c response matrix")
        FDhistogram_efficiency_jetpt_unfolded = GetUnfoldedSpectrum("{0}_{1}_unfolded_c".format(ts, FDhistogram_efficiency_jetpt_detector.GetName()), FDhistogram_efficiency_jetpt_detector, cResponseFile_efficiency, FDhistogram_fineBins_efficiency_jetpt, spectrumName, unfolding_debug)
        FDhistogram_efficiency_jetpt_unfolded.SetName("Unfolded_c_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")
        ptdList[FDhistogram_efficiency_jetpt_unfolded.GetName()] = FDhistogram_efficiency_jetpt_unfolded

        print("Unfolding using the b response matrix")
        FDhistogram_efficiency_jetpt_unfolded_b = GetUnfoldedSpectrum("{0}_{1}_unfolded_b".format(ts, FDhistogram_efficiency_jetpt_detector.GetName()), FDhistogram_efficiency_jetpt_detector, bResponseFile_efficiency, FDhistogram_fineBins_efficiency_jetpt, spectrumName, unfolding_debug)
        FDhistogram_efficiency_jetpt_unfolded_b.SetName("Unfolded_b_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide")
        ptdList[FDhistogram_efficiency_jetpt_unfolded_b.GetName()] = FDhistogram_efficiency_jetpt_unfolded_b

    return result

def GenerateUnfoldingEngine(name, responseFile, prior, spectrumName):
    analysis = dict()
    analysis["name"] = name
    analysis["d_meson"] = "_d_meson_"
    analysis["spectrum"] = "_spectrum_"
    analysis["jet_type"] = "Charged"
    analysis["jet_radius"] = "R040"
    analysis["active"] = True
    analysis["d_meson_response"] = "D0"
    analysis["spectrum_response"] = spectrumName
    analysis["priors"] = ["GeneratedSpectrum"]  # , "ResponseTruth", "PowerLaw_3", "PowerLaw_7"]
    # analysis["priors"] = ["ResponseTruth"] #, "ResponseTruth", "PowerLaw_3", "PowerLaw_7"]
    analysis["default_prior"] = "GeneratedSpectrum"
    analysis["default_method"] = "Svd"
    svd = dict()
    default_reg_svd = {"GeneratedSpectrum" : 6}  # , "PowerLaw_3" : 5, "PowerLaw_7" : 5}
    svd["default_reg"] = default_reg_svd
    bayes = dict()
    default_reg_bayes = {"GeneratedSpectrum" : 6}  # , "PowerLaw_3" : 5, "PowerLaw_7" : 5}
    bayes["default_reg"] = default_reg_bayes
    bayes["iter_min"] = 1
    bayes["iter_max"] = 6
    bayes["iter_step"] = 1
    binBybin = dict()
    binBybin["default_reg"] = None
    analysis["methods"] = {"Svd" : svd, "Bayes" : bayes, "BinByBin" : binBybin}
    unf = DMesonJetUnfolding.DMesonJetUnfoldingEngine(analysis)
    unf.fDoErrorCrossChecks = False
    unf.SetCustomPrior("GeneratedSpectrum", prior)
    unf.fUseOverflow = True
    unf.LoadResponse(responseFile, False)
    return unf

def GetUnfoldedSpectrum(name, histogram, responseFile, prior, spectrumName, unfolding_debug):
    unf = GenerateUnfoldingEngine(name, responseFile, prior, spectrumName)
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
    unf = GenerateUnfoldingEngine(truth.GetName(), responseFile, prior, spectrumName)
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
        print("Could not get list {0}".format(rlistName))
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
    file = ROOT.TFile(file_name)
    if not file or file.IsZombie():
        print("Could not open input file {0}".format(file_name))
        exit(1)
    else:
        print("File {0} open".format(file_name))
    dlist = file.Get("D0_MCTruth")
    jlist = dlist.FindObject("Charged_R040")
    slist = jlist.FindObject("D0_MCTruth_Charged_R040_JetPtDPtSpectrum")
    if not slist:
        print("Could not find FD histogram (jet pt)!")
        jlist.Print()
        exit(1)
    jetpt_hist = slist.FindObject("D0_MCTruth_Charged_R040_JetPtDPtSpectrum")
    if not jetpt_hist:
        print("Could not find FD histogram (jet pt)!")
        slist.Print()
        exit(1)
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

    jetpt_hist.GetZaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} #times #Delta#it{p}_{T} (mb)")
    dpt_hist.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} #times #Delta#it{p}_{T} (mb)")
    return jetpt_hist, dpt_hist

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
