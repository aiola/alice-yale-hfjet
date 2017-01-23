#!/usr/bin/env python
# python script to do extract B feed down correction factors

import yaml
import IPython
import ROOT
import DMesonJetUtils
import RawYieldSpectrumLoader
import subprocess
import DMesonJetCompare
import numpy

globalList = []

input_path = "/Volumes/DATA/ALICE/JetResults"
ana = "Jets_EMC_pp_823_824_825_826/LHC10_Train823_efficiency"

def SideBandFinalRawYieldUnc():
    fname = "{0}/{1}/RawYieldUnc_refl_DoubleGaus/DistributionOfFinalYields_SBApproach_Dzero_AfterDbinSum.root".format(input_path, ana)
    file = ROOT.TFile(fname)
    canvas = file.Get("cDistr_Dzero_SideBand")
    histos = []
    for obj in canvas.GetListOfPrimitives():
        if isinstance(obj, ROOT.TH1):
            n = len(histos)
            h_copy = obj.Clone("var{0}".format(n))
            h_copy.SetTitle("Variation {0}".format(n))
            histos.append(h_copy)
    file.Close()
    fname = "{0}/{1}/RawYieldUnc_refl_DoubleGaus/TrialExpoFreeS_Dzero_SideBand.root".format(input_path, ana)
    file = ROOT.TFile(fname)
    baseline = file.Get("fJetSpectrSBDef")
    baseline_copy = baseline.Clone("baseline")
    baseline_copy.SetTitle("TrialExpoFreeS")
    file.Close()
    globalList.append(baseline_copy)
    globalList.extend(histos)
    cname = "CompareRawYieldUncVariations_AfterDbinSum"
    comp = DMesonJetCompare.DMesonJetCompare(cname)
    comp.fOptRatio = "hist"
    r = comp.CompareSpectra(baseline_copy, histos)
    for obj in r:
        if not obj in globalList:
            globalList.append(obj)

def SideBandRawYieldUnc():
    for ibin in range(0, 9):
        fname = "{0}/{1}/RawYieldUnc_refl_DoubleGaus/DistributionOfFinalYields_SBApproach_Dzero_SideBand_{2}.root".format(input_path, ana, ibin)
        file = ROOT.TFile(fname)
        canvas = file.Get("cDistr_Dzero_SideBand_{0}".format(ibin))
        histos = []
        for obj in canvas.GetListOfPrimitives():
            if isinstance(obj, ROOT.TH1):
                n = len(histos)
                h_copy = obj.Clone("bin{0}_var{1}".format(ibin, n))
                h_copy.SetTitle("Variation {0}".format(n))
                histos.append(h_copy)
        file.Close()
        fname = "{0}/{1}/RawYieldUnc_refl_DoubleGaus/TrialExpoFreeS_Dzero_SideBand_{2}.root".format(input_path, ana, ibin)
        file = ROOT.TFile(fname)
        baseline = file.Get("hjetpt{0}".format(ibin))
        baseline_copy = baseline.Clone("bin{0}_baseline".format(ibin))
        baseline_copy.SetTitle("TrialExpoFreeS")
        file.Close()
        globalList.append(baseline_copy)
        globalList.extend(histos)
        cname = "CompareRawYieldUncVariations_{0}".format(ibin)
        comp = DMesonJetCompare.DMesonJetCompare(cname)
        comp.fOptRatio = "hist"
        r = comp.CompareSpectra(baseline_copy, histos)
        for obj in r:
            if not obj in globalList:
                globalList.append(obj)

def SideBandRawYieldReflUnc():
    reflVar = ["DoubleGaus_15", "DoubleGaus_5", "gaus", "pol3", "pol6"]
    for ibin in range(0, 9):
        histos = []
        for var in reflVar:
            fname = "{0}/{1}/RawYieldUnc_refl_{2}/TrialExpoFreeS_Dzero_SideBand_{3}.root".format(input_path, ana, var, ibin)
            file = ROOT.TFile(fname)
            h = file.Get("hjetpt{0}".format(ibin))
            h_copy = h.Clone("bin{0}_var{1}".format(ibin, var))
            h_copy.SetTitle("TrialExpoFreeS, refl={0}".format(var))
            file.Close()
            histos.append(h_copy)
        fname = "{0}/{1}/RawYieldUnc_refl_DoubleGaus/TrialExpoFreeS_Dzero_SideBand_{2}.root".format(input_path, ana, ibin)
        file = ROOT.TFile(fname)
        baseline = file.Get("hjetpt{0}".format(ibin))
        baseline_copy = baseline.Clone("bin{0}_baseline".format(ibin))
        baseline_copy.SetTitle("TrialExpoFreeS, refl=DoubleGaus")
        file.Close()
        globalList.append(baseline_copy)
        globalList.extend(histos)
        cname = "CompareRawYieldUncReflVariations_{0}".format(ibin)
        comp = DMesonJetCompare.DMesonJetCompare(cname)
        comp.fOptRatio = "hist"
        r = comp.CompareSpectra(baseline_copy, histos)
        for obj in r:
            if not obj in globalList:
                globalList.append(obj)

def main():
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    SideBandRawYieldUnc()
    SideBandRawYieldReflUnc()
    SideBandFinalRawYieldUnc()

    for obj in globalList:
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("{0}/{1}.pdf".format(input_path, obj.GetName()))

if __name__ == '__main__':
    main()

    IPython.embed()
