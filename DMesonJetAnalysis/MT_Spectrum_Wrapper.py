#!/usr/bin/env python
# python script to do extract B feed down correction factors

import ROOT

def GetDefaultSpectrumInvMassFitFromMultiTrial(config):
    h = ROOT.TH1D("TrialExpoFreeSigmaInvMassFit", "Trial Expo Free Sigma Inv. Mass Fit", len(ptJetbins) - 1, numpy.array(ptJetbins, dtype=numpy.float64))
    for ibin, (ptmin, ptmax) in enumerate(zip(ptJetbins[:-1], ptJetbins[1:])):
        fname = "{0}/{1}/{2}/RawYieldUnc/RawYieldVariations_Dzero_InvMassFit_{3:.1f}to{4:.1f}.root".format(config["input_path"], config["train"], config["name"], ptmin, ptmax)
        file = ROOT.TFile(fname)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(fname))
            file.ls()
            exit(1)
        hry = file.Get("hRawYieldTrialExpoFreeS")
        if not hry:
            print("Could not find histogram {0} in file {1}".format("hRawYieldTrialExpoFreeS", fname))
            file.ls()
            exit(1)
        ry = hry.GetBinContent(1)  # this should be with the normal inv mass binning
        ery = hry.GetBinError(1)
        h.SetBinContent(ibin + 1, ry)
        h.SetBinError(ibin + 1, ery)
    return h

def GetDefaultSpectrumSideBandFromMultiTrial(config):
    fname = "{0}/{1}/{2}/RawYieldUnc/TrialExpoFreeS_Dzero_SideBand.root".format(config["input_path"], config["train"], config["name"])
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        file.ls()
        exit(1)
    h = file.Get("fJetSpectrSBDef")
    if not h:
        print("Could not find histogram {0} in file {1}".format("fJetSpectrSBDef", fname))
        file.ls()
        exit(1)
    h_copy = h.Clone("TrialExpoFreeSigmaSideBand")
    h_copy.SetTitle("Trial Expo Free Sigma Inv. Mass Fit")
    return h

def GetDefaultSpectrumFromMultiTrial(method, config, fd=False, fd_error_band=0):
    if method == "SideBand": h = GetDefaultSpectrumSideBandFromMultiTrial(config)
    elif method == "InvMassFit": h = GetDefaultSpectrumInvMassFitFromMultiTrial(config)
    else:
        print("Method {0} unknown!".format(method))
        exit(1)
    if fd: ApplyFDCorrection(h, config, fd_error_band)
    return h

def GetAverageSpectrum(method, config, fd=False, fd_error_band=0):
    fname = "{0}/{1}/{2}/RawYieldUnc/FinalRawYieldCentralPlusSystUncertainty_Dzero_{3}.root".format(config["input_path"], config["train"], config["name"], method)
    file = ROOT.TFile(fname)
    if not file or file.IsZombie():
        print("Could not open file {0}".format(fname))
        file.ls()
        exit(1)
    h = file.Get("JetRawYieldCentral")
    if not h:
        print("Could not find histogram {0} in file {1}".format("JetRawYieldCentral", fname))
        file.ls()
        exit(1)
    h_copy = h.Clone("{0}_copy".format(h.GetName()))
    if fd: ApplyFDCorrection(h_copy, config, fd_error_band)
    return h_copy

def GetNumberOfEvents(config):
    file = ROOT.TFile("{0}/{1}/{2}.root".format(config["input_path"], config["train"], config["name"]))
    dataMesonList = dataFile.Get("D0")
    if not dataMesonList:
        print("Could not find list {0} in file {1}". format("D0", dataFile.GetName()))
        exit(1)
    dataJetListName = "_".join(["Charged", "R040"])
    dataJetList = dataMesonList.FindObject(dataJetListName)
    if not dataJetList:
        print("Could not find list {0}/{1} in file {2}". format("D0", dataJetListName, dataFile.GetName()))
        dataMesonList.Print()
        exit(1)
    dataListName = "{0}_{1}_{2}".format("D0", dataJetListName, self.fSpectrumName)
    dataList = dataJetList.FindObject(dataListName)
    if not dataList:
        print("Could not find list {0}/{1}/{2} in file {3}". format("D0", dataJetListName, dataListName, dataFile.GetName()))
        dataJetList.Print()
        exit(1)
    inputSpectrumName = "{0}_{1}_{2}".format("D0", dataJetListName, "JetPtSpectrum_DPt_30_SideBand")
    inputSpectrum = dataList.FindObject(inputSpectrumName)
    if not inputSpectrum:
        print("Could not find histogram {0} in list {1} in file {2}". format(inputSpectrumName, dataListName, dataFile.GetName()))
        exit(1)
    normInputSpectrumName = "{0}_Normalized".format(inputSpectrumName)
    normInputSpectrum = dataList.FindObject(normInputSpectrumName)
    if not normInputSpectrum:
        print("Could not find histogram {0} in list {1}". format(normInputSpectrumName, dataList.GetName()))
        dataList.Print()
        exit(1)
    temp = self.fInputSpectrum.Clone("temp")
    temp.Scale(1, "width")
    temp.Divide(normInputSpectrum)
    return temp.GetBinContent(1)

def ApplyFDCorrection(h, config, error_band=0):
    fdConfig = dict()
    fdConfig["file_name"] = "BFeedDown.root"
    fdConfig["central_points"] = "default"
    if "efficiency" in  config["name"]:
        fdConfig["spectrum"] = "DetectorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide"
        print("Using FD correction with efficiency")
    else:
        fdConfig["spectrum"] = "DetectorLevel_JetPtSpectrum_bEfficiencyMultiply"
        print("Using FD correction without efficiency")

    fdCorrection = DMesonJetFDCorrection.DMesonJetFDCorrection(fdConfig, "JetPtSpectrum_DPt_30_SideBand", config["inputPath"], "D0", "Charged", "R040")
    fdHist = fdCorrection.fFDHistogram.Clone("{0}_FD".format(h.GetName()))
    fdSyst = fdCorrection.fFDSystUncHistogram.Clone("{0}_FDSystUnc".format(h.GetName()))

    events = GetNumberOfEvents(config)
    crossSection = 62.3  # mb CINT1
    branchingRatio = 0.0393  # D0->Kpi
    fdHist.Scale(events / crossSection * branchingRatio)
    fdSyst.Scale(events / crossSection * branchingRatio)
    h.Add(fdHist, -1)
    if error_band > 0: h.Add(fdSyst, error_band)
