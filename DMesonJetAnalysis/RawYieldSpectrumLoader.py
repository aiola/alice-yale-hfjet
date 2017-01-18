#!/usr/bin/env python
# python script to do extract B feed down correction factors

import ROOT
import numpy
import DMesonJetFDCorrection

ptDbins = [3, 4, 5, 6, 7, 8, 10, 12, 16, 30]
ptJetbins = [5, 6, 8, 10, 14, 20, 30]  # used for eff.scale approach, but also in sideband approach to define the bins of the output jet spectrum

class RawYieldSpectrumLoader:
    def __init__(self, input_path=None, train=None, ana_name=None, refl=False, events=None):
        self.fInputPath = input_path
        self.fTrain = train
        self.fAnalysisName = ana_name
        self.fUseReflections = refl
        self.fEvents = events
        self.fDMeson = None
        self.fJetType = None
        self.fJetRadius = None
        self.fSpectrumName = None
        self.fDataList = None
        self.fDataFile = None
        self.fMultiTrialSubDir = None

    def LoadDataFileFromDMesonJetAnalysis(self):
        self.fDataFile = ROOT.TFile("{0}/{1}/{2}.root".format(self.fInputPath, self.fTrain, self.fAnalysisName))
        return self.fDataFile

    def LoadDataListFromDMesonJetAnalysis(self):
        if not self.fDataFile: self.LoadDataFileFromDMesonJetAnalysis()
        dataMesonList = self.fDataFile.Get(self.fDMeson)
        if not dataMesonList:
            print("Could not find list {0} in file {1}". format(self.fDMeson, self.fDataFile.GetName()))
            exit(1)
        dataJetListName = "_".join([self.fJetType, self.fJetRadius])
        dataJetList = dataMesonList.FindObject(dataJetListName)
        if not dataJetList:
            print("Could not find list {0}/{1} in file {2}". format(self.fDMeson, dataJetListName, self.fDataFile.GetName()))
            dataMesonList.Print()
            exit(1)
        dataListName = "{0}_{1}_{2}".format(self.fDMeson, dataJetListName, self.fSpectrumName)
        self.fDataList = dataJetList.FindObject(dataListName)
        if not self.fDataList:
            print("Could not find list {0}/{1}/{2} in file {3}". format(self.fDMeson, dataJetListName, dataListName, self.fDataFile.GetName()))
            dataJetList.Print()
            exit(1)
        return self.fDataList

    def GetDefaultSpectrumFromDMesonJetAnalysis(self, method, fd=False, fd_error_band=0, ry_error_band=0):
        if not (method == "SideBand" or method == "InvMassFit"):
            print("Method {0} unknown!".format(method))
            exit(1)
        if self.fUseReflections:
            print("****Attention Attention Attention****")
            print("You asked for reflections, but reflections are not available outside of the multi-trial code!")
            print("The reflection option will be ignored!")
        if not self.fDataList: self.LoadDataListFromDMesonJetAnalysis()
        inputSpectrumName = "_".join([self.fDMeson, self.fJetType, self.fJetRadius, self.fSpectrumName])
        h = self.fDataList.FindObject(inputSpectrumName)
        if not h:
            print("Could not find histogram {0} in list {1}". format(inputSpectrumName, self.fDataList.GetName()))
            exit(1)
        if fd: self.ApplyFDCorrection(h, fd_error_band)
        if ry_error_band <> 0: self.ApplyRawYieldSystFromMultiTrial(h, method, ry_error_band)
        return h

    def GetDefaultSpectrumInvMassFitFromMultiTrial(self):
        if self.fUseReflections:
            self.fMultiTrialSubDir = "RawYieldUnc_refl"
        else:
            self.fMultiTrialSubDir = "RawYieldUnc"
        h = ROOT.TH1D("TrialExpoFreeSigmaInvMassFit", "Trial Expo Free Sigma Inv. Mass Fit", len(ptJetbins) - 1, numpy.array(ptJetbins, dtype=numpy.float64))
        for ibin, (ptmin, ptmax) in enumerate(zip(ptJetbins[:-1], ptJetbins[1:])):
            fname = "{0}/{1}/{2}/{3}/RawYieldVariations_Dzero_InvMassFit_{4:.1f}to{5:.1f}.root".format(self.fInputPath, self.fTrain, self.fAnalysisName, self.fMultiTrialSubDir, ptmin, ptmax)
            file = ROOT.TFile(fname)
            if not file or file.IsZombie():
                print("Could not open file {0}".format(fname))
                file.ls()
                exit(1)
            else:
                print("File {0} open successfully".format(fname))
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

    def GetDefaultSpectrumSideBandFromMultiTrial(self):
        if self.fUseReflections:
            self.fMultiTrialSubDir = "RawYieldUnc_refl"
        else:
            self.fMultiTrialSubDir = "RawYieldUnc"
        fname = "{0}/{1}/{2}/{3}/TrialExpoFreeS_Dzero_SideBand.root".format(self.fInputPath, self.fTrain, self.fAnalysisName, self.fMultiTrialSubDir)
        file = ROOT.TFile(fname)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(fname))
            file.ls()
            exit(1)
        else:
            print("File {0} open successfully".format(fname))
        h = file.Get("fJetSpectrSBDef")
        if not h:
            print("Could not find histogram {0} in file {1}".format("fJetSpectrSBDef", fname))
            file.ls()
            exit(1)
        h_copy = h.Clone("TrialExpoFreeSigmaSideBand")
        h_copy.SetTitle("Trial Expo Free Sigma Inv. Mass Fit")
        return h

    def GetDefaultSpectrumFromMultiTrial(self, method, fd=False, fd_error_band=0, ry_error_band=0):
        if method == "SideBand": h = self.GetDefaultSpectrumSideBandFromMultiTrial()
        elif method == "InvMassFit": h = self.GetDefaultSpectrumInvMassFitFromMultiTrial()
        else:
            print("Method {0} unknown!".format(method))
            exit(1)
        if fd: self.ApplyFDCorrection(h, fd_error_band)
        if ry_error_band <> 0: self.ApplyRawYieldSystFromMultiTrial(h, method, ry_error_band)
        return h

    def GetAverageSpectrumFromMultiTrial(self, method, fd=False, fd_error_band=0):
        fname = "{0}/{1}/{2}/RawYieldUnc/FinalRawYieldCentralPlusSystUncertainty_Dzero_{3}.root".format(self.fInputPath, self.fTrain, self.fAnalysisName, method)
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
        if fd: self.ApplyFDCorrection(h_copy, fd_error_band)
        return h_copy

    def ApplyFDCorrection(self, h, error_band=0):
        print("Applying FD correction with error band point: {0} (0 = central point, +/- = upper/lower error band)".format(error_band))
        fdConfig = dict()
        fdConfig["file_name"] = "BFeedDown.root"
        fdConfig["central_points"] = "default"
        if "efficiency" in self.fAnalysisName:
            fdConfig["spectrum"] = "DetectorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide"
            print("Using FD correction with efficiency")
        else:
            fdConfig["spectrum"] = "DetectorLevel_JetPtSpectrum_bEfficiencyMultiply"
            print("Using FD correction without efficiency")

        fdCorrection = DMesonJetFDCorrection.DMesonJetFDCorrection(fdConfig, "JetPtSpectrum_DPt_30", self.fInputPath, "D0", "Charged", "R040")
        fdHist = fdCorrection.fFDHistogram.Clone("{0}_FD".format(h.GetName()))
        fdSyst = fdCorrection.fFDSystUncHistogram.Clone("{0}_FDSystUnc".format(h.GetName()))

        crossSection = 62.3  # mb CINT1
        branchingRatio = 0.0393  # D0->Kpi
        fdHist.Scale(self.fEvents / crossSection * branchingRatio)
        fdSyst.Scale(self.fEvents / crossSection * branchingRatio)
        h.Add(fdHist, -1)
        if error_band <> 0: h.Add(fdSyst, error_band)

    def ApplyRawYieldSystFromMultiTrial(self, h, method, error_band):
        if error_band == 0:
            print("No raw yield systematic to apply!")
            return
        print("Applying raw yield systematic with error band point: {0} (0 = central point, +/- = upper/lower error band)".format(error_band))

        fname = "{0}/{1}/{2}/RawYieldUnc/FinalRawYieldUncertainty_Dzero_{3}.root".format(self.fInputPath, self.fTrain, self.fAnalysisName, method)
        file = ROOT.TFile(fname)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(fname))
            file.ls()
            exit(1)
        rySyst = file.Get("JetRawYieldUncert")
        if not rySyst:
            print("Could not find histogram {0} in file {1}".format("JetRawYieldUncert", fname))
            file.ls()
            exit(1)

        h.Add(rySyst, error_band)
