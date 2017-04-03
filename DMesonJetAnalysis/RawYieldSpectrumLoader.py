#!/usr/bin/env python
# python script to do extract B feed down correction factors

import ROOT
import numpy
import DMesonJetFDCorrection
import DMesonJetUtils

class RawYieldSpectrumLoader:
    def __init__(self, input_path=None, train=None, ana_name=None):
        self.fInputPath = input_path
        self.fTrain = train
        self.fAnalysisName = ana_name
        self.fDPtBins = [3, 4, 5, 6, 7, 8, 10, 12, 16, 30]
        self.fJetPtBins = [5, 6, 8, 10, 14, 20, 30]
        self.fUseReflections = True
        self.fReflectionFit = "DoubleGaus"
        self.fReflectionRoS = 0
        self.fEvents = None
        self.fDMeson = None
        self.fJetType = None
        self.fJetRadius = None
        self.fSpectrumName = None
        self.fKinematicCuts = None
        self.fRawYieldMethod = None
        self.fDataSpectrumList = None
        self.fDataFile = None
        self.fMultiTrialSubDir = None
        self.fDataJetList = None
        self.fDataMesonList = None
        self.LoadNumberOfEvents = self.LoadNumberOfEventsNew

    def LoadNumberOfEventsOld(self):
        if not self.fDataSpectrumList: self.LoadDataListFromDMesonJetAnalysis()
        inputSpectrumName = "_".join([s for s in [self.fDMeson, self.fJetType, self.fJetRadius, self.fSpectrumName, self.fKinematicCuts, self.fRawYieldMethod] if s])
        inputSpectrum = self.fDataSpectrumList.FindObject(inputSpectrumName)
        if not inputSpectrum:
            print("Could not find histogram {0} in list {1}". format(inputSpectrumName, self.fDataSpectrumList.GetName()))
            self.fDataSpectrumList.Print()
            exit(1)
        normInputSpectrumName = "{0}_Normalized".format(inputSpectrumName)
        normInputSpectrum = self.fDataSpectrumList.FindObject(normInputSpectrumName)
        if not normInputSpectrum:
            print("Could not find histogram {0} in list {1}". format(normInputSpectrumName, self.fDataSpectrumList.GetName()))
            self.fDataSpectrumList.Print()
            exit(1)
        temp = inputSpectrum.Clone("temp")
        temp.Scale(1, "width")
        temp.Divide(normInputSpectrum)
        self.fEvents = temp.GetBinContent(1)
        return self.fEvents

    def LoadNumberOfEventsNew(self):
        periods = ["LHC10b", "LHC10c", "LHC10d", "LHC10e"]
        nRecoVertVz_LT10 = 0
        nRecoVertVz_GT10 = 0
        nNoVert = 0
        nVzSPD = 0
        for p in periods:
            fname = "{input_path}/{train}/{period}/merge/AnalysisResults.root".format(input_path=self.fInputPath, train=self.fTrain, period=p)
            file = ROOT.TFile(fname)
            if not file or file.IsZombie():
                print("Could not open file {0}".format(fname))
                exit(1)
            EMCALrejEvents = DMesonJetUtils.GetObject(file, "AliAnalysisTaskDmesonJets_AnyINT_histos/fHistEventRejection")
            RDHFrejEvents = DMesonJetUtils.GetObject(file, "AliAnalysisTaskDmesonJets_AnyINT_histos/histosAliAnalysisTaskDmesonJets_AnyINT/{dmeson}/fHistEventRejectionReasons".format(dmeson=self.fDMeson))
            accEvents = DMesonJetUtils.GetObject(file, "AliAnalysisTaskDmesonJets_AnyINT_histos/histosAliAnalysisTaskDmesonJets_AnyINT/{dmeson}/fHistNEvents".format(dmeson=self.fDMeson))
            nRecoVertVz_LT10 += accEvents.GetBinContent(accEvents.GetXaxis().FindBin("Accepted"))
            nRecoVertVz_GT10 += EMCALrejEvents.GetBinContent(EMCALrejEvents.GetXaxis().FindBin("Vz"))
            nNoVert += EMCALrejEvents.GetBinContent(EMCALrejEvents.GetXaxis().FindBin("vertex contr."))
            nNoVert += RDHFrejEvents.GetBinContent(RDHFrejEvents.GetXaxis().FindBin("NoVertex"))
            nVzSPD += EMCALrejEvents.GetBinContent(EMCALrejEvents.GetXaxis().FindBin("VzSPD"))
            file.Close()
        print("Total number of accepted events (Vz < 10 cm): {0:e}".format(nRecoVertVz_LT10))
        print("Total number of rejected events with reconstructed vertex Vz > 10 cm: {0:e}".format(nRecoVertVz_GT10))
        print("Total number of rejected events with reconstructed vertex |Vz_V0 - Vz_SPD| > 0.5 cm: {0:e}".format(nVzSPD))
        print("Total number of rejected events without reconstructed vertex: {0:e}".format(nNoVert))
        nRecoVert = nRecoVertVz_LT10 + nRecoVertVz_GT10
        normEvents = nRecoVertVz_LT10 + nNoVert * (1 - nRecoVertVz_GT10 / nRecoVert)
        print("Number of events with Vz < 10 cm (regardless of whether the vertex was reconstructed): {0:e}".format(normEvents))
        print("Corrective factor effective / actual accepted events: {0:e}".format(normEvents / nRecoVertVz_LT10))
        print("Fraction of events cut because of SPD cut: {0:e}".format(nVzSPD / normEvents))
        self.fEvents = normEvents
        return self.fEvents

    def LoadDataFileFromDMesonJetAnalysis(self):
        self.fDataFile = ROOT.TFile("{0}/{1}/{2}.root".format(self.fInputPath, self.fTrain, self.fAnalysisName))
        return self.fDataFile

    def LoadDataListFromDMesonJetAnalysis(self,):
        if not self.fDataFile: self.LoadDataFileFromDMesonJetAnalysis()
        if self.fDMeson:
            self.fDataMesonList = self.fDataFile.Get(self.fDMeson)
            if not self.fDataMesonList:
                print("Could not find list {0} in file {1}". format(self.fDMeson, self.fDataFile.GetName()))
                exit(1)
            if self.fJetType and self.fJetRadius:
                dataJetListName = "_".join([self.fJetType, self.fJetRadius])
                self.fDataJetList = self.fDataMesonList.FindObject(dataJetListName)
                if not self.fDataJetList:
                    print("Could not find list {0}/{1} in file {2}". format(self.fDMeson, dataJetListName, self.fDataFile.GetName()))
                    self.fDataMesonList.Print()
                    exit(1)
                if self.fSpectrumName:
                    dataListName = "_".join([s for s in [self.fDMeson, dataJetListName, self.fSpectrumName, self.fKinematicCuts, self.fRawYieldMethod] if s])
                    self.fDataSpectrumList = self.fDataJetList.FindObject(dataListName)
                    if not self.fDataSpectrumList:
                        print("Could not find list {0}/{1}/{2} in file {3}". format(self.fDMeson, dataJetListName, dataListName, self.fDataFile.GetName()))
                        self.fDataJetList.Print()
                        exit(1)
            if self.fSpectrumName and not (self.fJetType and self.fJetRadius):
                dataListName = "_".join([s for s in [self.fDMeson, self.fSpectrumName, self.fKinematicCuts, self.fRawYieldMethod] if s])
                self.fDataSpectrumList = self.fDataMesonList.FindObject(dataListName)
                if not self.fDataSpectrumList:
                    print("Could not find list {0}/{1} in file {2}". format(self.fDMeson, dataListName, self.fDataFile.GetName()))
                    self.fDataMesonList.Print()
                    exit(1)
        return self.fDataSpectrumList

    def GetDefaultSpectrumFromDMesonJetAnalysis(self, fd=False, fd_error_band=0, ry_error_band=0):
        if not self.fSpectrumName:
            print("No spectrum name provided!")
            exit(1)
        if self.fUseReflections:
            print("****Attention Attention Attention****")
            print("You asked for reflections, but reflections are not available in DMesonJetAnalysis!")
            exit(1)
        if not self.fDataSpectrumList: self.LoadDataListFromDMesonJetAnalysis()
        inputSpectrumName = "_".join([s for s in [self.fDMeson, self.fJetType, self.fJetRadius, self.fSpectrumName, self.fKinematicCuts, self.fRawYieldMethod] if s])
        h = self.fDataSpectrumList.FindObject(inputSpectrumName)
        if not h:
            print("Could not find histogram {0} in list {1}". format(inputSpectrumName, self.fDataSpectrumList.GetName()))
            self.fDataSpectrumList.Print()
            exit(1)
        if fd: self.ApplyFDCorrection(h, fd_error_band)
        if ry_error_band <> 0: self.ApplyRawYieldSystFromMultiTrial(h, ry_error_band)
        return h

    def GetDefaultSpectrumInvMassFitFromMultiTrial(self):
        if not self.fSpectrumName:
            print("No spectrum name provided!")
            exit(1)
        if self.fUseReflections:
            self.fMultiTrialSubDir = "RawYieldUnc_refl_{0}".format(self.fReflectionFit)
            if not self.fReflectionRoS == 0: self.fMultiTrialSubDir += "_{0}".format(self.fReflectionRoS)
        else:
            self.fMultiTrialSubDir = "RawYieldUnc"
        h = ROOT.TH1D("TrialExpoFreeSigmaInvMassFit", "Trial Expo Free Sigma Inv. Mass Fit", len(self.fJetPtBins) - 1, numpy.array(self.fJetPtBins, dtype=numpy.float64))
        for ibin, (ptmin, ptmax) in enumerate(zip(self.fJetPtBins[:-1], self.fJetPtBins[1:])):
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
        if not self.fSpectrumName:
            print("No spectrum name provided!")
            exit(1)
        if self.fUseReflections:
            self.fMultiTrialSubDir = "RawYieldUnc_refl_{0}".format(self.fReflectionFit)
            if not self.fReflectionRoS == 0: self.fMultiTrialSubDir += "_{0}".format(self.fReflectionRoS)
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

    def GetDefaultSpectrumFromMultiTrial(self, fd=False, fd_error_band=0, ry_error_band=0):
        if self.fRawYieldMethod == "SideBand": h = self.GetDefaultSpectrumSideBandFromMultiTrial()
        elif self.fRawYieldMethod == "InvMassFit": h = self.GetDefaultSpectrumInvMassFitFromMultiTrial()
        else:
            print("Method {0} unknown!".format(self.fRawYieldMethod))
            exit(1)
        if fd: self.ApplyFDCorrection(h, fd_error_band)
        if ry_error_band <> 0: self.ApplyRawYieldSystFromMultiTrial(h, ry_error_band)
        return h

    def GetAverageSpectrumFromMultiTrial(self, fd=False, fd_error_band=0):
        if not self.fSpectrumName:
            print("No spectrum name provided!")
            exit(1)
        if self.fUseReflections:
            self.fMultiTrialSubDir = "RawYieldUnc_refl_{0}".format(self.fReflectionFit)
            if not self.fReflectionRoS == 0: self.fMultiTrialSubDir += "_{0}".format(self.fReflectionRoS)
        else:
            self.fMultiTrialSubDir = "RawYieldUnc"
        fname = "{0}/{1}/{2}/{3}/FinalRawYieldCentralPlusSystUncertainty_Dzero_{4}.root".format(self.fInputPath, self.fTrain, self.fAnalysisName, self.fMultiTrialSubDir, self.fRawYieldMethod)
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
        fdHist = self.GetFDCorrection(-error_band)
        h.Add(fdHist, -1)

    def GetFDCorrection(self, error_band=0):
        fdConfig = dict()
        fdConfig["file_name"] = "BFeedDown.root"
        fdConfig["central_points"] = "default"
        if "JetPt" in self.fSpectrumName:
            if "efficiency" in self.fAnalysisName:
                fdConfig["spectrum"] = "DetectorLevel_JetPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide"
                print("Using FD correction with efficiency")
            else:
                fdConfig["spectrum"] = "DetectorLevel_JetPtSpectrum_bEfficiencyMultiply"
                print("Using FD correction without efficiency")
        elif "DPt" in self.fSpectrumName:
            if "efficiency" in self.fAnalysisName:
                fdConfig["spectrum"] = "GeneratorLevel_DPtSpectrum_bEfficiencyMultiply_cEfficiencyDivide"
                print("Using FD correction with efficiency")
            else:
                fdConfig["spectrum"] = "DetectorLevel_DPtSpectrum_bEfficiencyMultiply"
                print("Using FD correction without efficiency")
        else:
            print("ApplyFDCorrection: Not jet pt nor d pt?")
            exit(1)

        spectrumName = self.fSpectrumName
        if self.fKinematicCuts: spectrumName += "_{0}".format(self.fKinematicCuts)
        fdCorrection = DMesonJetFDCorrection.DMesonJetFDCorrection(fdConfig, spectrumName, self.fInputPath, "D0", "Charged", "R040")
        fdHist = fdCorrection.fFDHistogram.Clone("FD")
        fdSyst = dict()
        fdSyst[1] = fdCorrection.fFDUpSystUncHistogram.Clone("FDUpSystUnc")
        fdSyst[-1] = fdCorrection.fFDLowSystUncHistogram.Clone("FDLowSystUnc")
        fdSystGraph = fdCorrection.fFDSystUncGraph.Clone("FDSystUnc")

        crossSection = 62.2  # mb CINT1
        branchingRatio = 0.0393  # D0->Kpi
        if self.fEvents is None: self.LoadNumberOfEvents()
        norm = self.fEvents / crossSection * branchingRatio
        fdHist.Scale(norm)
        fdSyst[1].Scale(norm)
        fdSyst[-1].Scale(norm)
        if error_band == "graph":
            for ibin in xrange(0, fdSystGraph.GetN()):
                fdSystGraph.SetPoint(ibin, fdSystGraph.GetX()[ibin], fdSystGraph.GetY()[ibin] * norm * fdHist.GetXaxis().GetBinWidth(ibin + 1))
                fdSystGraph.SetPointEYlow(ibin, fdSystGraph.GetErrorYlow(ibin) * norm * fdHist.GetXaxis().GetBinWidth(ibin + 1))
                fdSystGraph.SetPointEYhigh(ibin, fdSystGraph.GetErrorYhigh(ibin) * norm * fdHist.GetXaxis().GetBinWidth(ibin + 1))
            return fdSystGraph
        elif error_band <> 0:
            fdHist.Add(fdSyst[error_band], error_band)
        return fdHist

    def ApplyRawYieldSystFromMultiTrial(self, h, error_band):
        if error_band == 0:
            print("No raw yield systematic to apply!")
            return
        print("Applying raw yield systematic with error band point: {0} (0 = central point, +/- = upper/lower error band)".format(error_band))
        if self.fUseReflections:
            self.fMultiTrialSubDir = "RawYieldUnc_refl_{0}".format(self.fReflectionFit)
            if not self.fReflectionRoS == 0: self.fMultiTrialSubDir += "_{0}".format(self.fReflectionRoS)
        else:
            self.fMultiTrialSubDir = "RawYieldUnc"
        fname = "{0}/{1}/{2}/{3}/FinalRawYieldUncertainty_Dzero_{4}.root".format(self.fInputPath, self.fTrain, self.fAnalysisName, self.fMultiTrialSubDir, self.fRawYieldMethod)
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
