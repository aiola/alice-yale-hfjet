#!/usr/bin/env python
#python program to perform a D meson jet unfolding

import ROOT
import math
import DMesonJetProjectors
from DMesonJetBase import *
import array
import copy
import DMesonJetUtils
from collections import OrderedDict

globalList = []

class DMesonJetUnfoldingEngine:
    def __init__(self, config):
        self.fDMeson = config["d_meson"]
        self.fDMesonTruth = config["d_meson_truth"]
        self.fDMesonResponse = config["d_meson_response"]
        self.fJetDefinition = config["jet"]
        self.fSpectrumName = config["spectrum"]
        self.fPriors = config["priors"]
        self.UnfoldingConfig = config["methods"]
        self.DefaultPrior = config["default_prior"]
        self.DefaultMethod = config["default_method"]
        self.fResponseMatrices = OrderedDict()
        self.fUnfoldedSpectra = OrderedDict()
        self.fRefoldedSpectra = OrderedDict()
        self.fPriorSpectra = OrderedDict()
        self.fCovarianceMatrices = OrderedDict()
        self.fSvdDvectors = OrderedDict()
        self.fPearsonMatrices = OrderedDict()

    def LoadData(self, dataFile, responseFile):
        dataListName = "{0}_{1}".format(self.fDMeson, self.fSpectrumName)
        dataList = dataFile.Get(dataListName)
        if not dataList:
            print("Could not find list {0} in file {1}". format(dataListName, dataFile.GetName()))
            return False
        self.fInputSpectrum = dataList.FindObject("{0}_{1}".format(self.fDMeson, self.fSpectrumName))

        if self.fDMesonTruth:
            dataTruthListName = "{0}_{1}".format(self.fDMesonTruth, self.fSpectrumName)
            dataTruthList = dataFile.Get(dataTruthListName)
            if not dataTruthList:
                print("Could not find list {0} in file {1}". format(dataTruthListName, dataFile.GetName()))
                return False
            self.fTruthSpectrum = dataTruthList.FindObject("{0}_{1}".format(self.fDMesonTruth, self.fSpectrumName))
        else:
            self.fTruthSpectrum = None

        responseListName = "{0}_{1}_{2}".format(self.fDMesonResponse, self.fJetDefinition, self.fSpectrumName)
        responseList = responseFile.Get(responseListName)
        if not responseList:
            print("Could not find list {0} in file {1}". format(responseListName, responseFile.GetName()))
            return False
        self.fDetectorResponse = responseList.FindObject("{0}_{1}_{2}_DetectorResponse".format(self.fDMesonResponse, self.fJetDefinition, self.fSpectrumName))
        self.fDetectorEfficiency = responseList.FindObject("{0}_{1}_{2}_DetectorEfficiency".format(self.fDMesonResponse, self.fJetDefinition, self.fSpectrumName))
        self.fDetectorTrainTruth = responseList.FindObject("{0}_{1}_{2}_Truth".format(self.fDMesonResponse, self.fJetDefinition, self.fSpectrumName))
        self.fDetectorRecontructedTruth = responseList.FindObject("{0}_{1}_{2}_ReconstructedTruth".format(self.fDMesonResponse, self.fJetDefinition, self.fSpectrumName))
        return True

    def Start(self):
        self.GenerateRooUnfoldResponse()
        for method, config in self.UnfoldingConfig.iteritems():
            if method == "Svd":
                self.UnfoldSvd()
            elif method == "Bayes":
                self.UnfoldBayes(config["iter_min"], config["iter_max"], config["iter_step"])
            else:
                print("Unfolding method {0} not known!".format(method))
        self.GeneratePearsonMatrices()
        self.Plot()
        
    def GeneratePearsonMatrices(self):
        for id, cov in self.fCovarianceMatrices.iteritems():
            pearson = cov.Clone(cov.GetName().replace("CovMat", "Pearson"))
            pearson.Reset()
            for x in range(1, pearson.GetNbinsX()+1):
                if not cov.GetBinContent(x,x) > 0:
                    continue
                sigmax = math.sqrt(cov.GetBinContent(x,x))
                for y in range(1, pearson.GetNbinsY()+1):
                    if not cov.GetBinContent(y,y) > 0:
                        continue
                    sigmay = math.sqrt(cov.GetBinContent(y,y))
                    pearson.SetBinContent(x,y,cov.GetBinContent(x,y) / sigmax / sigmay)
            self.fPearsonMatrices[id] = pearson

    def Plot(self):
        if len(self.UnfoldingConfig) > 1 or self.fTruthSpectrum:
            self.CompareMethods()
        self.UnfoldingSummary()
        self.PlotPearsonMatrices()
        self.PlotSvdDvectors()

    def PlotPearsonMatrices(self):
        for (method, reg, prior), cov in self.fPearsonMatrices.iteritems():
            cname = cov.GetName()
            c = ROOT.TCanvas(cname, cname) 
            c.cd()
            globalList.append(c)
            cov_copy = cov.DrawCopy("colz")
            cov_copy.GetZaxis().SetRangeUser(-1,1)
            globalList.append(cov_copy)

    def PlotSvdDvectors(self):
        for prior, d in self.fSvdDvectors.iteritems():
            cname = d.GetName()
            c = ROOT.TCanvas(cname, cname) 
            c.cd()
            c.SetLogy()
            globalList.append(c)
            d_copy = d.DrawCopy()
            globalList.append(d_copy)

    def UnfoldingSummary(self):
        for method, config in self.UnfoldingConfig.iteritems():
            self.UnfoldingSummaryForMethod(method, config["default_reg"], self.DefaultPrior)
            for prior in self.fPriors:
                self.RegularizationComparisonForMethod(method, prior)

        if self.fTruthSpectrum:
            ratioTruthOverMeas = self.fTruthSpectrum.Clone("{0}_TruthOverMeasured".format(self.fInputSpectrum.GetName()))
            cnameTOM = "Unfolding_TruthOverMeasured_{0}".format(self.fInputSpectrum.GetName())
            ratioTruthOverMeas.GetYaxis().SetTitle("Truth / Measured")
            globalList.append(ratioTruthOverMeas)
            cTOM = ROOT.TCanvas(cnameTOM, cnameTOM)
            cTOM.cd()
            globalList.append(cTOM)
            ratioTruthOverMeas.Divide(self.fInputSpectrum)
            ratioTruthOverMeas.SetFillColor(ROOT.kRed-10)
            ratioTruthOverMeas.SetFillStyle(1001)
            ratioTruthOverMeas.SetLineColor(ROOT.kRed+2)
            ratioTruthOverMeas.SetMarkerColor(ROOT.kRed+2)
            ratioTruthOverMeas.SetMarkerStyle(ROOT.kFullCircle)
            ratioTruthOverMeas.SetMarkerSize(0.9)
            ratioTruthOverMeas.Draw()

    def RegularizationComparisonForMethod(self, method, prior):
        if self.fTruthSpectrum:
            baselineId = None
            baseline = self.fTruthSpectrum.Clone("{0}_copy".format(self.fTruthSpectrum.GetName()))
            baseline.SetTitle("MC Truth")
        else:
            baselineId = (method, self.UnfoldingConfig[method]["default_reg"], prior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
            baseline.SetTitle("Reg = {0}".format(self.UnfoldingConfig[method]["default_reg"]))
        globalList.append(baseline)
        spectra = []
        for (methodIt, regIt, priorIt), unfolded in self.fUnfoldedSpectra.iteritems():
            if not (methodIt == method and  priorIt == prior):
                continue
            if baselineId == (methodIt, regIt, priorIt):
                continue
            spectrum = unfolded.Clone("{0}_copy".format(unfolded.GetName()))
            spectrum.SetTitle("Reg = {0}".format(regIt))
            spectra.append(spectrum)
            globalList.append(spectrum)
        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "UnfoldingRegularization_{0}_{1}_{2}".format(method, prior, self.fInputSpectrum.GetName()))
        for obj in r:
            globalList.append(obj)

    def UnfoldingSummaryForMethod(self, method, reg, prior):
        cname = "Unfolding_{0}_{1}_{2}_{3}".format(self.fInputSpectrum.GetName(), method, reg, prior)
        c = ROOT.TCanvas(cname, cname)
        c.cd()
        c.SetLogy()
        globalList.append(c)

        leg = ROOT.TLegend(0.55, 0.70, 0.85, 0.87)
        globalList.append(leg)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(43)
        leg.SetTextSize(16)

        if self.fTruthSpectrum:
            truth = self.fTruthSpectrum.DrawCopy("hist")
            globalList.append(truth)
            truth.SetName("{0}_{1}_{2}_{3}".format(self.fUnfoldedSpectra[(method, reg, prior)].GetName(), method, reg, prior))
            truth.SetMarkerStyle(ROOT.kOpenSquare)
            truth.SetMarkerColor(ROOT.kOrange+2)
            truth.SetMarkerSize(1.1)
            truth.SetFillColor(ROOT.kOrange-9)
            truth.SetFillStyle(1001)
            truth.SetLineColor(ROOT.kOrange+2)
            leg.AddEntry(truth, "Truth", "f")
        else:
            truth = None

        if truth:
            meas = self.fInputSpectrum.DrawCopy("same")
        else:
            meas = self.fInputSpectrum.DrawCopy()
        globalList.append(meas)
        meas.SetName("{0}_{1}_{2}_{3}".format(meas, method, reg, prior))
        meas.SetMarkerStyle(ROOT.kFullCircle)
        meas.SetMarkerColor(ROOT.kRed+2)
        meas.SetMarkerSize(0.9)
        meas.SetLineColor(ROOT.kRed+2)
        leg.AddEntry(meas, "Measured", "pe")

        unfolded = self.fUnfoldedSpectra[(method, reg, prior)].DrawCopy("same")
        globalList.append(unfolded)
        unfolded.SetName("{0}_{1}_{2}_{3}".format(unfolded.GetName(), method, reg, prior))
        unfolded.SetMarkerStyle(ROOT.kFullSquare)
        unfolded.SetMarkerColor(ROOT.kBlue+2)
        unfolded.SetMarkerSize(0.9)
        unfolded.SetLineColor(ROOT.kBlue+2)
        leg.AddEntry(unfolded, "{0}, reg={1}, prior={2}".format(method, reg, prior), "pe")

        refolded = self.fRefoldedSpectra[(method, reg, prior)].DrawCopy("same")
        globalList.append(refolded)
        refolded.SetName("{0}_{1}_{2}_{3}".format(refolded.GetName(), method, reg, prior))
        refolded.SetMarkerStyle(ROOT.kOpenCircle)
        refolded.SetMarkerColor(ROOT.kGreen+2)
        refolded.SetMarkerSize(1.1)
        refolded.SetLineColor(ROOT.kGreen+2)
        leg.AddEntry(refolded, "Refolded", "l")

        leg.Draw()

        #ratios
        # unfolded/measured
        ratioUnfoldedOverMeas = unfolded.Clone("{0}_UnfoldedOverMeasured".format(meas.GetName()))
        cnameUOM = "Unfolding_UnfoldedOverMeasured_{0}_{1}_{2}_{3}".format(self.fInputSpectrum.GetName(), method, reg, prior)
        ratioUnfoldedOverMeas.GetYaxis().SetTitle("Unfolded / Measured")
        globalList.append(ratioUnfoldedOverMeas)
        cUOM = ROOT.TCanvas(cnameUOM, cnameUOM)
        cUOM.cd()
        globalList.append(cUOM)
        ratioUnfoldedOverMeas.Divide(meas)
        ratioUnfoldedOverMeas.Draw()

        # refolded/measured and unfolded/truth 
        cnameROM = "Unfolding_RefoldedOverMeasured_{0}_{1}_{2}_{3}".format(self.fInputSpectrum.GetName(), method, reg, prior)
        cROM = ROOT.TCanvas(cnameROM, cnameROM)
        cROM.cd()
        globalList.append(cROM)
        ratioRefoldedOverMeas = refolded.Clone("{0}_RefoldedOverMeasured".format(meas.GetName()))
        globalList.append(ratioRefoldedOverMeas)
        ratioRefoldedOverMeas.Divide(meas)
        ratioRefoldedOverMeas.Draw("hist")
        if truth:
            ratioRefoldedOverMeas.GetYaxis().SetTitle("Ratio")
            ratioUnfoldedOverTruth = unfolded.Clone("{0}_UnfoldedOverTruth".format(meas.GetName()))
            globalList.append(ratioUnfoldedOverTruth)
            ratioUnfoldedOverTruth.Divide(truth)
            minY = min(ratioUnfoldedOverTruth.GetMinimum(), ratioRefoldedOverMeas.GetMinimum()) * 0.8
            maxY = max(ratioUnfoldedOverTruth.GetMaximum(), ratioRefoldedOverMeas.GetMaximum()) * 1.2
            ratioRefoldedOverMeas.GetYaxis().SetRangeUser(minY, maxY)
            ratioUnfoldedOverTruth.Draw("same")
            legROM = ROOT.TLegend(0.55, 0.78, 0.85, 0.87)
            globalList.append(legROM)
            legROM.SetFillStyle(0)
            legROM.SetBorderSize(0)
            legROM.SetTextFont(43)
            legROM.SetTextSize(16)
            legROM.AddEntry(ratioRefoldedOverMeas, "Refolded / Measured", "l")
            legROM.AddEntry(ratioUnfoldedOverTruth, "Unfolded / Truth", "pe")
            legROM.Draw()
        else:
            ratioRefoldedOverMeas.GetYaxis().SetTitle("Refolded / Measured")

    def CompareMethods(self):
        if self.fTruthSpectrum:
            baselineId = None
            baseline = self.fTruthSpectrum.Clone("{0}_copy".format(self.fTruthSpectrum.GetName()))
            baseline.SetTitle("MC Truth")
        else:
            baselineId = (self.DefaultMethod, self.UnfoldingConfig[self.DefaultMethod]["default_reg"], self.DefaultPrior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
        globalList.append(baseline)
        spectra = []
        for method, config in self.UnfoldingConfig.iteritems():
            id = (method, config["default_reg"], self.DefaultPrior)
            if id == baselineId:
                continue
            h = self.fUnfoldedSpectra[id].Clone("{0}_copy".format(self.fUnfoldedSpectra[id].GetName()))
            spectra.append(h)
            globalList.append(h)

        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "UnfoldingMethod_{0}".format(self.fInputSpectrum.GetName()))
        for obj in r:
            globalList.append(obj)

    def UnfoldSvd(self):
        max_reg = self.fInputSpectrum.GetNbinsX()
        for prior, resp in self.fResponseMatrices.iteritems():
            for reg in range(1, max_reg+1):
                unfold = ROOT.RooUnfoldSvd(resp, self.fInputSpectrum, reg)
                unfolded = unfold.Hreco()
                if not prior in self.fSvdDvectors.keys():
                    dvector = unfold.Impl().GetD().Clone("SvdDvector_{0}_{1}".format(self.fInputSpectrum.GetName(), prior))
                    dvector.GetXaxis().SetTitle("k")
                    dvector.GetYaxis().SetTitle("d")
                    self.fSvdDvectors[prior] = dvector
                cov = ROOT.TH2D(unfold.Ereco())
                cov.SetName("CovMat_{0}_{1}_{2}_{3}".format(self.fInputSpectrum.GetName(), "Svd", reg, prior))
                cov.GetXaxis().SetTitle("bin number")
                cov.GetYaxis().SetTitle("bin number")
                self.fCovarianceMatrices["Svd", reg, prior] = cov
                unfolded.SetName("{0}_{1}_{2}_{3}".format(self.fInputSpectrum.GetName(), "Svd", reg, prior))
                unfolded.SetTitle("{0}, k={1}, prior={2}".format("Svd", reg, prior))
                self.fUnfoldedSpectra["Svd", reg, prior] = unfolded
                refolded = resp.ApplyToTruth(unfolded)
                refolded.SetName("{0}_{1}_{2}_{3}_Refolded".format(self.fInputSpectrum.GetName(), "Svd", reg, prior))
                refolded.SetTitle("{0}, k={1}, prior={2}".format("Svd", reg, prior))
                self.fRefoldedSpectra["Svd", reg, prior] = refolded
    
    def UnfoldBayes(self, iter_min, iter_max, iter_step):
        for prior, resp in self.fResponseMatrices.iteritems():
            for reg in range(iter_min, iter_max, iter_step):
                unfold = ROOT.RooUnfoldBayes(resp, self.fInputSpectrum, reg)
                unfolded = unfold.Hreco()
                cov = ROOT.TH2D(unfold.Ereco())
                cov.SetName("CovMat_{0}_{1}_{2}_{3}".format(self.fInputSpectrum.GetName(), "Bayes", reg, prior))
                cov.GetXaxis().SetTitle("bin number")
                cov.GetYaxis().SetTitle("bin number")
                self.fCovarianceMatrices["Bayes", reg, prior] = cov
                unfolded.SetName("{0}_{1}_{2}_{3}".format(self.fInputSpectrum.GetName(), "Bayes", reg, prior))
                unfolded.SetTitle("{0}, iter={1}, prior={2}".format("Bayes", reg, prior))
                self.fUnfoldedSpectra["Bayes", reg, prior] = unfolded
                refolded = resp.ApplyToTruth(unfolded)
                refolded.SetName("{0}_{1}_{2}_{3}_Refolded".format(self.fInputSpectrum.GetName(), "Svd", reg, prior))
                refolded.SetTitle("{0}, k={1}, prior={2}".format("Svd", reg, prior))
                self.fRefoldedSpectra["Bayes", reg, prior] = refolded

    def GenerateRooUnfoldResponse(self):
        for prior in self.fPriors:
            if prior == "train_truth":
                rooUnfoldResp = ROOT.RooUnfoldResponse(ROOT.nullptr, self.fDetectorTrainTruth, self.fDetectorResponse)
                priorSpectrum = self.fDetectorTrainTruth.Clone()
                
            else:
                (detResp, trainTruth) = self.NormalizeResponseMatrix(self.fDetectorResponse, prior, self.fDetectorEfficiency)
                rooUnfoldResp = ROOT.RooUnfoldResponse(ROOT.nullptr, trainTruth, detResp)
                priorSpectrum = trainTruth.Clone()
            priorSpectrum.Scale(1. / priorSpectrum.Integral())
            priorSpectrum.SetName("Prior_{0}".format(prior))
            self.fPriorSpectra[prior] = priorSpectrum
            self.fResponseMatrices[prior] = rooUnfoldResp

class DMesonJetUnfolding:
    def __init__(self, name, input_path, dataTrain, data, responseTrain, response):
        self.fName = name
        self.fInputPath = input_path
        self.fDataTrain = dataTrain
        self.fData = data
        self.fResponseTrain = responseTrain
        self.fResponse = response
        self.fUnfoldingEngine = []
        self.OpenFiles()

    def OpenFiles(self):
        self.fDataFile = ROOT.TFile("{0}/{1}/{2}.root".format(self.fInputPath, self.fDataTrain, self.fData))
        self.fResponseFile = ROOT.TFile("{0}/{1}/{2}.root".format(self.fInputPath, self.fResponseTrain, self.fResponse))

    def StartUnfolding(self, config):
        eng = DMesonJetUnfoldingEngine(config)
        self.fUnfoldingEngine.append(eng)
        r = eng.LoadData(self.fDataFile, self.fResponseFile)
        if r:
            eng.Start()
