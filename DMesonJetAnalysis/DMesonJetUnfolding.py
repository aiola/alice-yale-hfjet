#!/usr/bin/env python
#python program to perform a D meson jet unfolding

import ROOT
import math
import DMesonJetProjectors
from DMesonJetBase import *
import array
import copy
import DMesonJetUtils
import os
from collections import OrderedDict

globalList = []

class ResponseMatrix:
    def __init__(self, name, rooUnfoldResponse, response, truth, smallBinResponse, smallBinTruth, normResp):
        self.fName = name
        self.fRooUnfoldResponse = rooUnfoldResponse

        self.fResponse = response.Clone("{0}_ResponseMatrix".format(self.fName))
        self.fResponse.SetTitle("{0} Response Matrix".format(self.fName))

        self.fNormResponse = normResp.Clone("{0}_NormalizedResponseMatrix".format(self.fName))
        self.fNormResponse.SetTitle("{0} Normalized Response Matrix".format(self.fName))

        self.fTruth = truth.Clone("{0}_Truth".format(self.fName))
        self.fTruth.SetTitle("{0} Truth".format(self.fName))

        self.fSmallBinResponse = smallBinResponse.Clone("{0}_SmallBinResponseMatrix".format(self.fName))
        self.fSmallBinResponse.SetTitle("{0} Small-Bin Response Matrix".format(self.fName))

        self.fSmallBinTruth = smallBinTruth.Clone("{0}_SmallBinTruth".format(self.fName))
        self.fSmallBinTruth.SetTitle("{0} Small-Bin Truth".format(self.fName))

        self.GenerateEfficiencies()

        self.GenerateResponseUncertainty()

    def GenerateResponseUncertainty(self):
        self.fResponseUncertainty = ROOT.TH2D("{0}_ResponseMatrixUncertainty".format(self.fName), "{0} Response Matrix Uncertainty".format(self.fName), self.fResponse.GetNbinsX(), self.fResponse.GetXaxis().GetXbins().GetArray(), self.fResponse.GetNbinsY(), self.fResponse.GetYaxis().GetXbins().GetArray())
        self.fResponseUncertainty.GetXaxis().SetTitle(self.fResponse.GetXaxis().GetTitle())
        self.fResponseUncertainty.GetYaxis().SetTitle(self.fResponse.GetYaxis().GetTitle())
        self.fResponseUncertainty.GetZaxis().SetTitle("relative statistical uncertainty")
        for x in range(0, self.fResponseUncertainty.GetNbinsX()+2):
            for y in range(0, self.fResponseUncertainty.GetNbinsY()+2):
                self.fResponseUncertainty.SetBinContent(x,y,self.fResponse.GetBinError(x,y) / self.fResponse.GetBinContent(x,y))

    def GenerateEfficiencies(self):
        # total efficiency
        self.fTotalEfficiency =  self.fResponse.ProjectionY("{0}_TotalEfficiency".format(self.fName), 1, self.fResponse.GetXaxis().GetNbins())
        self.fTotalEfficiency.Divide(self.fTruth)
        self.fTotalEfficiency.SetTitle("Total Efficiency")
        self.fTotalEfficiency.GetYaxis().SetTitle("Total Efficiency")

        # reconstruction efficiency
        self.fRecoEfficiency =  self.fResponse.ProjectionY("{0}_RecoEfficiency".format(self.fName), 0, -1)
        self.fRecoEfficiency.Divide(self.fTruth)
        self.fRecoEfficiency.SetTitle("Reconstruction Efficiency")
        self.fRecoEfficiency.GetYaxis().SetTitle("Reconstruction Efficiency")
        
        # kinematic efficiency
        self.fKineEfficiency =  self.fResponse.ProjectionY("{0}_KineEfficiency".format(self.fName), 1, self.fResponse.GetXaxis().GetNbins())
        self.fKineEfficiency.Divide(self.fResponse.ProjectionY("temp", 0, -1))
        self.fKineEfficiency.SetTitle("Kinematic Efficiency")
        self.fKineEfficiency.GetYaxis().SetTitle("Kinematic Efficiency")

    def GenerateRootList(self):
        rlist = ROOT.TList()
        rlist.SetName("Response_{0}".format(self.fName))
        rlist.Add(self.fResponse)
        rlist.Add(self.fNormResponse)
        rlist.Add(self.fResponseUncertainty)
        rlist.Add(self.fTruth)
        rlist.Add(self.fSmallBinResponse)
        rlist.Add(self.fSmallBinTruth)
        rlist.Add(self.fTotalEfficiency)
        rlist.Add(self.fRecoEfficiency)
        rlist.Add(self.fKineEfficiency)
        return rlist

class DMesonJetUnfoldingEngine:
    def __init__(self, config):
        self.fDMeson = config["d_meson"]
        self.fDMesonTruth = config["d_meson_truth"]
        self.fDMesonResponse = config["d_meson_response"]
        self.fJetDefinition = config["jet"]
        self.fSpectrumName = config["spectrum"]
        self.fSpectrumResponseName = config["spectrum_response"]
        self.fPriors = config["priors"]
        self.fUnfoldingConfig = config["methods"]
        self.fDefaultPrior = config["default_prior"]
        self.fDefaultMethod = config["default_method"]
        self.fName = config["name"]
        self.fResponseMatrices = OrderedDict()
        self.fUnfoldedSpectra = OrderedDict()
        self.fUnfoldedSpectraErrors = OrderedDict()
        self.fRefoldedSpectra = OrderedDict()
        self.fCovarianceMatrices = OrderedDict()
        self.fSvdDvectors = OrderedDict()
        self.fPearsonMatrices = OrderedDict()
        self.fBinByBinCorrectionFactors = OrderedDict()
        self.fCanvases = []

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))

    def GenerateRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        rlist.Add(self.fInputSpectrum)
        if self.fTruthSpectrum:
            rlist.Add(self.fTruthSpectrum)
        rlist.Add(self.fDetectorResponse)
        rlist.Add(self.fDetectorTrainTruth)
        for r in self.fResponseMatrices.itervalues():
            rlist.Add(r.GenerateRootList())
        unfoldingHistos = dict()
        for method in self.fUnfoldingConfig.iterkeys():
            rMethodList = ROOT.TList()
            rMethodList.SetName(method)
            rlist.Add(rMethodList)
            unfoldingHistos[method] = rMethodList
        for (method,reg,prior),h in self.fUnfoldedSpectra.iteritems():
            unfoldingHistos[method].Add(h)            
        for (method,reg,prior),h in self.fRefoldedSpectra.iteritems():
            unfoldingHistos[method].Add(h)
        for (method,reg,prior),h in self.fCovarianceMatrices.iteritems():
            unfoldingHistos[method].Add(h)
        for h in self.fSvdDvectors.itervalues():
            unfoldingHistos["Svd"].Add(h)
        for (method,reg,prior),h in self.fPearsonMatrices.iteritems():
            unfoldingHistos[method].Add(h)
        for (method,reg,prior),(h1,h2,h3) in self.fUnfoldedSpectraErrors.iteritems():
            unfoldingHistos[method].Add(h1)
            unfoldingHistos[method].Add(h2)
            unfoldingHistos[method].Add(h3)    
        for h in self.fBinByBinCorrectionFactors.itervalues():
            unfoldingHistos["BinByBin"].Add(h)
        return rlist

    def LoadData(self, dataFile, responseFile, eff):
        print("Now loading data")
        dataListName = "{0}_{1}".format(self.fDMeson, self.fSpectrumName)
        dataList = dataFile.Get(dataListName)
        if not dataList:
            print("Could not find list {0} in file {1}". format(dataListName, dataFile.GetName()))
            return False
        inputSpectrumName = "{0}_{1}".format(self.fDMeson, self.fSpectrumName)
        inputSpectrum = dataList.FindObject(inputSpectrumName)
        if not inputSpectrum:
            print("Could not find histogrm {0} in list {1} in file {2}". format(inputSpectrumName, dataListName, dataFile.GetName()))
            return False
        self.fInputSpectrum = inputSpectrum.Clone("{0}_InputSpectrum".format(self.fName))
        self.fInputSpectrum.SetTitle("{0} Input Spectrum".format(self.fName))

        if self.fDMesonTruth:
            dataTruthListName = "{0}_{1}".format(self.fDMesonTruth, self.fSpectrumName)
            dataTruthList = dataFile.Get(dataTruthListName)
            if not dataTruthList:
                print("Could not find list {0} in file {1}". format(dataTruthListName, dataFile.GetName()))
                return False
            truthSpectrumName = "{0}_{1}".format(self.fDMesonTruth, self.fSpectrumName)
            truthSpectrum = dataTruthList.FindObject(truthSpectrumName)
            if not truthSpectrum:
                print("Could not find histogrm {0} in list {1} in file {2}". format(truthSpectrumName, dataTruthListName, dataFile.GetName()))
                return False
            self.fTruthSpectrum = truthSpectrum.Clone("{0}_TruthSpectrum".format(self.fName))
            self.fTruthSpectrum.SetTitle("{0} Truth Spectrum".format(self.fName))
        else:
            self.fTruthSpectrum = None

        responseListName = "{0}_{1}_{2}".format(self.fDMesonResponse, self.fJetDefinition, self.fSpectrumResponseName)
        responseList = responseFile.Get(responseListName)
        if not responseList:
            print("Could not find list {0} in file {1}". format(responseListName, responseFile.GetName()))
            return False
        detectorResponseName = "{0}_{1}_{2}_DetectorResponse".format(self.fDMesonResponse, self.fJetDefinition, self.fSpectrumResponseName)
        detectorResponse = responseList.FindObject(detectorResponseName)
        if not detectorResponse:
            print("Could not find histogrm {0} in list {1} in file {2}". format(detectorResponseName, responseListName, responseFile.GetName()))
            return False
        self.fDetectorResponse = detectorResponse.Clone("{0}_DetectorResponse".format(self.fName))
        self.fDetectorResponse.SetTitle("{0} Detector Response".format(self.fName))
        if eff:
            detTrainTruthName = "{0}_{1}_{2}_Truth".format(self.fDMesonResponse, self.fJetDefinition, self.fSpectrumResponseName)
            detTrainTruth = responseList.FindObject(detTrainTruthName)
            if not detTrainTruth:
                print("Could not find histogrm {0} in list {1} in file {2}". format(detTrainTruthName, responseListName, responseFile.GetName()))
                return False
            self.fDetectorTrainTruth = detTrainTruth.Clone("{0}_ResponseTruth".format(self.fName))
        else:
            self.fDetectorTrainTruth = self.fDetectorResponse.ProjectionY("{0}_ResponseTruth".format(self.fName), 0, -1)
        self.fDetectorTrainTruth.SetTitle("{0} Response Truth".format(self.fName))
        return True

    def Start(self):
        self.GenerateResponse()
        for method, config in self.fUnfoldingConfig.iteritems():
            if method == "Svd":
                self.UnfoldSvd()
            elif method == "Bayes":
                self.UnfoldBayes(config["iter_min"], config["iter_max"], config["iter_step"])
            elif method == "BinByBin":
                self.UnfoldBinByBin()
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
        #ROOT.gStyle.SetPalette(55,ROOT.nullptr)
        if len(self.fUnfoldingConfig) > 1:
            self.CompareMethods()
        if len(self.fPriors) > 1:
            self.ComparePriors()
        self.UnfoldingSummary()
        self.PlotPearsonMatrices()
        self.PlotSvdDvectors()
        self.PlotResponses()
        self.PlotPriors()
        self.PlotBinByBinCorrectionFactors()
        self.PlotUnfoldingErrors()

    def PlotBinByBinCorrectionFactors(self):
        baselineId = self.fDefaultPrior
        baseline = self.fBinByBinCorrectionFactors[baselineId].Clone("{0}_copy".format(self.fBinByBinCorrectionFactors[baselineId].GetName()))
        baseline.SetTitle(self.fDefaultPrior)
        yaxisRatio = "Prior x / {0}".format(self.fDefaultPrior)
        globalList.append(baseline)
        spectra = []
        for priorIt, corrFact in self.fBinByBinCorrectionFactors.iteritems():
            if baselineId == priorIt:
                continue
            cf = corrFact.Clone("{0}_copy".format(corrFact.GetName()))
            cf.SetTitle(priorIt)
            spectra.append(cf)
            globalList.append(cf)
        if len(spectra) < 1:
            return
        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "{0}_BinByBinCorrectionFactors".format(self.fName), "hist", "hist", yaxisRatio, "lineary", "lineary")
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def PlotPriors(self):
        default_reg = self.GetDefaultRegularization(self.fDefaultMethod, self.fDefaultPrior)
        default_unfolded = self.fUnfoldedSpectra[self.fDefaultMethod, default_reg, self.fDefaultPrior]
        if self.fTruthSpectrum:
            baselineId = None
            baseline = self.fTruthSpectrum.Clone("{0}_copy".format(self.fTruthSpectrum.GetName()))
            baseline.SetTitle("MC Truth")
            yaxisRatio = "Prior / Truth"
        elif default_unfolded:
            baselineId = None
            baseline = default_unfolded.Clone("{0}_copy".format(default_unfolded.GetName()))
            yaxisRatio = "Prior / Unfolded"
        else:
            baselineId = self.fDefaultPrior
            baseline = self.fResponseMatrices[baselineId].fTruth.Clone("{0}_copy".format(self.fResponseMatrices[baselineId].fTruth.GetName()))
            baseline.SetTitle(self.fDefaultPrior)
            yaxisRatio = "Prior x / {0}".format(self.fDefaultPrior)
        baseline.Scale(1. / baseline.Integral())
        globalList.append(baseline)
        spectra = []
        for priorIt, resp in self.fResponseMatrices.iteritems():
            if baselineId == priorIt:
                continue
            spectrum = resp.fTruth.Clone("{0}_copy".format(resp.fTruth.GetName()))
            spectrum.Scale(1. / spectrum.Integral())
            spectrum.SetTitle(priorIt)
            spectra.append(spectrum)
            globalList.append(spectrum)
        if len(spectra) < 1:
            return
        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "{0}_Priors".format(self.fName), "hist", "hist", yaxisRatio, "logy", "logy")
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def PlotResponses(self):
        for name, resp in self.fResponseMatrices.iteritems():
            cname = "{0}_Response_Prior{1}".format(self.fName, name)
            c = ROOT.TCanvas(cname, cname, 850, 750)
            self.fCanvases.append(c)
            globalList.append(c)
            c.Divide(2,2)

            pad = c.cd(1)
            pad.SetLogz()
            pad.SetLeftMargin(0.16)
            pad.SetRightMargin(0.24)
            pad.SetBottomMargin(0.15)
            rm = resp.fNormResponse.DrawCopy("col2z")
            rm.GetXaxis().SetTitleFont(43)
            rm.GetXaxis().SetTitleSize(19)
            rm.GetXaxis().SetTitleOffset(2.6)
            rm.GetXaxis().SetLabelFont(43)
            rm.GetXaxis().SetLabelSize(16)
            rm.GetYaxis().SetTitleFont(43)
            rm.GetYaxis().SetTitleSize(19)
            rm.GetYaxis().SetTitleOffset(2.6)
            rm.GetYaxis().SetLabelFont(43)
            rm.GetYaxis().SetLabelSize(16)
            rm.GetZaxis().SetTitleFont(43)
            rm.GetZaxis().SetTitleSize(16)
            rm.GetZaxis().SetTitleOffset(2.8)
            rm.GetZaxis().SetLabelFont(43)
            rm.GetZaxis().SetLabelSize(15)
            rm.SetMaximum(1)
            globalList.append(rm)

            pad = c.cd(2)
            pad.SetLeftMargin(0.14)
            pad.SetBottomMargin(0.15)
            self.PlotEfficiency(resp.fTotalEfficiency)

            pad = c.cd(3)
            pad.SetLeftMargin(0.14)
            pad.SetBottomMargin(0.15)
            self.PlotEfficiency(resp.fRecoEfficiency)

            pad = c.cd(4)
            pad.SetLeftMargin(0.14)
            pad.SetBottomMargin(0.15)
            self.PlotEfficiency(resp.fKineEfficiency)

    def PlotEfficiency(self, eff):
        eff_copy = eff.DrawCopy("hist")
        if eff_copy.GetMinimum() > 0.99: # fix ugly plotting when the efficiency is 1
            eff_copy.SetMaximum(1.2)
            eff_copy.SetMinimum(0.8)
        eff_copy.GetXaxis().SetTitleFont(43)
        eff_copy.GetXaxis().SetTitleSize(19)
        eff_copy.GetXaxis().SetTitleOffset(2.6)
        eff_copy.GetXaxis().SetLabelFont(43)
        eff_copy.GetXaxis().SetLabelSize(16)
        eff_copy.GetYaxis().SetTitleFont(43)
        eff_copy.GetYaxis().SetTitleSize(19)
        eff_copy.GetYaxis().SetTitleOffset(2.6)
        eff_copy.GetYaxis().SetLabelFont(43)
        eff_copy.GetYaxis().SetLabelSize(16)
        globalList.append(eff_copy)

    def GetDefaultRegularization(self, method, prior):
        if isinstance(self.fUnfoldingConfig[method]["default_reg"], dict):
            reg = self.fUnfoldingConfig[method]["default_reg"][prior]
        else:
            reg = self.fUnfoldingConfig[method]["default_reg"]
        return reg

    def FilterObjectsFromDict(self, objects, method, reg, prior):
        for (methodIt, regIt, priorIt), cov in objects.iteritems():
            if method and not methodIt == method:
                continue
            if prior and not priorIt == prior:
                continue
            if reg and not regIt == reg:
                continue
            yield (methodIt, regIt, priorIt), cov

    def PlotUnfoldingErrors(self):
        for method, config in self.fUnfoldingConfig.iteritems():
            for prior in self.fPriors:
                self.PlotUnfoldingErrorsForMethod(method, prior)
        self.PlotUnfoldingErrorsMethodComparison()

    def PlotUnfoldingErrorsMethodComparison(self):
        reg = self.GetDefaultRegularization(self.fDefaultMethod, self.fDefaultPrior)
        baselineId = (self.fDefaultMethod, reg, self.fDefaultPrior)
        (baselineErr, baselineCov, baselineToyMC) = self.fUnfoldedSpectraErrors[baselineId]
        baseline = baselineErr.Clone("{0}_copy".format(baselineErr.GetName()))
        baseline.SetTitle(self.fDefaultMethod)
        globalList.append(baseline)

        spectra = []
        for method in self.fUnfoldingConfig.iterkeys():
            if method == self.fDefaultMethod:
                continue
            reg = self.GetDefaultRegularization(method, self.fDefaultPrior)
            id = (method, reg, self.fDefaultPrior)
            (hErr, hCov, hToyMC) = self.fUnfoldedSpectraErrors[id]
            h = hErr.Clone("{0}_copy".format(hErr.GetName()))
            h.SetTitle(method)
            spectra.append(h)
            globalList.append(h)

        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "{0}_UnfoldingStatisticalUncertaintyCompareMethods".format(self.fName), "hist", "", "", "lineary", False)
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def PlotUnfoldingErrorsForMethod(self, method, prior):
        #first compare different regularization strengths
        spectra = []
        for (methodIt, regIt, priorIt), unfolded in self.fUnfoldedSpectra.iteritems():
            if not (methodIt == method and  priorIt == prior):
                continue
            id = (method, regIt, prior)
            (hErr, hCov, hToyMC) = self.fUnfoldedSpectraErrors[id]
            h = hErr.Clone("{0}_copy".format(hErr.GetName()))
            h.SetTitle("Reg={0}".format(regIt))
            spectra.append(h)
            globalList.append(h)

        if len(spectra) > 1:
            r = DMesonJetUtils.CompareSpectra(spectra[0], spectra[1:], "{0}_UnfoldingStatisticalUncertainty_{1}_Prior{2}".format(self.fName, method, prior), "hist", "", "", "lineary", False)
            for obj in r:
                globalList.append(obj)
                if isinstance(obj, ROOT.TCanvas):
                    self.fCanvases.append(obj)

        #now compare different error estimation methods
        reg = self.GetDefaultRegularization(method, prior)
        baselineId = (method, reg, prior)
        (baselineErr, baselineCov, baselineToyMC) = self.fUnfoldedSpectraErrors[baselineId]
        baseline = baselineErr.Clone("{0}_copy".format(baselineErr.GetName()))
        baseline.SetTitle("kError")
        globalList.append(baseline)

        spectra = []
        hcov = baselineCov.Clone("{0}_copy".format(baselineCov.GetName()))
        hcov.SetTitle("kCovariance")
        globalList.append(hcov)
        spectra.append(hcov)
        htoy = baselineToyMC.Clone("{0}_copy".format(baselineToyMC.GetName()))
        htoy.SetTitle("kCovToy")
        globalList.append(htoy)
        spectra.append(htoy)

        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "{0}_UnfoldingStatisticalUncertaintyStrategy_{1}_Reg{2}_Prior{3}".format(self.fName, method, reg, prior), "hist", "", "", "lineary", False)
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def PlotPearsonMatrices(self):
        for method, config in self.fUnfoldingConfig.iteritems():
            for prior in self.fPriors:
                self.PlotPearsonMatricesForMethod(method, prior)

    def PlotPearsonMatricesForMethod(self, method, prior):
        pearsons = []
        filteredPearsons = self.FilterObjectsFromDict(self.fPearsonMatrices, method, None, prior)
        for (methodIt, reg, priorIt), cov in filteredPearsons:
            cov_copy = cov.Clone()
            if reg:
                cov_copy.SetTitle("reg={0}".format(reg))
            else:
                cov_copy.SetTitle("No regularization".format(reg))
            pearsons.append(cov_copy)
            globalList.append(cov_copy)
        if len(pearsons) < 1:
            return

        cname = "{0}_Pearson_{1}_Prior{2}".format(self.fName, method, prior)
        rows = int(math.floor(math.sqrt(len(pearsons))))
        cols = int(math.ceil(len(pearsons) / rows))
        c = ROOT.TCanvas(cname, cname, 360*cols, 350*rows)
        self.fCanvases.append(c)
        globalList.append(c)
        c.Divide(cols, rows)

        for i, cov in enumerate(pearsons):
            c.cd(i+1)
            cov.GetXaxis().SetTitleFont(43)
            cov.GetXaxis().SetTitleSize(16)
            cov.GetXaxis().SetTitleOffset(2)
            cov.GetXaxis().SetLabelFont(43)
            cov.GetXaxis().SetLabelSize(15)
            cov.GetYaxis().SetTitleFont(43)
            cov.GetYaxis().SetTitleSize(16)
            cov.GetYaxis().SetTitleOffset(2)
            cov.GetYaxis().SetLabelFont(43)
            cov.GetYaxis().SetLabelSize(15)
            if i % cols == cols-1:
                cov.Draw("col2z")
            else:
                cov.Draw("col2")
            cov.GetZaxis().SetRangeUser(-1,1)
            pave = ROOT.TPaveText(0, 1, 1, 0.9, "NB NDC")
            pave.SetBorderSize(0)
            pave.SetFillStyle(0)
            pave.SetTextFont(43)
            pave.SetTextSize(17)
            pave.SetTextAlign(22)
            pave.AddText(cov.GetTitle())
            pave.Draw()
            globalList.append(pave)

    def PlotSvdDvectors(self):
        for prior, d in self.fSvdDvectors.iteritems():
            cname = d.GetName()
            c = ROOT.TCanvas(cname, cname)
            self.fCanvases.append(c)
            c.cd()
            c.SetLogy()
            globalList.append(c)
            d_copy = d.DrawCopy()
            globalList.append(d_copy)

    def UnfoldingSummary(self):
        for method, config in self.fUnfoldingConfig.iteritems():
            reg = self.GetDefaultRegularization(method, self.fDefaultPrior)
            self.UnfoldingSummaryForMethod(method, reg, self.fDefaultPrior)
            for prior in self.fPriors:
                self.RegularizationComparisonForMethod(method, prior)

        if self.fTruthSpectrum:
            ratioTruthOverMeas = self.fTruthSpectrum.Clone("{0}_TruthOverMeasured".format(self.fName))
            cnameTOM = "{0}_UnfoldingSummary_TruthOverMeasured".format(self.fName)
            ratioTruthOverMeas.GetYaxis().SetTitle("Truth / Measured")
            globalList.append(ratioTruthOverMeas)
            cTOM = ROOT.TCanvas(cnameTOM, cnameTOM)
            self.fCanvases.append(cTOM)
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
            yaxisRatio = "Unfolded / Truth"
        else:
            reg = self.GetDefaultRegularization(method, prior)
            baselineId = (method, reg, prior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
            baseline.SetTitle("Reg = {0}".format(reg))
            yaxisRatio = "Reg = x / Reg = {0}".format(reg)
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
        if (self.fTruthSpectrum and len(spectra) < 2) or len(spectra) < 1:
            return
        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "{0}_UnfoldingRegularization_{1}_Prior{2}".format(self.fName, method, prior), "", "", yaxisRatio, "logy", "lineary")
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def UnfoldingSummaryForMethod(self, method, reg, prior):
        cname = "{0}_UnfoldingSummary_{1}".format(self.fName, method)
        c = ROOT.TCanvas(cname, cname)
        self.fCanvases.append(c)
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
        cnameUOM = "{0}_UnfoldingSummary_{1}_UnfoldedOverMeasured".format(self.fName, method)
        ratioUnfoldedOverMeas.GetYaxis().SetTitle("Unfolded / Measured")
        globalList.append(ratioUnfoldedOverMeas)
        cUOM = ROOT.TCanvas(cnameUOM, cnameUOM)
        self.fCanvases.append(cUOM)
        cUOM.cd()
        globalList.append(cUOM)
        ratioUnfoldedOverMeas.Divide(meas)
        ratioUnfoldedOverMeas.Draw()

        # refolded/measured and unfolded/truth 
        cnameROM = "{0}_UnfoldingSummary_{1}_RefoldedOverMeasured".format(self.fName, method)
        cROM = ROOT.TCanvas(cnameROM, cnameROM)
        self.fCanvases.append(cROM)
        cROM.cd()
        globalList.append(cROM)
        ratioRefoldedOverMeas = refolded.Clone("{0}_RefoldedOverMeasured".format(meas.GetName()))
        globalList.append(ratioRefoldedOverMeas)
        ratioRefoldedOverMeas.Divide(meas)
        ratioRefoldedOverMeas.Draw("hist")
        if truth:
            ratioRefoldedOverMeas.GetYaxis().SetTitle("ratio")
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
            yaxisRatio = "Unfolded / Truth"
        else:
            reg = self.GetDefaultRegularization(self.fDefaultMethod, self.fDefaultPrior)
            baselineId = (self.fDefaultMethod, reg, self.fDefaultPrior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
            yaxisRatio = "Other Method / {0}".format(self.fDefaultMethod)
        globalList.append(baseline)
        spectra = []
        for method, config in self.fUnfoldingConfig.iteritems():
            reg = self.GetDefaultRegularization(method, self.fDefaultPrior)
            id = (method, reg, self.fDefaultPrior)
            if id == baselineId:
                continue
            h = self.fUnfoldedSpectra[id].Clone("{0}_copy".format(self.fUnfoldedSpectra[id].GetName()))
            spectra.append(h)
            globalList.append(h)

        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "{0}_UnfoldingMethod".format(self.fName), "", "", yaxisRatio, "logy", "lineary")
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def ComparePriors(self):
        for method in self.fUnfoldingConfig.iterkeys():
            self.ComparePriorsForMethod(method)

    def ComparePriorsForMethod(self, method):
        config = self.fUnfoldingConfig[method]
        if self.fTruthSpectrum:
            baselineId = None
            baseline = self.fTruthSpectrum.Clone("{0}_copy".format(self.fTruthSpectrum.GetName()))
            baseline.SetTitle("MC Truth")
            yaxisRatio = "Unfolded / Truth"
        else:
            reg = self.GetDefaultRegularization(method, self.fDefaultPrior)
            baselineId = (method, reg, self.fDefaultPrior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
            yaxisRatio = "Prior = x / Prior = {0}".format(self.fDefaultPrior)
        globalList.append(baseline)
        spectra = []
        for prior in self.fPriors:
            reg = self.GetDefaultRegularization(method, prior)
            id = (method, reg, prior)
            if id == baselineId:
                continue
            h = self.fUnfoldedSpectra[id].Clone("{0}_copy".format(self.fUnfoldedSpectra[id].GetName()))
            spectra.append(h)
            globalList.append(h)

        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "{0}_UnfoldingPrior_{1}".format(self.fName, method), "", "", yaxisRatio, "logy", "lineary")
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def UnfoldSvd(self):
        max_reg = self.fInputSpectrum.GetNbinsX()
        for prior, resp in self.fResponseMatrices.iteritems():
            for reg in range(1, max_reg+1):
                print("Unfolding {0}, reg={1}, prior={2}".format("SVD", reg, prior))
                unfold = ROOT.RooUnfoldSvd(resp.fRooUnfoldResponse, self.fInputSpectrum, reg)
                unfolded = unfold.Hreco()
                if not prior in self.fSvdDvectors.keys():
                    dvector = unfold.Impl().GetD().Clone("{0}_SvdDvector_Prior{1}".format(self.fName, prior))
                    dvector.GetXaxis().SetTitle("k")
                    dvector.GetYaxis().SetTitle("d")
                    self.fSvdDvectors[prior] = dvector

                unfolded.SetName("{0}_UnfoldedSpectrum_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                unfolded.SetTitle("{0}, k={1}, prior={2}".format("Svd", reg, prior))
                self.fUnfoldedSpectra["Svd", reg, prior] = unfolded

                #refolded
                refolded = resp.fRooUnfoldResponse.ApplyToTruth(unfolded)
                refolded.SetName("{0}_RefoldedSpectrum_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                refolded.SetTitle("{0}, k={1}, prior={2}".format("Svd", reg, prior))
                self.fRefoldedSpectra["Svd", reg, prior] = refolded

                #covariance
                cov = ROOT.TH2D(unfold.Ereco())
                cov.SetName("{0}_CovMat_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                cov.GetXaxis().SetTitle("bin number")
                cov.GetYaxis().SetTitle("bin number")
                self.fCovarianceMatrices["Svd", reg, prior] = cov

                #errors
                unfoldingErr = ROOT.TH1D(unfold.ErecoV(ROOT.RooUnfold.kErrors))
                unfoldingErr.SetName("{0}_UnfoldErr_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                unfoldingErr.SetTitle("{0} Unfolding Errors {1} Reg={2} Prior={3}".format(self.fName, "Svd", reg, prior))
                DMesonJetUtils.DivideNoErrors(unfoldingErr, unfolded)
                diagCov = ROOT.TH1D(unfold.ErecoV(ROOT.RooUnfold.kCovariance))
                diagCov.SetName("{0}_DiagCov_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))
                diagCov.SetTitle("{0} Unfolding Errors (covariance) {1} Reg={2} Prior={3}".format(self.fName, "Svd", reg, prior))
                DMesonJetUtils.DivideNoErrors(diagCov, unfolded)
                errToy = ROOT.TH1D(unfold.ErecoV(ROOT.RooUnfold.kCovToy))
                errToy.SetName("{0}_MCToyErr_{1}_Reg{2}_Prior{3}".format(self.fName, "Svd", reg, prior))        
                errToy.SetTitle("{0} Unfolding Errors (MC-toy) {1} Reg={2} Prior={3}".format(self.fName, "Svd", reg, prior))
                DMesonJetUtils.DivideNoErrors(errToy, unfolded)
                self.fUnfoldedSpectraErrors["Svd", reg, prior] = unfoldingErr, diagCov, errToy

    def UnfoldBayes(self, iter_min, iter_max, iter_step):
        for prior, resp in self.fResponseMatrices.iteritems():
            for reg in range(iter_min, iter_max+1, iter_step):
                print("Unfolding {0}, reg={1}, prior={2}".format("Bayes", reg, prior))
                unfold = ROOT.RooUnfoldBayes(resp.fRooUnfoldResponse, self.fInputSpectrum, reg)
                unfolded = unfold.Hreco()

                unfolded.SetName("{0}_UnfoldedSpectrum_{1}_Reg{2}_Prior{3}".format(self.fName, "Bayes", reg, prior))
                unfolded.SetTitle("{0}, iter={1}, prior={2}".format("Bayes", reg, prior))
                self.fUnfoldedSpectra["Bayes", reg, prior] = unfolded

                #refolded
                refolded = resp.fRooUnfoldResponse.ApplyToTruth(unfolded)
                refolded.SetName("{0}_RefoldedSpectrum_{1}_Reg{2}_Prior{3}".format(self.fName, "Bayes", reg, prior))
                refolded.SetTitle("{0}, iter={1}, prior={2}".format("Bayes", reg, prior))
                self.fRefoldedSpectra["Bayes", reg, prior] = refolded

                #covariance
                cov = ROOT.TH2D(unfold.Ereco())
                cov.SetName("{0}_CovMat_{1}_{2}_{3}".format(self.fName, "Bayes", reg, prior))
                cov.GetXaxis().SetTitle("bin number")
                cov.GetYaxis().SetTitle("bin number")
                self.fCovarianceMatrices["Bayes", reg, prior] = cov

                #errors
                unfoldingErr = ROOT.TH1D(unfold.ErecoV(ROOT.RooUnfold.kErrors))
                unfoldingErr.SetName("{0}_UnfoldErr_{1}_Reg{2}_Prior{3}".format(self.fName, "Bayes", reg, prior))
                unfoldingErr.SetTitle("{0} Unfolding Errors {1} Reg={2} Prior={3}".format(self.fName, "Bayes", reg, prior))
                DMesonJetUtils.DivideNoErrors(unfoldingErr, unfolded)
                diagCov = ROOT.TH1D(unfold.ErecoV(ROOT.RooUnfold.kCovariance))
                diagCov.SetName("{0}_DiagCov_{1}_Reg{2}_Prior{3}".format(self.fName, "Bayes", reg, prior))
                diagCov.SetTitle("{0} Unfolding Errors (covariance) {1} Reg={2} Prior={3}".format(self.fName, "Bayes", reg, prior))
                DMesonJetUtils.DivideNoErrors(diagCov, unfolded)
                errToy = ROOT.TH1D(unfold.ErecoV(ROOT.RooUnfold.kCovToy))
                errToy.SetName("{0}_MCToyErr_{1}_Reg{2}_Prior{3}".format(self.fName, "Bayes", reg, prior))        
                errToy.SetTitle("{0} Unfolding Errors (MC-toy) {1} Reg={2} Prior={3}".format(self.fName, "Bayes", reg, prior))
                DMesonJetUtils.DivideNoErrors(errToy, unfolded)
                self.fUnfoldedSpectraErrors["Bayes", reg, prior] = unfoldingErr, diagCov, errToy

    def UnfoldBinByBin(self):
        for prior, resp in self.fResponseMatrices.iteritems():
            print("Unfolding {0}, prior={1}".format("BinByBin", prior))
            binBybinCorrFactors = resp.fResponse.ProjectionY("{0}_BinByBinCorrectionFactors_Prior{1}".format(self.fName, prior))
            binBybinCorrFactors.SetTitle("{0} Correction Factors, prior={1}".format("BinByBin", prior))
            binBybinCorrFactors.Divide(resp.fResponse.ProjectionX())
            self.fBinByBinCorrectionFactors[prior] = binBybinCorrFactors

            unfold = ROOT.RooUnfoldBinByBin(resp.fRooUnfoldResponse, self.fInputSpectrum)
            unfolded = unfold.Hreco()

            unfolded.SetName("{0}_UnfoldedSpectrum_{1}_Prior{2}".format(self.fName, "BinByBin", prior))
            unfolded.SetTitle("{0}, prior={1}".format("BinByBin", prior))
            self.fUnfoldedSpectra["BinByBin", None, prior] = unfolded

            #refolded
            refolded = resp.fRooUnfoldResponse.ApplyToTruth(unfolded)
            #refolded = unfolded.Clone("refolded")
            #refolded.Divide(binBybinCorrFactors)
            refolded.SetName("{0}_RefoldedSpectrum_{1}_Prior{2}".format(self.fName, "BinByBin", prior))
            refolded.SetTitle("{0}, prior={1}".format("BinByBin", prior))
            self.fRefoldedSpectra["BinByBin", None, prior] = refolded

            #covariance
            cov = ROOT.TH2D(unfold.Ereco())
            cov.SetName("{0}_CovMat_{1}_{2}".format(self.fName, "BinByBin", prior))
            cov.GetXaxis().SetTitle("bin number")
            cov.GetYaxis().SetTitle("bin number")
            self.fCovarianceMatrices["BinByBin", None, prior] = cov

            #errors
            unfoldingErr = ROOT.TH1D(unfold.ErecoV(ROOT.RooUnfold.kErrors))
            unfoldingErr.SetName("{0}_UnfoldErr_{1}_Prior{2}".format(self.fName, "BinByBin", prior))
            unfoldingErr.SetTitle("{0} Unfolding Errors {1} Prior={2}".format(self.fName, "BinByBin", prior))
            DMesonJetUtils.DivideNoErrors(unfoldingErr, unfolded)
            diagCov = ROOT.TH1D(unfold.ErecoV(ROOT.RooUnfold.kCovariance))
            diagCov.SetName("{0}_DiagCov_{1}_Prior{2}".format(self.fName, "BinByBin", prior))
            diagCov.SetTitle("{0} Unfolding Errors (covariance) {1} Prior={2}".format(self.fName, "BinByBin", prior))
            DMesonJetUtils.DivideNoErrors(diagCov, unfolded)
            errToy = ROOT.TH1D(unfold.ErecoV(ROOT.RooUnfold.kCovToy))
            errToy.SetName("{0}_MCToyErr_{1}_Prior{2}".format(self.fName, "BinByBin", prior))        
            errToy.SetTitle("{0} Unfolding Errors (MC-toy) {1} Prior={2}".format(self.fName, "BinByBin", prior))
            DMesonJetUtils.DivideNoErrors(errToy, unfolded)
            self.fUnfoldedSpectraErrors["BinByBin", None, prior] = unfoldingErr, diagCov, errToy

    def GenerateResponse(self):
        for prior in self.fPriors:
            print("Generating response for prior {0}".format(prior))
            (detResp, trainTruth), (smallBinDetResp, smallBinTrainTruth) = self.NormalizeAndRebinResponseMatrix(prior)
            detResp.GetZaxis().SetTitle("counts")
            rooUnfoldResp = ROOT.RooUnfoldResponse(ROOT.nullptr, trainTruth, detResp)
            detRespNorm = self.Normalize2D(detResp, self.GenerateFlatPrior())
            detRespNorm.GetZaxis().SetTitle("Probability Density")
            self.fResponseMatrices[prior] = ResponseMatrix("Prior{0}".format(prior), rooUnfoldResp, detResp, trainTruth, smallBinDetResp, smallBinTrainTruth, detRespNorm)

    def NormalizeAndRebinResponseMatrix(self, prior): 
        (normResp, priorHist) = self.NormalizeResponseMatrix(prior)
        coarseResp = self.Rebin(normResp)
        coarsePrior = self.Rebin(priorHist)
        return (coarseResp, coarsePrior), (normResp, priorHist)

    def Rebin(self, hist):
        if isinstance(hist, ROOT.TH2):
            return self.Rebin2D(hist, self.fInputSpectrum, self.fInputSpectrum)
        elif isinstance(hist, ROOT.TH1):
            return self.Rebin1D(hist, self.fInputSpectrum)
        else:
            print("Object {0} of unrecognized type!".format(hist))
            return None # will fail

    def Rebin1D(self, hist, template):
        r = ROOT.TH1D(hist.GetName(), hist.GetTitle(), template.GetXaxis().GetNbins(), template.GetXaxis().GetXbins().GetArray())
        r.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
        r.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
        for xbin in range(0, hist.GetXaxis().GetNbins()+2):
            xbinCenter = hist.GetXaxis().GetBinCenter(xbin)
            rxbin = r.GetXaxis().FindBin(xbinCenter)
            binValue = hist.GetBinContent(xbin) + r.GetBinContent(rxbin)
            binError = math.sqrt(hist.GetBinError(xbin)**2 + r.GetBinError(rxbin)**2)
            r.SetBinContent(rxbin, binValue)
            r.SetBinError(rxbin, binError)
        return r

    def Rebin2D(self, hist, templatex, templatey):
        r = ROOT.TH2D(hist.GetName(), hist.GetTitle(), templatex.GetXaxis().GetNbins(), templatex.GetXaxis().GetXbins().GetArray(), templatey.GetXaxis().GetNbins(), templatey.GetXaxis().GetXbins().GetArray())
        r.GetXaxis().SetTitle(hist.GetXaxis().GetTitle())
        r.GetYaxis().SetTitle(hist.GetYaxis().GetTitle())
        for xbin in range(0, hist.GetXaxis().GetNbins()+2):
            for ybin in range(0, hist.GetYaxis().GetNbins()+2):
                xbinCenter = hist.GetXaxis().GetBinCenter(xbin)
                ybinCenter = hist.GetYaxis().GetBinCenter(ybin)
                rxbin = r.GetXaxis().FindBin(xbinCenter)
                rybin = r.GetYaxis().FindBin(ybinCenter)
                binValue = hist.GetBinContent(xbin, ybin) + r.GetBinContent(rxbin, rybin)
                binError = math.sqrt(hist.GetBinError(xbin, ybin)**2 + r.GetBinError(rxbin, rybin)**2)
                r.SetBinContent(rxbin, rybin, binValue)
                r.SetBinError(rxbin, rybin, binError)
                
        for xbin in range(0, r.GetXaxis().GetNbins()+2):
            for ybin in range(0, r.GetYaxis().GetNbins()+2):
                binValue = r.GetBinContent(xbin, ybin)
                binError = r.GetBinError(xbin, ybin)
                if binValue > 0:
                    relErr = binError / binValue
                    if relErr > 0.9:
                        print("Bin ({0},{1}) has rel stat err = {2}. This is VERY dangerous!".format(xbin,ybin,relErr))
        return r

    def NormalizeResponseMatrix(self, prior):
        if prior == "ResponseTruth":
            return self.fDetectorResponse, self.fDetectorTrainTruth
        elif prior == "Flat":
            priorHist = self.GenerateFlatPrior()
        elif "PowerLaw" in prior:
            priorHist = self.GeneratePowerLawPrior(-int(prior.replace("PowerLaw_","")))
        else:
            print("Prior {0} not implemented!!".format(prior))
            return None # will fail

        priorHist.Scale(self.fDetectorTrainTruth.Integral() / priorHist.Integral())
        priorEffHist = priorHist.Clone("priorEffHist")
        priorEffHist.Multiply(self.fDetectorResponse.ProjectionY())
        priorEffHist.Divide(self.fDetectorTrainTruth)

        resp = self.Normalize2D(self.fDetectorResponse, priorEffHist)

        return resp, priorHist

    def GeneratePowerLawPrior(self, a, peak=3):
        print("Generating power law prior with index {0}".format(a))
        priorHist = self.fDetectorTrainTruth.Clone("myprior")
        priorHist.Reset()
        #this is a modified power law, where at low pT an exponential dumps the function to zero, to avoid infinities
        #other protections are added via the [2] parameter to avoid intermediate infinities in the function evaluation
        #the resulting modified power law has a local maximum, which is set to 3 GeV/c by default
        f = ROOT.TF1("f", "(x>[2])*TMath::Max([2],x)^[0]*TMath::Exp([1]/TMath::Max([2],x))", 0, priorHist.GetXaxis().GetXmax()*2)
        f.SetParameter(0, a)
        f.SetParameter(1, a*peak)
        f.SetParameter(2, 1e-6)
        priorHist.Add(f, 1, "i")
        return priorHist

    def GenerateFlatPrior(self):
        priorHist = self.fDetectorTrainTruth.Clone("myprior")
        for ibin in range(0, priorHist.GetNbinsX()+2):
            priorHist.SetBinContent(ibin, 1)
            priorHist.SetBinError(ibin, 0)
        return priorHist

    def Normalize2D(self, hist, norm):
        resp = hist.Clone("myresp")
        for ybin in range(0, resp.GetYaxis().GetNbins()+2):
            inty = resp.Integral(0, resp.GetXaxis().GetNbins(), ybin, ybin)
            if inty == 0:
                continue
            scaling = norm.GetBinContent(ybin)/inty
            for xbin in range(0, resp.GetXaxis().GetNbins()+2):
                binValue = resp.GetBinContent(xbin, ybin)*scaling
                binErr = resp.GetBinError(xbin, ybin)*scaling
                resp.SetBinContent(xbin, ybin, binValue)
                resp.SetBinError(xbin, ybin, binErr)
        return resp

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

    def StartUnfolding(self, config, eff):
        eng = DMesonJetUnfoldingEngine(config)
        self.fUnfoldingEngine.append(eng)
        r = eng.LoadData(self.fDataFile, self.fResponseFile, eff)
        if r:
            eng.Start()

    def SaveRootFile(self, path):
        if not os.path.isdir(path):
            os.makedirs(path)
        file = ROOT.TFile("{0}/{1}.root".format(path, self.fName), "recreate")
        file.cd()
        for eng in self.fUnfoldingEngine:
            rlist = eng.GenerateRootList()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)
        file.Close()

    def SavePlots(self, path, format):
        if not os.path.isdir(path):
            os.makedirs(path)
        for eng in self.fUnfoldingEngine:
            eng.SavePlots(path, format)
