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
    def __init__(self, name, rooUnfoldResponse, response, truth, prior, normResp):
        self.fName = name
        self.fRooUnfoldResponse = rooUnfoldResponse
        self.fResponse = response
        self.fNormResponse = normResp
        self.fTruth = truth
        self.fPrior = prior
        self.GenerateEfficiencies()

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
        rlist.Add(self.fTruth)
        rlist.Add(self.fPrior)
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
        self.fPriors = config["priors"]
        self.UnfoldingConfig = config["methods"]
        self.DefaultPrior = config["default_prior"]
        self.DefaultMethod = config["default_method"]
        self.fName = config["name"]
        self.fResponseMatrices = OrderedDict()
        self.fUnfoldedSpectra = OrderedDict()
        self.fRefoldedSpectra = OrderedDict()
        self.fCovarianceMatrices = OrderedDict()
        self.fSvdDvectors = OrderedDict()
        self.fPearsonMatrices = OrderedDict()
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
        for method in self.UnfoldingConfig.iterkeys():
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
        return rlist

    def LoadData(self, dataFile, responseFile, eff):
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
        if eff:
            self.fDetectorTrainTruth = responseList.FindObject("{0}_{1}_{2}_Truth".format(self.fDMesonResponse, self.fJetDefinition, self.fSpectrumName))
        else:
            self.fDetectorTrainTruth = self.fDetectorResponse.ProjectionY("{0}_{1}_{2}_Truth".format(self.fDMesonResponse, self.fJetDefinition, self.fSpectrumName), 0, -1)
        return True

    def Start(self):
        self.GenerateResponse()
        for method, config in self.UnfoldingConfig.iteritems():
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
        ROOT.gStyle.SetPalette(55,ROOT.nullptr)
        if len(self.UnfoldingConfig) > 1:
            self.CompareMethods()
        if len(self.fPriors) > 1:
            self.ComparePriors()
        self.UnfoldingSummary()
        self.PlotPearsonMatrices()
        self.PlotSvdDvectors()
        self.PlotResponses()

    def PlotResponses(self):
        for name, resp in self.fResponseMatrices.iteritems():
            cname = "Response_{0}_{1}".format(self.fName, name)
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

    def FilterObjectsFromDict(self, objects, method, reg, prior):
        for (methodIt, regIt, priorIt), cov in objects.iteritems():
            if method and not methodIt == method:
                continue
            if prior and not priorIt == prior:
                continue
            if reg and not regIt == reg:
                continue
            yield (methodIt, regIt, priorIt), cov

    def PlotPearsonMatrices(self):
        for method, config in self.UnfoldingConfig.iteritems():
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

        cname = "Pearson_{0}_{1}_{2}".format(self.fName, method, prior)
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
        for method, config in self.UnfoldingConfig.iteritems():
            if isinstance(config["default_reg"], dict):
                reg = config["default_reg"][self.DefaultPrior]
            else:
                reg = config["default_reg"]
            self.UnfoldingSummaryForMethod(method, reg, self.DefaultPrior)
            for prior in self.fPriors:
                self.RegularizationComparisonForMethod(method, prior)

        if self.fTruthSpectrum:
            ratioTruthOverMeas = self.fTruthSpectrum.Clone("{0}_TruthOverMeasured".format(self.fName))
            cnameTOM = "Unfolding_TruthOverMeasured_{0}".format(self.fName)
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
        else:
            if isinstance(self.UnfoldingConfig[method]["default_reg"], dict):
                reg = self.UnfoldingConfig[method]["default_reg"][prior]
            else:
                reg = self.UnfoldingConfig[method]["default_reg"]
            baselineId = (method, reg, prior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
            baseline.SetTitle("Reg = {0}".format(reg))
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
        if len(spectra) < 1:
            return
        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "UnfoldingRegularization_{0}_{1}_{2}".format(method, prior, self.fName))
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def UnfoldingSummaryForMethod(self, method, reg, prior):
        cname = "Unfolding_{0}_{1}_{2}_{3}".format(self.fName, method, reg, prior)
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
        cnameUOM = "Unfolding_UnfoldedOverMeasured_{0}_{1}_{2}_{3}".format(self.fName, method, reg, prior)
        ratioUnfoldedOverMeas.GetYaxis().SetTitle("Unfolded / Measured")
        globalList.append(ratioUnfoldedOverMeas)
        cUOM = ROOT.TCanvas(cnameUOM, cnameUOM)
        self.fCanvases.append(cUOM)
        cUOM.cd()
        globalList.append(cUOM)
        ratioUnfoldedOverMeas.Divide(meas)
        ratioUnfoldedOverMeas.Draw()

        # refolded/measured and unfolded/truth 
        cnameROM = "Unfolding_RefoldedOverMeasured_{0}_{1}_{2}_{3}".format(self.fName, method, reg, prior)
        cROM = ROOT.TCanvas(cnameROM, cnameROM)
        self.fCanvases.append(cROM)
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
            if isinstance(self.UnfoldingConfig[method]["default_reg"], dict):
                reg = self.UnfoldingConfig[method]["default_reg"][self.DefaultPrior]
            else:
                reg = self.UnfoldingConfig[method]["default_reg"]
            baselineId = (self.DefaultMethod, reg, self.DefaultPrior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
        globalList.append(baseline)
        spectra = []
        for method, config in self.UnfoldingConfig.iteritems():
            if isinstance(config["default_reg"], dict):
                reg = config["default_reg"][self.DefaultPrior]
            else:
                reg = config["default_reg"]
            id = (method, reg, self.DefaultPrior)
            if id == baselineId:
                continue
            h = self.fUnfoldedSpectra[id].Clone("{0}_copy".format(self.fUnfoldedSpectra[id].GetName()))
            spectra.append(h)
            globalList.append(h)

        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "UnfoldingMethod_{0}".format(self.fName))
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

    def ComparePriors(self):
        for method in self.UnfoldingConfig.iterkeys():
            self.ComparePriorsForMethod(method)

    def ComparePriorsForMethod(self, method):
        config = self.UnfoldingConfig[method]
        if self.fTruthSpectrum:
            baselineId = None
            baseline = self.fTruthSpectrum.Clone("{0}_copy".format(self.fTruthSpectrum.GetName()))
            baseline.SetTitle("MC Truth")
        else:
            if isinstance(config["default_reg"], dict):
                reg = config["default_reg"][self.DefaultPrior]
            else:
                reg = config["default_reg"]
            baselineId = (self.DefaultMethod, reg, self.DefaultPrior)
            baseline = self.fUnfoldedSpectra[baselineId].Clone("{0}_copy".format(self.fUnfoldedSpectra[baselineId].GetName()))
        globalList.append(baseline)
        spectra = []
        for prior in self.fPriors:
            if isinstance(config["default_reg"], dict):
                reg = config["default_reg"][prior]
            else:
                reg = config["default_reg"]
            id = (method, reg, prior)
            if id == baselineId:
                continue
            h = self.fUnfoldedSpectra[id].Clone("{0}_copy".format(self.fUnfoldedSpectra[id].GetName()))
            spectra.append(h)
            globalList.append(h)

        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "UnfoldingPrior_{0}_{1}".format(self.fName, method))
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
                    dvector = unfold.Impl().GetD().Clone("SvdDvector_{0}_{1}".format(self.fName, prior))
                    dvector.GetXaxis().SetTitle("k")
                    dvector.GetYaxis().SetTitle("d")
                    self.fSvdDvectors[prior] = dvector
                cov = ROOT.TH2D(unfold.Ereco())
                cov.SetName("CovMat_{0}_{1}_{2}_{3}".format(self.fName, "Svd", reg, prior))
                cov.GetXaxis().SetTitle("bin number")
                cov.GetYaxis().SetTitle("bin number")
                self.fCovarianceMatrices["Svd", reg, prior] = cov
                unfolded.SetName("{0}_{1}_{2}_{3}".format(self.fName, "Svd", reg, prior))
                unfolded.SetTitle("{0}, k={1}, prior={2}".format("Svd", reg, prior))
                self.fUnfoldedSpectra["Svd", reg, prior] = unfolded
                refolded = resp.fRooUnfoldResponse.ApplyToTruth(unfolded)
                refolded.SetName("{0}_{1}_{2}_{3}_Refolded".format(self.fName, "Svd", reg, prior))
                refolded.SetTitle("{0}, k={1}, prior={2}".format("Svd", reg, prior))
                self.fRefoldedSpectra["Svd", reg, prior] = refolded
    
    def UnfoldBayes(self, iter_min, iter_max, iter_step):
        for prior, resp in self.fResponseMatrices.iteritems():
            for reg in range(iter_min, iter_max+1, iter_step):
                print("Unfolding {0}, reg={1}, prior={2}".format("Bayes", reg, prior))
                unfold = ROOT.RooUnfoldBayes(resp.fRooUnfoldResponse, self.fInputSpectrum, reg)
                unfolded = unfold.Hreco()
                cov = ROOT.TH2D(unfold.Ereco())
                cov.SetName("CovMat_{0}_{1}_{2}_{3}".format(self.fName, "Bayes", reg, prior))
                cov.GetXaxis().SetTitle("bin number")
                cov.GetYaxis().SetTitle("bin number")
                self.fCovarianceMatrices["Bayes", reg, prior] = cov
                unfolded.SetName("{0}_{1}_{2}_{3}".format(self.fName, "Bayes", reg, prior))
                unfolded.SetTitle("{0}, iter={1}, prior={2}".format("Bayes", reg, prior))
                self.fUnfoldedSpectra["Bayes", reg, prior] = unfolded
                refolded = resp.fRooUnfoldResponse.ApplyToTruth(unfolded)
                refolded.SetName("{0}_{1}_{2}_{3}_Refolded".format(self.fName, "Bayes", reg, prior))
                refolded.SetTitle("{0}, iter={1}, prior={2}".format("Bayes", reg, prior))
                self.fRefoldedSpectra["Bayes", reg, prior] = refolded

    def UnfoldBinByBin(self):
        for prior, resp in self.fResponseMatrices.iteritems():
            print("Unfolding {0}, prior={1}".format("BinByBin", prior))
            unfold = ROOT.RooUnfoldBinByBin(resp.fRooUnfoldResponse, self.fInputSpectrum)
            unfolded = unfold.Hreco()
            cov = ROOT.TH2D(unfold.Ereco())
            cov.SetName("CovMat_{0}_{1}_{2}".format(self.fName, "BinByBin", prior))
            cov.GetXaxis().SetTitle("bin number")
            cov.GetYaxis().SetTitle("bin number")
            self.fCovarianceMatrices["BinByBin", None, prior] = cov
            unfolded.SetName("{0}_{1}_{2}".format(self.fName, "BinByBin", prior))
            unfolded.SetTitle("{0}, prior={1}".format("BinByBin", prior))
            self.fUnfoldedSpectra["BinByBin", None, prior] = unfolded
            refolded = resp.fRooUnfoldResponse.ApplyToTruth(unfolded)
            refolded.SetName("{0}_{1}_{2}_Refolded".format(self.fName, "BinByBin", prior))
            refolded.SetTitle("{0}, prior={1}".format("BinByBin", prior))
            self.fRefoldedSpectra["BinByBin", None, prior] = refolded

    def GenerateResponse(self):
        for prior in self.fPriors:
            (detResp, trainTruth) = self.NormalizeAndRebinResponseMatrix(prior)
            detResp.GetZaxis().SetTitle("yield")
            rooUnfoldResp = ROOT.RooUnfoldResponse(ROOT.nullptr, trainTruth, detResp)
            priorSpectrum = trainTruth.Clone("Prior_{0}".format(prior))
            priorSpectrum.Scale(1. / priorSpectrum.Integral())
            detRespNorm = self.Normalize2D(detResp, self.GenerateFlatPrior())
            detRespNorm.GetZaxis().SetTitle("Probability Density")
            self.fResponseMatrices[prior] = ResponseMatrix(prior, rooUnfoldResp, detResp, trainTruth, priorSpectrum, detRespNorm)

    def NormalizeAndRebinResponseMatrix(self, prior): 
        (normResp, priorHist) = self.NormalizeResponseMatrix(prior)
        coarseResp = self.Rebin(normResp)
        coarsePrior = self.Rebin(priorHist)
        return coarseResp, coarsePrior

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
        return r

    def NormalizeResponseMatrix(self, prior):
        if prior == "train_truth":
            return self.fDetectorResponse, self.fDetectorTrainTruth
        elif prior == "flat":
            priorHist = self.GenerateFlatPrior()
        else:
            print("Prior {0} not implemented!!".format(prior))
            return None # will fail

        priorEffHist = priorHist.Clone()
        priorEffHist.Multiply(self.fDetectorResponse.ProjectionY())
        priorEffHist.Divide(self.fDetectorTrainTruth)

        resp = self.Normalize2D(self.fDetectorResponse, priorEffHist)

        return resp, priorHist

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
                resp.SetBinContent(xbin, ybin, resp.GetBinContent(xbin, ybin)*scaling)
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