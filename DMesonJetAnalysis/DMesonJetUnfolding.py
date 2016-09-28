#!/usr/bin/env python
#python program to perform a D meson jet unfolding

import ROOT
import math
import DMesonJetProjectors
from DMesonJetBase import *
import array
import copy
import DMesonJetUtils
import collections

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
        print(self.UnfoldingConfig)
        self.DefaultPrior = config["default_prior"]
        self.DefaultMethod = config["default_method"]
        self.fResponseMatrices = dict()
        self.fUnfoldedSpectra = dict()

    def Start(self):
        self.GenerateRooUnfoldResponse()
        for method, config in self.UnfoldingConfig.iteritems():
            if method == "SVD":
                self.UnfoldSVD()
            elif method == "Bayes":
                self.UnfoldBayes(config["iter_min"], config["iter_max"], config["iter_step"])
            else:
                print("Unfolding method {0} not known!".format(method))
        self.Plot()

    def Plot(self):
        if len(self.UnfoldingConfig) > 1:
            self.CompareMethods()

    def CompareMethods(self):
        if self.fTruthSpectrum:
            baselineId = None
            baseline = self.fTruthSpectrum
        else:
            baselineId = (self.DefaultMethod, self.UnfoldingConfig[self.DefaultMethod]["default_reg"], self.DefaultPrior)
            baseline = self.fUnfoldedSpectra[baselineId]

        spectra = []
        for method, config in self.UnfoldingConfig.iteritems():
            id = (method, config["default_reg"], self.DefaultPrior)
            if id == baselineId:
                continue
            spectra.append(self.fUnfoldedSpectra[id])

        r = DMesonJetUtils.CompareSpectra(baseline, spectra, "UnfoldingMethod")
        for obj in r:
            globalList.append(obj)

    def UnfoldSVD(self):
        max_reg = self.fInputSpectrum.GetNbinsX()
        for prior, resp in self.fResponseMatrices.iteritems():
            for reg in range(1, max_reg+1):
                unfold = ROOT.RooUnfoldSvd(resp, self.fInputSpectrum, reg)
                unfolded = unfold.Hreco()
                unfolded.SetName("{0}_{1}_{2}_{3}".format(self.fInputSpectrum.GetName(), "SVD", reg, prior))
                unfolded.SetTitle("{0}, k={1}, prior={2}".format("SVD", reg, prior))
                self.fUnfoldedSpectra["SVD", reg, prior] = unfolded
    
    def UnfoldBayes(self, iter_min, iter_max, iter_step):
        for prior, resp in self.fResponseMatrices.iteritems():
            for reg in range(iter_min, iter_max, iter_step):
                unfold = ROOT.RooUnfoldBayes(resp, self.fInputSpectrum, reg)
                unfolded = unfold.Hreco()
                unfolded.SetName("{0}_{1}_{2}_{3}".format(self.fInputSpectrum.GetName(), "Bayes", reg, prior))
                unfolded.SetTitle("{0}, iter={1}, prior={2}".format("Bayes", reg, prior))
                self.fUnfoldedSpectra["Bayes", reg, prior] = unfolded

    def GenerateRooUnfoldResponse(self):
        for prior in self.fPriors:
            if prior == "train_truth":
                rooUnfoldResp = ROOT.RooUnfoldResponse(ROOT.nullptr, self.fDetectorTrainTruth, self.fDetectorResponse)
            else:
                (detResp, trainTruth) = self.NormalizeResponseMatrix(self.fDetectorResponse, prior, self.fDetectorEfficiency)
                rooUnfoldResp = ROOT.RooUnfoldResponse(ROOT.nullptr, trainTruth, detResp)
            self.fResponseMatrices[prior] = rooUnfoldResp

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
