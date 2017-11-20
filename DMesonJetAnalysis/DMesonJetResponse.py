#!/usr/bin/env python
# python program to generate a D meson jet response

import math
import os

import ROOT

import DetectorResponse
import DMesonJetProjectors
import DMesonJetCompare
import DMesonJetUtils
import Axis

globalList = []

class DMesonJetResponseContainer:
    def __init__(self, trigger, DMesonDef, jetDefinitions, respDefinitions):
        self.fResponseObjects = dict()
        for jetDef in jetDefinitions:
            jetName = "Jet_AKT{0}{1}_pt_scheme".format(jetDef["type"], jetDef["radius"])
            for axisName, (axisDef, efficiency, applyEffTruth, cuts) in respDefinitions.iteritems():
                respName = "_".join([obj for obj in [trigger, DMesonDef, jetName, axisName] if obj])
                if efficiency:
                    eff_file_name = efficiency["file_name"]
                    eff_list_name = "_".join([obj for obj in [trigger, DMesonDef, jetName, efficiency["list_name"]] if obj])
                    eff_obj_name = "_".join([obj for obj in [trigger, DMesonDef, jetName, efficiency["list_name"], efficiency["object_name"]] if obj])
                    effWeight = DMesonJetProjectors.EfficiencyWeightCalculator(eff_file_name, eff_list_name, eff_obj_name)
                else:
                    effWeight = DMesonJetProjectors.SimpleWeight()
                resp = DetectorResponse.DetectorResponse(respName, jetName, axisDef, cuts, effWeight, applyEffTruth)
                resp.GenerateHistograms()
                self.fResponseObjects[respName] = resp

    def Fill(self, event, eventWeight):
        for r in self.fResponseObjects.itervalues():
            r.Fill(event, eventWeight)

class DMesonJetResponseEngine:
    def __init__(self, trigger, dmeson, d_meson_cuts, jets, axis, projector):
        self.fDMeson = dmeson
        self.fDMesonCuts = d_meson_cuts
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fAxis = axis
        self.fResponses = None
        self.fCanvases = []
        self.fTrigger = trigger

    def SaveRootFile(self, file):
        for resp in self.fResponses.itervalues():
            rlist = resp.GenerateRootList()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))

    def ProjectResponse(self):
        response_container = DMesonJetResponseContainer(self.fTrigger, "{}_{}".format(self.fDMeson, self.fDMesonCuts), self.fJetDefinitions, self.fAxis)
        self.fProjector.StartProjection(self.fTrigger, "{}_kSignalOnly".format(self.fDMeson), self.fDMesonCuts, response_container)
        self.fResponses = response_container.fResponseObjects

    def Start(self):
        self.ProjectResponse()
        self.GenerateEfficiency()
        self.GenerateResolution()
        self.GenerateResponseUncertainty()
        self.PlotResponse()

    def GenerateEfficiency(self):
        for resp in self.fResponses.itervalues():
            resp.GenerateEfficiency()

    def GenerateResolution(self):
        for resp in self.fResponses.itervalues():
            resp.GenerateResolution()

    def GenerateResponseUncertainty(self):
        for resp in self.fResponses.itervalues():
            resp.GenerateResponseUncertainty()

    def PlotResponse(self):
        for resp in self.fResponses.itervalues():
            self.PlotResponseMatrix(resp)
            self.PlotResponseMatrixUncertainty(resp)
            self.PlotEfficiency(resp)
            self.PlotResolution(resp)
            self.PlotEnergyScaleShift(resp)
            self.PlotDetectorResponseSummary(resp)

    def PlotDetectorResponseSummary(self, resp):
        if resp.fStatistics:
            c = DMesonJetUtils.GenerateMultiCanvas(resp.fStatistics.fName, len(resp.fStatistics.fStatisticSets))
            globalList.append(c)
            self.fCanvases.append(c)
            for i, s in enumerate(resp.fStatistics.fStatisticSets):
                if not s.fHistogram:
                    continue
                pad = c.cd(i + 1)
                # pad.SetLogy()
                pad.SetLeftMargin(0.14)
                pad.SetRightMargin(0.05)
                pad.SetTopMargin(0.08)
                pad.SetBottomMargin(0.18)
                pad.SetGridx()
                pad.SetGridy()
                h = s.fHistogram.DrawCopy()
                h.GetYaxis().SetTitle("Probability density")
                if resp.fStatistics.fAxis.fName == "jet_pt":
                    h.GetXaxis().SetRangeUser(-1, 0.6)
                elif resp.fStatistics.fAxis.fName == "d_z":
                    h.GetXaxis().SetRangeUser(-0.6, 1)
                h.SetMarkerColor(ROOT.kRed + 2)
                h.SetMarkerStyle(ROOT.kFullCircle)
                h.SetMarkerSize(0.8)
                h.SetLineColor(ROOT.kRed + 2)
                globalList.append(h)
                h.GetXaxis().SetTitleFont(43)
                h.GetXaxis().SetTitleOffset(2.8)
                h.GetXaxis().SetTitleSize(19)
                h.GetXaxis().SetLabelFont(43)
                h.GetXaxis().SetLabelOffset(0.009)
                h.GetXaxis().SetLabelSize(18)
                h.GetYaxis().SetTitleFont(43)
                h.GetYaxis().SetTitleOffset(2.6)
                h.GetYaxis().SetTitleSize(19)
                h.GetYaxis().SetLabelFont(43)
                h.GetYaxis().SetLabelOffset(0.009)
                h.GetYaxis().SetLabelSize(18)
                # h.SetMaximum(h.GetMaximum()*1.3)
                h.SetMinimum(0)
                htitle = ROOT.TPaveText(0.12, 0.99, 0.95, 0.93, "NB NDC")
                htitle.SetBorderSize(0)
                htitle.SetFillStyle(0)
                htitle.SetTextFont(43)
                htitle.SetTextSize(18)
                htitle.AddText(h.GetTitle())
                htitle.Draw()
                globalList.append(htitle)

    def PlotResolution(self, resp):
        if resp.fResolution:
            c = ROOT.TCanvas("{0}_canvas".format(resp.fResolution.GetName()), resp.fResolution.GetTitle())
            c.cd()
            c.SetLeftMargin(0.18)
            c.SetBottomMargin(0.12)
            c.SetTopMargin(0.08)
            c.SetRightMargin(0.08)
            c.SetGridx()
            c.SetGridy()
            h = resp.fResolution.DrawCopy("")
            h.GetXaxis().SetTitleFont(43)
            h.GetXaxis().SetTitleOffset(1.2)
            h.GetXaxis().SetTitleSize(19)
            h.GetXaxis().SetLabelFont(43)
            h.GetXaxis().SetLabelOffset(0.009)
            h.GetXaxis().SetLabelSize(18)
            h.GetYaxis().SetTitleFont(43)
            h.GetYaxis().SetTitleOffset(1.8)
            h.GetYaxis().SetTitleSize(19)
            h.GetYaxis().SetLabelFont(43)
            h.GetYaxis().SetLabelOffset(0.009)
            h.GetYaxis().SetLabelSize(18)
            # h.GetYaxis().SetRangeUser(0.04, 0.2)
            h.SetMarkerStyle(ROOT.kFullCircle)
            h.SetMarkerSize(0.8)
            h.SetMarkerColor(ROOT.kBlue + 2)
            h.SetLineColor(ROOT.kBlue + 2)
            globalList.append(h)
            globalList.append(c)
            self.fCanvases.append(c)

    def PlotEnergyScaleShift(self, resp):
        if resp.fEnergyScaleShift:
            cname = "{0}_canvas".format(resp.fEnergyScaleShift.GetName())
            comp = DMesonJetCompare.DMesonJetCompare(cname)
            comp.fDoSpectraPlot = "lineary"
            comp.fDoRatioPlot = False
            comp.fColors = [ROOT.kBlue + 2, ROOT.kRed + 2]
            comp.fMarkers = [ROOT.kFullCircle, ROOT.kOpenSquare]
            comp.fMinimumLimit = float("-inf")
            comp.fOptSpectrum = "e0"

            hMean = resp.fEnergyScaleShift.Clone("{}_copy".format(resp.fEnergyScaleShift.GetName()))
            globalList.append(hMean)
            hMean.SetTitle("Mean")

            hMed = resp.fEnergyScaleShiftMedian.Clone("{}_copy".format(resp.fEnergyScaleShiftMedian.GetName()))
            globalList.append(hMed)
            hMed.SetTitle("Median")

            results = comp.CompareSpectra(hMean, [hMed])

            comp.fCanvasSpectra.SetGridx()
            comp.fCanvasSpectra.SetGridy()

            for obj in results:
                globalList.append(obj)

            self.fCanvases.append(comp.fCanvasSpectra)

    def PlotResponseMatrix(self, resp):
        if len(resp.fAxis) == 1:
            c = ROOT.TCanvas("{0}_canvas".format(resp.fResponseMatrix.GetName()), resp.fResponseMatrix.GetTitle())
            c.cd()
            c.SetRightMargin(0.17)
            c.SetTopMargin(0.08)
            c.SetBottomMargin(0.14)
            c.SetLeftMargin(0.12)
            c.SetLogz()
            h = resp.fResponseMatrix.DrawCopy("colz")
            find_minimum = DMesonJetUtils.FindMinimum(h)
            if find_minimum is not None: h.SetMinimum(DMesonJetUtils.FindMinimum(h) / 2)
            h.GetXaxis().SetTitleFont(43)
            h.GetXaxis().SetTitleOffset(1.3)
            h.GetXaxis().SetTitleSize(21)
            h.GetXaxis().SetLabelFont(43)
            h.GetXaxis().SetLabelOffset(0.009)
            h.GetXaxis().SetLabelSize(19)
            h.GetYaxis().SetTitleFont(43)
            h.GetYaxis().SetTitleOffset(1.2)
            h.GetYaxis().SetTitleSize(21)
            h.GetYaxis().SetLabelFont(43)
            h.GetYaxis().SetLabelOffset(0.009)
            h.GetYaxis().SetLabelSize(19)
            h.GetZaxis().SetTitleFont(43)
            h.GetZaxis().SetTitleOffset(1.2)
            h.GetZaxis().SetTitleSize(21)
            h.GetZaxis().SetLabelFont(43)
            h.GetZaxis().SetLabelOffset(0.009)
            h.GetZaxis().SetLabelSize(19)
            globalList.append(h)
            globalList.append(c)
            self.fCanvases.append(c)

    def PlotResponseMatrixUncertainty(self, resp):
        if resp.fResponseMatrixUncertainty:
            c = ROOT.TCanvas("{0}_canvas".format(resp.fResponseMatrixUncertainty.GetName()), resp.fResponseMatrixUncertainty.GetTitle())
            c.cd()
            c.SetRightMargin(0.17)
            c.SetTopMargin(0.08)
            c.SetBottomMargin(0.14)
            c.SetLeftMargin(0.12)
            h = resp.fResponseMatrixUncertainty.DrawCopy("colz")
            h.GetXaxis().SetTitleFont(43)
            h.GetXaxis().SetTitleOffset(1.3)
            h.GetXaxis().SetTitleSize(21)
            h.GetXaxis().SetLabelFont(43)
            h.GetXaxis().SetLabelOffset(0.009)
            h.GetXaxis().SetLabelSize(19)
            h.GetYaxis().SetTitleFont(43)
            h.GetYaxis().SetTitleOffset(1.2)
            h.GetYaxis().SetTitleSize(21)
            h.GetYaxis().SetLabelFont(43)
            h.GetYaxis().SetLabelOffset(0.009)
            h.GetYaxis().SetLabelSize(19)
            h.GetZaxis().SetTitleFont(43)
            h.GetZaxis().SetTitleOffset(1.2)
            h.GetZaxis().SetTitleSize(21)
            h.GetZaxis().SetLabelFont(43)
            h.GetZaxis().SetLabelOffset(0.009)
            h.GetZaxis().SetLabelSize(19)
            globalList.append(h)
            globalList.append(c)
            self.fCanvases.append(c)

    def PlotEfficiency(self, resp):
        if len(resp.fAxis) == 1:
            c = ROOT.TCanvas("{0}_canvas".format(resp.fEfficiency.GetName()), resp.fEfficiency.GetTitle())
            c.cd()
            c.SetLeftMargin(0.12)
            c.SetBottomMargin(0.12)
            c.SetTopMargin(0.08)
            c.SetRightMargin(0.08)
            c.SetGridx()
            c.SetGridy()
            h = resp.fEfficiency.Clone()
            h.SetMarkerStyle(ROOT.kFullCircle)
            h.SetMarkerSize(0.9)
            h.SetMarkerColor(ROOT.kBlue + 2)
            h.SetLineColor(ROOT.kBlue + 2)
            h.GetXaxis().SetTitleFont(43)
            h.GetXaxis().SetTitleOffset(1.2)
            h.GetXaxis().SetTitleSize(19)
            h.GetXaxis().SetLabelFont(43)
            h.GetXaxis().SetLabelOffset(0.009)
            h.GetXaxis().SetLabelSize(18)
            h.GetYaxis().SetTitleFont(43)
            h.GetYaxis().SetTitleOffset(1.2)
            h.GetYaxis().SetTitleSize(19)
            h.GetYaxis().SetLabelFont(43)
            h.GetYaxis().SetLabelOffset(0.009)
            h.GetYaxis().SetLabelSize(18)
            if isinstance(h, ROOT.TGraph):
                h.Draw("AP")
            else:
                h.Draw()
            c.Update()
            if c.GetUymax() > 1.5: h.GetYaxis().SetRangeUser(c.GetUymin(), 1.5)
        elif len(resp.fAxis) == 2:
            c = ROOT.TCanvas("{0}_canvas".format(resp.fEfficiency.GetName()), resp.fEfficiency.GetTitle())
            c.SetRightMargin(0.17)
            c.SetTopMargin(0.08)
            c.SetBottomMargin(0.14)
            c.SetLeftMargin(0.12)
            c.cd()
            h = resp.fEfficiency.DrawCopy("colz")
            h.GetXaxis().SetTitleFont(43)
            h.GetXaxis().SetTitleOffset(1.3)
            h.GetXaxis().SetTitleSize(21)
            h.GetXaxis().SetLabelFont(43)
            h.GetXaxis().SetLabelOffset(0.009)
            h.GetXaxis().SetLabelSize(19)
            h.GetYaxis().SetTitleFont(43)
            h.GetYaxis().SetTitleOffset(1.2)
            h.GetYaxis().SetTitleSize(21)
            h.GetYaxis().SetLabelFont(43)
            h.GetYaxis().SetLabelOffset(0.009)
            h.GetYaxis().SetLabelSize(19)
            h.GetZaxis().SetTitleFont(43)
            h.GetZaxis().SetTitleOffset(1.2)
            h.GetZaxis().SetTitleSize(21)
            h.GetZaxis().SetLabelFont(43)
            h.GetZaxis().SetLabelOffset(0.009)
            h.GetZaxis().SetLabelSize(19)
            c.Update()
            if h.GetZaxis().GetXmax() > 1.5: h.SetMaximum(1.5)

            self.PlotPartialMultiEfficiency(resp)

        globalList.append(h)
        globalList.append(c)
        self.fCanvases.append(c)

    def GenerateMultiCanvas(self, name, n):
        rows = int(math.floor(math.sqrt(n)))
        cols = int(math.ceil(float(n) / rows))
        cname = "{0}_{1}".format(self.fDMeson, name)
        c = ROOT.TCanvas(cname, cname, cols * 400, rows * 400)
        c.Divide(cols, rows)
        globalList.append(c)
        self.fCanvases.append(c)
        return c

    def PlotPartialMultiEfficiency(self, resp):
        self.PlotMultiHistogram(resp, "PartialEfficiency", "Efficiency", "fEfficiency1D")

    def PlotMultiHistogram(self, resp, name, yaxisTitle, listName):
        histList = getattr(resp, listName)
        if not histList: return

        cname = "{0}_{1}".format(resp.fName, name)
        comp = DMesonJetCompare.DMesonJetCompare(cname)
        comp.fDoSpectraPlot = "lineary"
        comp.fDoRatioPlot = "lineary"

        colors = [ROOT.kBlue + 2, ROOT.kRed + 2, ROOT.kGreen + 2, ROOT.kOrange + 2, ROOT.kMagenta + 2, ROOT.kAzure + 2, ROOT.kPink + 2, ROOT.kBlack]
        markers = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond, ROOT.kFullStar, ROOT.kFullCross, ROOT.kOpenCircle]
        comp.fColors = colors[:len(histList) - 1] + colors[-1:]
        comp.fMarkers = markers[:len(histList) - 1] + markers[-1:]

        results = comp.CompareSpectra(histList[0], histList[1:])

        comp.fCanvasSpectra.SetGridx()
        comp.fCanvasSpectra.SetGridy()
        comp.fCanvasRatio.SetGridx()
        comp.fCanvasRatio.SetGridy()

        comp.fCanvasSpectra.Update()
        comp.fCanvasRatio.Update()

        if comp.fCanvasSpectra.GetUymax() > 1.5: comp.fMainHistogram.SetMaximum(1.5)
        if comp.fCanvasRatio.GetUymax() > 1.5: comp.fMainRatioHistogram.SetMaximum(1.5)

        for obj in results:
            globalList.append(obj)

        self.fCanvases.append(comp.fCanvasSpectra)
        self.fCanvases.append(comp.fCanvasRatio)

class DMesonJetResponse:
    def __init__(self, name):
        self.fName = name
        self.fResponseEngine = []

    def SetProjector(self, projector):
        self.fProjector = projector

    def StartAnalysis(self, config):
        axis = dict()
        for resp in config["response"]:
            if not resp["active"]:
                continue
            if "efficiency" in resp and resp["efficiency"]:
                resp["efficiency"]["file_name"] = "{}/{}".format(self.fProjector.fInputPath, resp["efficiency"]["file_name"])
            else:
                resp["efficiency"] = False
            if resp["efficiency"] and "apply_to_truth" in resp["efficiency"]:
                apply_to_truth = resp["efficiency"]["apply_to_truth"]
            else:
                apply_to_truth = False
            axis_list = []
            for axis_name, bins in resp["bins"].iteritems():
                if "bins" in bins:
                    a = DetectorResponse.ResponseAxis(Axis.Axis(axis_name, bins["bins"], "reco"), Axis.Axis(axis_name, bins["bins"], "truth"))
                else:
                    a = DetectorResponse.ResponseAxis(Axis.Axis.fromLimits(axis_name, bins["min"], bins["max"], bins["width"], "reco"), Axis.Axis.fromLimits(axis_name, bins["min"], bins["max"], bins["width"], "truth"))
                if "coarse" in bins:
                    a.SetCoarseAxis(Axis.Axis(axis_name, bins["coarse"], "reco"), Axis.Axis(axis_name, bins["coarse"], "truth"))
                if axis_name == "jet_pt":
                    axis_list.insert(0, a)
                else:
                    axis_list.append(a)

            if "cuts" in resp:
                cut_list = resp["cuts"]
            else:
                cut_list = []

            axis[resp["name"]] = axis_list, resp["efficiency"], apply_to_truth, cut_list

        for trigger in config["trigger"]:
            for d_meson in config["d_meson"]:
                for d_meson_cuts in config["d_meson_cuts"]:
                    eng = DMesonJetResponseEngine(trigger, d_meson, d_meson_cuts, config["jets"], axis, self.fProjector)
                    self.fResponseEngine.append(eng)
                    eng.Start()

    def SaveRootFile(self, path):
        fname = "{0}/{1}.root".format(path, self.fName)
        file = ROOT.TFile.Open(fname, "recreate")
        if not file or file.IsZombie():
            print("Could not create file '{}'!".format(fname))
            return
        file.cd()
        for eng in self.fResponseEngine:
            eng.SaveRootFile(file)
        file.Close()

    def SavePlots(self, path, format):
        fullPath = "{0}/{1}/{2}".format(path, self.fName, format)
        if not os.path.isdir(fullPath):
            os.makedirs(fullPath)
        for eng in self.fResponseEngine:
            eng.SavePlots(fullPath, format)
