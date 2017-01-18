#!/usr/bin/env python
# python program to generate a D meson jet response

import math
import os

import ROOT

import DetectorResponse
import DMesonJetProjectors
import DMesonJetUtils
import Axis

globalList = []

class DMesonJetResponseEngine:
    def __init__(self, dmeson, jets, axis, projector):
        self.fDMeson = dmeson
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fAxis = axis
        self.fResponses = None
        self.fCanvases = []

    def SaveRootFile(self, file):
        for resp in self.fResponses.itervalues():
            rlist = resp.GenerateRootList()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))

    def ProjectResponse(self):
        self.fResponses = self.fProjector.GetDetectorResponse(self.fAxis, self.fDMeson, self.fJetDefinitions)

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
                h = s.fHistogram.DrawCopy()
                h.GetYaxis().SetTitle("Probability density")
                h.GetXaxis().SetRangeUser(-1, 0.5)
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
            h.GetYaxis().SetRangeUser(0.04, 0.2)
            h.SetMarkerStyle(ROOT.kFullCircle)
            h.SetMarkerSize(0.8)
            h.SetMarkerColor(ROOT.kBlue + 2)
            h.SetLineColor(ROOT.kBlue + 2)
            globalList.append(h)
            globalList.append(c)
            self.fCanvases.append(c)

    def PlotEnergyScaleShift(self, resp):
        if resp.fEnergyScaleShift:
            c = ROOT.TCanvas("{0}_canvas".format(resp.fEnergyScaleShift.GetName()), resp.fEnergyScaleShift.GetTitle())
            c.cd()
            c.SetLeftMargin(0.18)
            c.SetBottomMargin(0.12)
            c.SetTopMargin(0.08)
            c.SetRightMargin(0.08)
            h = resp.fEnergyScaleShift.DrawCopy("")
            h.GetYaxis().SetRangeUser(-0.1, 0.08)
            h.GetYaxis().SetTitle(resp.fStatistics.fVariableName)
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
            h.SetMarkerStyle(ROOT.kFullCircle)
            h.SetMarkerSize(0.8)
            h.SetMarkerColor(ROOT.kBlue + 2)
            h.SetLineColor(ROOT.kBlue + 2)
            globalList.append(h)

            hMed = resp.fEnergyScaleShiftMedian.DrawCopy("same")
            hMed.SetMarkerStyle(ROOT.kFullSquare)
            hMed.SetMarkerSize(0.8)
            hMed.SetMarkerColor(ROOT.kRed + 2)
            hMed.SetLineColor(ROOT.kRed + 2)

            leg = ROOT.TLegend(0.8, 0.85, 0.9, 0.75)
            leg.SetTextFont(43)
            leg.SetTextSize(16)
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.AddEntry(h, "Mean", "pe")
            leg.AddEntry(hMed, "Median", "pe")
            leg.Draw()

            globalList.append(c)
            globalList.append(leg)
            self.fCanvases.append(c)

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
            h.SetMinimum(DMesonJetUtils.FindMinimum(h) / 2)
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
            c.SetLogz()
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
        self.PlotMultiHistogram(resp, "PartialEfficiencyRatios", "Ratio", "fEfficiency1DRatios")

    def PlotMultiHistogram(self, resp, name, yaxisTitle, listName):
        histList = getattr(resp, listName)
        if not histList:
            return
        cname = "{0}_{1}".format(resp.fName, name)
        c = ROOT.TCanvas(cname, cname)
        c.SetLeftMargin(0.12)
        c.SetBottomMargin(0.12)
        c.SetTopMargin(0.08)
        c.SetRightMargin(0.08)
        c.cd()
        globalList.append(c)
        self.fCanvases.append(c)
        blank = ROOT.TH1D("blankHist", "blankHist;{0};{1}".format(resp.fAxis[1].fTruthAxis.GetTitle(), yaxisTitle), 100, resp.fAxis[1].fTruthAxis.fBins[0], resp.fAxis[1].fTruthAxis.fBins[-1])
        blank.GetXaxis().SetTitleFont(43)
        blank.GetXaxis().SetTitleOffset(1.2)
        blank.GetXaxis().SetTitleSize(19)
        blank.GetXaxis().SetLabelFont(43)
        blank.GetXaxis().SetLabelOffset(0.009)
        blank.GetXaxis().SetLabelSize(18)
        blank.GetYaxis().SetTitleFont(43)
        blank.GetYaxis().SetTitleOffset(1.2)
        blank.GetYaxis().SetTitleSize(19)
        blank.GetYaxis().SetLabelFont(43)
        blank.GetYaxis().SetLabelOffset(0.009)
        blank.GetYaxis().SetLabelSize(18)
        blank.Draw()
        globalList.append(blank)
        colors = [ROOT.kBlue + 2, ROOT.kRed + 2, ROOT.kGreen + 2, ROOT.kOrange + 2, ROOT.kMagenta + 2, ROOT.kAzure + 2, ROOT.kPink + 2, ROOT.kBlack]
        markers = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond, ROOT.kFullStar, ROOT.kFullCross, ROOT.kOpenCircle]
        max = 0;
        leg = ROOT.TLegend(0.15, 0.90, 0.45, 0.65)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(43)
        leg.SetTextSize(16)
        for color, marker, eff in zip(colors[:len(histList) - 1] + colors[-1:], markers[:len(histList) - 1] + markers[-1:], histList):
            h = eff.Clone()
            globalList.append(h)
            h.SetMarkerStyle(marker)
            h.SetMarkerSize(0.9)
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            leg.AddEntry(h, h.GetTitle(), "pe")
            if isinstance(h, ROOT.TGraph):
                h.Draw("P same")
                for i in range(0, h.GetN()):
                    y = h.GetY()[i]
                    if y > max:
                        max = y
            else:
                h.Draw("same")
                for i in range(1, h.GetNbinsX() + 1):
                    y = h.GetBinContent(i)
                    if y > max:
                        max = y

        leg.Draw()
        globalList.append(leg)
        blank.SetMaximum(max * 1.5)

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
            if resp["efficiency"]:
                effWeight = DMesonJetProjectors.EfficiencyWeightCalculator("{0}/{1}".format(self.fProjector.fInputPath, resp["efficiency"]["file_name"]), resp["efficiency"]["list_name"], resp["efficiency"]["object_name"])
            else:
                effWeight = DMesonJetProjectors.SimpleWeight()
            if "apply_to_truth" in resp["efficiency"]:
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

            axis[resp["name"]] = axis_list, effWeight, apply_to_truth, cut_list

        for d_meson in config["d_meson"]:
            eng = DMesonJetResponseEngine(d_meson, config["jets"], axis, self.fProjector)
            self.fResponseEngine.append(eng)
            eng.Start()

    def SaveRootFile(self, path):
        file = ROOT.TFile.Open("{0}/{1}.root".format(path, self.fName), "recreate")
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
