#!/usr/bin/env python
#python program to generate a D meson jet response

import ROOT
import math
import DMesonJetProjectors
from DMesonJetBase import *
from array import array
from copy import deepcopy

globalList = []

class DMesonJetResponseEngine:
    def __init__(self, dmeson, jets, axis, projector):
        self.fDMeson = dmeson
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fAxis = axis
        self.fResponses = None
        
    def SaveRootFile(self, file):
        for resp in self.fResponses.itervalues():
            rlist = resp.GenerateRootList()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)
    
    def ProjectResponse(self):
        self.fResponses = self.fProjector.GetDetectorResponse(self.fAxis, self.fDMeson, self.fJetDefinitions)
        
    def Start(self):
        self.ProjectResponse()
        self.GenerateEfficiency()
        self.PlotResponse()

    def GenerateEfficiency(self):
        for resp in self.fResponses.itervalues():
            resp.GenerateFoldedTruth()
            resp.GenerateEfficiency()

    def PlotResponse(self):
        for resp in self.fResponses.itervalues():
            self.PlotResponseMatrix(resp)
            self.PlotEfficiency(resp)

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
            h.SetMarkerColor(ROOT.kBlue+2)
            h.SetLineColor(ROOT.kBlue+2)
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
        
    def GenerateMultiCanvas(self, name, n):
        rows = int(math.floor(math.sqrt(n)))
        cols = int(math.ceil(float(n) / rows))
        cname = "{0}_{1}".format(self.fDMeson, name)
        c = ROOT.TCanvas(cname, cname, cols*400, rows*400)
        c.Divide(cols, rows)
        globalList.append(c)
        return c
    
    def PlotPartialMultiEfficiency(self, resp):
        cname = "{0}_{1}_PartialEfficiency".format(self.fDMeson, resp.fName)
        c = ROOT.TCanvas(cname, cname)
        c.SetLeftMargin(0.12)
        c.SetBottomMargin(0.12)
        c.SetTopMargin(0.08)
        c.SetRightMargin(0.08)
        c.cd()
        globalList.append(c)
        blank = ROOT.TH1D("blankHist", "blankHist;{0};Efficiency".format(resp.fAxis[1].fTruthAxis.GetTitle()), 100, resp.fAxis[1].fTruthAxis.fBins[0], resp.fAxis[1].fTruthAxis.fBins[-1])
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
        colors = [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kMagenta+2]
        markers = [ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond]
        max = 0;
        leg = ROOT.TLegend(0.15, 0.90, 0.45, 0.65)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(43)
        leg.SetTextSize(16)
        for color,marker,eff in zip(colors,markers,resp.fEfficiency1D):
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
                for i in range(1, h.GetNbinsX()+1):
                    y = h.GetBinContent(i)
                    if y > max:
                        max = y

        leg.Draw()
        globalList.append(leg)
        blank.SetMaximum(max*1.5)

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
            axis_list = []
            axis_names = ["jet_pt", "d_pt", "d_z"]
            for axis_name in axis_names:
                if resp[axis_name]:
                    axis_list.append(ResponseAxis(Axis(axis_name, resp[axis_name]["reco"], "reco"), Axis(axis_name, resp[axis_name]["truth"], "truth")))

            axis[resp["name"]] = axis_list, effWeight

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
