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
        
    def ProjectResponse(self):
        self.fResponses = self.fProjector.GetDetectorResponse(self.fAxis, self.fDMeson, self.fJetDefinitions)
        
    def Start(self):
        self.ProjectResponse()
        self.GenerateEfficiency()
        self.PlotResponse()

    def GenerateEfficiency(self):
        for resp in self.fResponses.itervalues():
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
            h.Draw("AP")
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
            
        globalList.append(h)
        globalList.append(c)

class DMesonJetResponse:
    def __init__(self, name):
        self.fName = name
        self.fResponseEngine = []

    def SetProjector(self, projector):
        self.fProjector = projector
        
    def StartAnalysis(self, config):
        axis = dict()
        for resp in config["response"]:
            axis[resp["name"]] = []
            axis_names = ["jet_pt", "d_pt", "d_z"]
            for axis_name in axis_names:
                if resp[axis_name]:
                    axis[resp["name"]].append(ResponseAxis(Axis(axis_name, resp[axis_name]["reco"], "reco"), Axis(axis_name, resp[axis_name]["truth"], "truth")))
        
        for d_meson in config["d_meson"]:
            eng = DMesonJetResponseEngine(d_meson, config["jets"], axis, self.fProjector)
            self.fResponseEngine.append(eng)
            eng.Start()
