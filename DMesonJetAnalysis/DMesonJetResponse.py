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
        self.fResponses = dict()
        
    def ProjectResponse(self):
        for name, axis in self.fAxis.iteritems():
            self.fResponses[name] = self.fProjector.GetDetectorResponse(name, self.fDMeson, self.fJetDefinitions, axis)
        
    def Start(self):
        self.ProjectResponse()
        self.GenerateEfficiency()
        self.PlotResponse()

    def GenerateEfficiency(self):
        for jets in self.fResponses.itervalues():
            for resp in jets.itervalues():
                resp.GenerateEfficiency()

    def PlotResponse(self):
        for jets in self.fResponses.itervalues():
            for resp in jets.itervalues():
                self.PlotResponseMatrix(resp)
                self.PlotEfficiency(resp)

    def PlotResponseMatrix(self, resp):
        if len(resp.fAxis) == 1:
            c = ROOT.TCanvas("{0}_canvas".format(resp.fResponseMatrix.GetName()), resp.fResponseMatrix.GetTitle())
            c.cd()
            c.SetLogz()
            h = resp.fResponseMatrix.DrawCopy("colz")
            globalList.append(h)
            globalList.append(c)

    def PlotEfficiency(self, resp):
        if len(resp.fAxis) == 1:
            c = ROOT.TCanvas("{0}_canvas".format(resp.fEfficiency.GetName()), resp.fEfficiency.GetTitle())
            c.cd()
            h = resp.fEfficiency.Clone()
            h.Draw()
        elif len(resp.fAxis) == 2:
            c = ROOT.TCanvas("{0}_canvas".format(resp.fEfficiency.GetName()), resp.fEfficiency.GetTitle())
            c.cd()
            h = resp.fEfficiency.DrawCopy("colz")
            
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
