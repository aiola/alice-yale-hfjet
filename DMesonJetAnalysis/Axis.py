#!/usr/local/bin/python
# python base classes and utilities for D Meson jet analysis

import ROOT
import array
import os
import math
import copy
import DMesonJetUtils
import DMesonJetProjectors
import numpy
from enum import Enum
import collections
import sys

class Axis:
    def __init__(self, name, bins, label="", charged_jet=True):
        self.fName = name
        self.fBins = bins
        self.fLabel = label
        self.fChargedJet = charged_jet

    @classmethod
    def fromLimits(cls, name, start, stop, step, label="", charged_jet=True):
        bins = []
        bins.extend(numpy.linspace(start, stop, (stop - start) / step + 1, True))
        return cls(name, bins, label, charged_jet)

    def FindBin(self, x):
        min = self.fBins[0]
        for i in xrange(1, len(self.fBins)):
            max = self.fBins[i]
            if x >= min and x < max: return i - 1
            min = max
        return -1

    def GetNbins(self):
        return len(self.fBins) - 1

    def GetTitle(self, label=""):
        varName = self.GetVariableName(label)
        units = self.GetVariableUnits()

        if units:
            title = "{0} ({1})".format(varName, units)
        else:
            title = varName

        return title

    def GetVariableUnits(self):
        if "pt" in self.fName:
            return "GeV/#it{c}"
        else:
            return ""

    def GetVariableName(self, label=""):
        if not label:
            label = self.fLabel

        if label == "nolabel":
            label = ""

        if self.fChargedJet:
            jetLabel = "ch jet"
        else:
            jetLabel = "jet"

        if self.fName == "jet_pt":
            if label:
                title = "#it{{p}}_{{T,{jetLabel}}}^{{{label}}}".format(jetLabel=jetLabel, label=label)
            else:
                title = "#it{{p}}_{{T,{jetLabel}}}".format(jetLabel=jetLabel)
        elif self.fName == "jet_corrpt":
            if label:
                title = "#it{{p}}_{{T,{jetLabel}}}^{{sub,{label}}}".format(jetLabel=jetLabel, label=label)
            else:
                title = "#it{{p}}_{{T,{jetLabel}}}^{{sub}}".format(jetLabel=jetLabel)
        elif self.fName == "jet_n":
            title = "N jet constituents"
        elif self.fName == "d_pt":
            if label:
                title = "#it{{p}}_{{T,D}}^{{{0}}}".format(label)
            else:
                title = "#it{p}_{T,D}"
        elif self.fName == "jet_eta":
            if label:
                title = "#it{{#eta}}_{{jet}}^{{{0}}}".format(label)
            else:
                title = "#it{#eta}_{jet}"
        elif self.fName == "d_eta":
            if label:
                title = "#it{{#eta}}_{{D}}^{{{0}}}".format(label)
            else:
                title = "#it{#eta}_{D}"
        elif self.fName == "d_z":
            if label:
                title = "#it{{z}}_{{||,D}}^{{{jetLabel},{label}}}".format(jetLabel=jetLabel, label=label)
            else:
                title = "#it{{z}}_{{||,D}}^{{{jetLabel}}}".format(jetLabel=jetLabel)

        return title
