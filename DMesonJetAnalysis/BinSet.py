#!/usr/bin/env python
# python base classes and utilities for D Meson jet analysis

import array
import copy
import collections
import numpy

import ROOT

import DMesonJetProjectors
import DMesonJetCuts
import Axis
import Spectrum
from DMesonJetBase import AnalysisType

class BinMultiSet:
    def __init__(self):
        self.fBinSets = collections.OrderedDict()

    def GenerateInvMassRootLists(self):
        for binSet in self.fBinSets.itervalues():
            yield binSet.GenerateInvMassRootList()

    def AddBinSet(self, binSet):
        self.fBinSets[binSet.fBinSetName] = binSet

    def Initialize(self, dmeson, jtype, jradius, jtitle, inputPath):
        self.fDMeson = dmeson
        self.fJetType = jtype
        self.fJetRadius = jradius
        for k, binSet in self.fBinSets.items():
            r = binSet.Initialize(dmeson, jtype, jradius, jtitle, inputPath)
            if not r:
                del self.fBinSets[k]

    def SetWeightEfficiency(self, weight):
        for binSet in self.fBinSets.itervalues():
            binSet.fWeightEfficiency = weight

    def FindBin(self, dmeson, jet, dmesonDef):
        for binSet in self.fBinSets.itervalues():
            if not dmesonDef in binSet.fNeedInvMass:
                continue
            bins = binSet.FindBin(dmeson, jet)
            for bin, w in bins:
                yield bin, w

    def FindSpectra(self, dmeson, jet):
        for binSet in self.fBinSets.itervalues():
            if not binSet.fCuts.ApplyCuts(dmeson, jet):
                continue
            for s in binSet.fSpectra.itervalues():
                if s.fAnalysisType != AnalysisType.Truth:
                    continue
                yield s, binSet.fWeightEfficiency.GetEfficiencyWeight(dmeson, jet)

    def FindAllSpectra(self):
        for binSet in self.fBinSets.itervalues():
            for s in binSet.fSpectra.itervalues():
                yield s

class BinSet:
    def __init__(self, name, title, need_inv_mass, limitSetList, spectra, axis, cutList=[], bin_count_axis=None, weight=None, fitOptions=""):
        self.fBinSetName = name
        self.fTitle = title
        self.fBins = []
        self.fSpectraConfigs = spectra
        self.fSpectra = collections.OrderedDict()
        self.fAxis = axis
        self.fCuts = DMesonJetCuts.DMesonJetCuts(cutList)
        self.fFitOptions = fitOptions
        self.fWeightEfficiency = weight
        self.fNeedInvMass = need_inv_mass
        self.fBinCountAnalysisAxis = None
        self.fBinCountAnalysisSecondAxis = None
        if bin_count_axis:
            self.fBinCountAnalysisAxis = Axis.Axis(bin_count_axis.keys()[0], bin_count_axis.values()[0], "", True)
            if len(bin_count_axis) > 1:
                self.fBinCountAnalysisSecondAxis = Axis.Axis(bin_count_axis.keys()[1], bin_count_axis.values()[1], "", True)
        limits = dict()
        self.AddBinsRecursive(limitSetList, limits)

    def AmIJetDependent(self):
        for axis in self.fAxis:
            if "jet" in axis.fName or axis.fName == "d_z":
                return True
        if self.fBinCountAnalysisAxis and ("jet" in self.fBinCountAnalysisAxis.fName or self.fBinCountAnalysisAxis.fName == "d_z"):
            return True
        if self.fBinCountAnalysisSecondAxis and ("jet" in self.fBinCountAnalysisSecondAxis.fName or self.fBinCountAnalysisSecondAxis.fName == "d_z"):
            return True
        for cut in self.fCuts.fCuts:
            if cut["object"] == "jet":
                return True
        return False

    def Initialize(self, dmeson, jtype, jradius, jtitle, inputPath):
        self.fName = "_".join(obj for obj in [dmeson, jtype, jradius, self.fBinSetName] if not obj == None)
        jetDep = self.AmIJetDependent()
        if jtype or jradius:
            if not jetDep:
                return False
        else:
            if jetDep:
                return False
        if jtype == "Full":
            if self.fBinCountAnalysisAxis:
                self.fBinCountAnalysisAxis.fChargedJet = False
            if self.fBinCountAnalysisSecondAxis:
                self.fBinCountAnalysisSecondAxis.fChargedJet = False
            for a in self.fAxis:
                a.fChargedJet = False
        for s in self.fSpectraConfigs:
            if not dmeson in s["active_mesons"]:
                continue
            if "efficiency" in s and s["efficiency"] and not "MCTruth" in dmeson:
                effWeight = DMesonJetProjectors.EfficiencyWeightCalculator("{0}/{1}".format(inputPath, s["efficiency"]["file_name"]), s["efficiency"]["list_name"], s["efficiency"]["object_name"])
            else:
                effWeight = DMesonJetProjectors.SimpleWeight()
            spectrum = Spectrum.Spectrum(s, dmeson, jtype, jradius, jtitle, self, effWeight)
            self.fSpectra[spectrum.fName] = spectrum
        if "MCTruth" in dmeson and not isinstance(self.fWeightEfficiency , DMesonJetProjectors.SimpleWeight):
            self.fWeightEfficiency = DMesonJetProjectors.SimpleWeight()
        return True

    def GenerateInvMassRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        for bin in self.fBins:
            if bin.fInvMassHisto:
                rlist.Add(bin.fInvMassHisto)
            if bin.fMassFitter:
                rlist.Add(bin.fMassFitter)
            if bin.fBinCountAnalysisHisto:
                rlist.Add(bin.fBinCountAnalysisHisto)
        return rlist

    def AddBinsRecursive(self, limitSetList, limits):
        if len(limitSetList) > 0:
            (limitSetName, limitSet) = limitSetList[0]
            for min, max in zip(limitSet[:-1], limitSet[1:]):
                limits[limitSetName] = min, max
                self.AddBinsRecursive(limitSetList[1:], limits)
        else:
            bin = BinLimits(limits)
            bin.fBinCountAnalysisAxis = self.fBinCountAnalysisAxis
            bin.fBinCountAnalysisSecondAxis = self.fBinCountAnalysisSecondAxis
            self.fBins.append(bin)

    def FindBin(self, dmeson, jet):
        if not self.fCuts.ApplyCuts(dmeson, jet):
            return
        else:
            w = self.fWeightEfficiency.GetEfficiencyWeight(dmeson, jet)
        for bin in self.fBins:
            if bin.IsInBinLimits(dmeson, jet):
                yield bin, w

class BinLimits:
    def __init__(self, limits=dict()):
        self.fLimits = copy.deepcopy(limits)
        self.fInvMassHisto = None
        self.fMassFitter = None
        self.fBinCountAnalysisAxis = None
        self.fBinCountAnalysisSecondAxis = None
        self.fBinCountAnalysisHisto = None
        self.fCounts = 0
        self.fSumw2 = 0

    def Fill(self, dmeson, jet, w):
        self.fCounts += w
        self.fSumw2 += w * w
        if self.fInvMassHisto:
            self.fInvMassHisto.Fill(dmeson.fInvMass, w)

        if self.fBinCountAnalysisHisto:
            obsValX = 0
            obsValY = 0
            if self.fBinCountAnalysisAxis:
                if self.fBinCountAnalysisAxis.fName == "jet_pt":
                    obsValX = jet.fPt
                elif self.fBinCountAnalysisAxis.fName == "d_z":
                    obsValX = jet.fZ
                elif self.fBinCountAnalysisAxis.fName == "jet_n":
                    obsValX = jet.fN
                elif self.fBinCountAnalysisAxis.fName == "jet_corrpt":
                    obsValX = jet.fCorrPt
            if self.fBinCountAnalysisSecondAxis:
                if self.fBinCountAnalysisSecondAxis.fName == "jet_pt":
                    obsValY = jet.fPt
                elif self.fBinCountAnalysisSecondAxis.fName == "d_z":
                    obsValY = jet.fZ
                elif self.fBinCountAnalysisSecondAxis.fName == "jet_n":
                    obsValY = jet.fN
                elif self.fBinCountAnalysisSecondAxis.fName == "jet_corrpt":
                    obsValX = jet.fCorrPt

            if self.fBinCountAnalysisHisto.GetDimension() == 2:
                self.fBinCountAnalysisHisto.Fill(dmeson.fInvMass, obsValX, w)
            elif self.fBinCountAnalysisHisto.GetDimension() == 3:
                self.fBinCountAnalysisHisto.Fill(dmeson.fInvMass, obsValX, obsValY, w)

    def AddFromAxis(self, axis, binIndex):
        if binIndex == 0:
            min = axis.fBins[0]
            max = axis.fBins[-1]
        else:
            min = axis.fBins[binIndex - 1]
            max = axis.fBins[binIndex]
        self.fLimits[axis.fName] = min, max

    def SetMassFitter(self, fitter):
        self.fMassFitter = fitter

    def IsSameOf(self, bin):
        for name, (min, max) in self.fLimits.iteritems():
            if bin.fLimits.has_key(name):
                if not bin.fLimits[name] == (min, max):
                    return False
            else:
                return False
        return True

    def SetJetCorrPtLimits(self, min, max):
        self.fLimits["jet_corrpt"] = min, max

    def SetJetPtLimits(self, min, max):
        self.fLimits["jet_pt"] = min, max

    def SetDPtLimits(self, min, max):
        self.fLimits["d_pt"] = min, max

    def SetDZLimits(self, min, max):
        self.fLimits["d_z"] = min, max

    def SetJetEtaLimits(self, min, max):
        self.fLimits["jet_eta"] = min, max

    def SetDEtaLimits(self, min, max):
        self.fLimits["d_eta"] = min, max

    def IsInBinLimits(self, dmeson, jet):
        for name, (min, max) in self.fLimits.iteritems():
            if not min < max:
                continue
            if name == "d_pt" and (dmeson.fPt < min or dmeson.fPt >= max):
                return False
            elif name == "jet_pt" and (jet.fPt < min or jet.fPt >= max):
                return False
            elif name == "jet_corrpt" and (jet.fCorrPt < min or jet.fCorrPt >= max):
                return False
            elif name == "d_eta" and (dmeson.fEta < min or dmeson.fEta >= max):
                return False
            elif name == "jet_eta" and (jet.fEta < min or jet.fEta >= max):
                return False
            elif name == "d_z" and (jet.fZ < min or jet.fZ >= max):
                return False

        return True

    def GetBinCenter(self, axis):
        if axis in self.fLimits:
            (min, max) = self.fLimits[axis]
            return (min + max) / 2
        else:
            return -1

    def GetName(self):
        name = ""
        for varName, (min, max) in self.fLimits.iteritems():
            if varName == "d_pt":
                name += "DPt_{0}_{1}_".format(int(min * 100), int(max * 100))
            elif varName == "jet_pt":
                name += "JetPt_{0}_{1}_".format(int(min * 100), int(max * 100))
            elif varName == "jet_corrpt":
                name += "JetCorrPt_{0}_{1}_".format(int(min * 100), int(max * 100)).replace("-", "N")
            elif varName == "d_eta":
                name += "DEta_{0}_{1}_".format(int(min * 10), int(max * 10)).replace("-", "N")
            elif varName == "jet_eta":
                name += "JetEta_{0}_{1}_".format(int(min * 10), int(max * 10)).replace("-", "N")
            elif varName == "d_z":
                name += "DZ_{0}_{1}_".format(int(min * 100), int(max * 100))

        # remove last "_"
        if name: name = name[:-1]

        return name

    def GetTitle(self):
        title = ""

        for varName, (min, max) in self.fLimits.iteritems():
            if varName == "d_pt":
                title += "{0:.1f} < #it{{p}}_{{T,D}} < {1:.1f} GeV/#it{{c}}, ".format(min, max)
            elif varName == "jet_pt":
                title += "{0:.0f} < #it{{p}}_{{T,jet}} < {1:.0f} GeV/#it{{c}}, ".format(min, max)
            elif varName == "jet_corrpt":
                title += "{0:.0f} < #it{{p}}_{{T,jet}}^{{sub}} < {1:.0f} GeV/#it{{c}}, ".format(min, max)
            elif varName == "d_eta":
                title += "{0:.1f} < #it{{#eta}}_{{D}} < {1:.1f}, ".format(min, max)
            elif varName == "jet_eta":
                title += "{0:.0f} < #it{{#eta}}_{{jet}} < {1:.0f}, ".format(min, max)
            elif varName == "d_z":
                title += "{0:.1f} < #it{{z}}_{{||, D}} < {1:.1f}, ".format(min, max)

        # remove last ", "
        if title: title = title[:-2]
        return title

    def Print(self):
        print(self.GetTitle())

    def CreateHistograms(self, trigger, DMesonDef, xAxis, yAxis, nMassBins, minMass, maxMass):
        self.CreateInvMassHisto(trigger, DMesonDef, xAxis, yAxis, nMassBins, minMass, maxMass)
        self.CreateQAHistos(trigger, DMesonDef, yAxis)

    def CreateQAHistos(self, trigger, DMesonDef, yAxis):
        pass

    def CreateInvMassHisto(self, trigger, DMesonDef, xAxis, yAxis, nMassBins, minMass, maxMass):
        if trigger:
            hname = "InvMass_{0}_{1}_{2}".format(trigger, DMesonDef, self.GetName())
            htitle = "{0} - {1} Invariant Mass: {2};{3};{4}".format(trigger, DMesonDef, self.GetTitle(), xAxis, yAxis)
            if self.fBinCountAnalysisAxis:
                hnameSB = "InvMassBinCounting_{0}_{1}_{2}".format(trigger, DMesonDef, self.GetName())
                if self.fBinCountAnalysisSecondAxis:
                    htitleSB = "{0} - {1} Invariant Mass: {2};{3};{4};{5}".format(trigger, DMesonDef, self.GetTitle(), xAxis, self.fBinCountAnalysisAxis.GetTitle(), self.fBinCountAnalysisSecondAxis.GetTitle())
                else:
                    htitleSB = "{0} - {1} Invariant Mass: {2};{3};{4};{5}".format(trigger, DMesonDef, self.GetTitle(), xAxis, self.fBinCountAnalysisAxis.GetTitle(), yAxis)
        else:
            hname = "InvMass_{0}_{1}".format(DMesonDef, self.GetName())
            htitle = "{0} Invariant Mass: {1};{2};{3}".format(DMesonDef, self.GetTitle(), xAxis, yAxis)
            if self.fBinCountAnalysisAxis:
                hnameSB = "InvMassBinCounting_{0}_{1}".format(DMesonDef, self.GetName())
                if self.fBinCountAnalysisSecondAxis:
                    htitleSB = "{0} Invariant Mass: {1};{2};{3};{4}".format(DMesonDef, self.GetTitle(), xAxis, self.fBinCountAnalysisAxis.GetTitle(), self.fBinCountAnalysisSecondAxis.GetTitle())
                else:
                    htitleSB = "{0} Invariant Mass: {1};{2};{3}".format(DMesonDef, self.GetTitle(), xAxis, self.fBinCountAnalysisAxis.GetTitle(), yAxis)

        if not "MCTruth" in DMesonDef:
            self.fInvMassHisto = ROOT.TH1D(hname, htitle, nMassBins, minMass, maxMass)
            self.fInvMassHisto.Sumw2()
            self.fInvMassHisto.SetMarkerSize(0.9)
            self.fInvMassHisto.SetMarkerStyle(ROOT.kFullCircle)
            self.fInvMassHisto.SetMarkerColor(ROOT.kBlue + 2)
            self.fInvMassHisto.SetLineColor(ROOT.kBlue + 2)

        if self.fBinCountAnalysisAxis:
            if self.fBinCountAnalysisSecondAxis:
                massAxis = numpy.linspace(minMass, maxMass, nMassBins + 1, True)
                self.fBinCountAnalysisHisto = ROOT.TH3D(hnameSB, htitleSB, nMassBins, array.array('d', massAxis), len(self.fBinCountAnalysisAxis.fBins) - 1, array.array('d', self.fBinCountAnalysisAxis.fBins), len(self.fBinCountAnalysisSecondAxis.fBins) - 1, array.array('d', self.fBinCountAnalysisSecondAxis.fBins))
                self.fBinCountAnalysisHisto.Sumw2()
            else:
                self.fBinCountAnalysisHisto = ROOT.TH2D(hnameSB, htitleSB, nMassBins, minMass, maxMass, len(self.fBinCountAnalysisAxis.fBins) - 1, array.array('d', self.fBinCountAnalysisAxis.fBins))
                self.fBinCountAnalysisHisto.Sumw2()
