#!/usr/bin/env python
#python base classes and utilities for D Meson jet analysis

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
import DMesonJetFDCorrection

class Axis:
    def __init__(self, name, bins, label = "", charged_jet = True):
        self.fName = name
        self.fBins = bins
        self.fLabel = label
        self.fChargedJet = charged_jet

    @classmethod
    def fromLimits(cls, name, start, stop, step, label = "", charged_jet = True):
        bins = []
        bins.extend(numpy.linspace(start, stop, (stop-start)/step+1, True))
        return cls(name, bins, label, charged_jet)

    def FindBin(self, x):
        for i,(min,max) in enumerate(zip(self.fBins[:-1], self.fBins[1:])):
            if x > min and x < max:
                return i
        return -1

    def GetNbins(self):
        return len(self.fBins)-1

    def GetTitle(self, label = ""):
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

    def GetVariableName(self, label = ""):
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

class AnalysisType(Enum):
    InvMassFit = 1
    SideBand = 2
    LikeSign = 3
    LikeSignFit = 4
    Truth = 5

class Spectrum:
    def __init__(self, config, dmeson, jtype, jradius, jtitle, binSet, effWeight, FD):
        self.fDMeson = dmeson
        self.fJetType = jtype
        self.fJetRadius = jradius
        self.fJetTitle = jtitle
        self.fFDCorrection = FD
        if "suffix" in config:
            suffix = config["suffix"]
        else:
            suffix = None 
        self.fName = '_'.join(obj for obj in [self.fDMeson, self.fJetType, self.fJetRadius, config["name"], suffix] if obj)
        self.fBinSet = binSet
        self.fHistogram = None
        self.fNormHistogram = None
        self.fNormFDCorrHistogram = None
        self.fUncertainty = None
        self.fMass = None
        self.fMassWidth = None
        self.fBackground = None
        self.fAxis = []
        self.fSideBandWindowInvMassHistos = dict()
        self.fSignalWindowInvMassHistos = dict()
        self.fSkipBins = None
        self.fEfficiencyWeight = effWeight
        self.fFDCorrHistogram = None
        self.fFDHistogram = None

        # S-B analysis
        self.fSideBandHistograms = None
        self.fSignalHistograms = None
        self.fSideBandWindowTotalHistogram = None
        self.fSignalWindowTotalHistogram = None

        # L-S analysis
        self.fLikeSignSubtractedBinSet = None
        self.fLikeSignHistograms = None
        self.fUnlikeSignHistograms = None
        self.fLikeSignTotalHistogram = None
        self.fUnlikeSignTotalHistogram = None

        self.fTitle = config["title"]
        if binSet.fTitle:
            if self.fTitle:
                self.fTitle += ", "
            self.fTitle += binSet.fTitle

        if "compare" in config:
            self.fCompare = config["compare"]
        else:
            self.fCompare = True

        if config["type"] == "side_band":
            self.fAnalysisType = AnalysisType.SideBand
            self.fSideBandMinSigmas = config["side_band"]["min_sigmas"]
            self.fSideBandMaxSigmas = config["side_band"]["max_sigmas"]
            self.fBinCountSignalSigmas = config["side_band"]["max_signal_sigmas"]
            self.fBackupSigma = config["side_band"]["backup_sigma"]
            self.fBackupMean = config["side_band"]["backup_mean"]
            if "skip_bins" in config["side_band"]:
                self.fSkipBins = config["side_band"]["skip_bins"]
        elif config["type"] == "like_sign":
            if config["like_sign"]["mode"] == "bin_count":
                self.fAnalysisType = AnalysisType.LikeSign
            elif config["like_sign"]["mode"] == "fit":
                self.fAnalysisType = AnalysisType.LikeSignFit
            else:
                print("Like Sign mode {0} not recognized!!".format(config["like_sign"]["mode"]))
                exit(1)
            self.fSideBandMinSigmas = config["like_sign"]["min_sigmas"]
            self.fSideBandMaxSigmas = config["like_sign"]["max_sigmas"]
            self.fBinCountSignalSigmas = config["like_sign"]["max_signal_sigmas"]
            self.fBinCountNormSignalSigmas = config["like_sign"]["max_signal_norm_sigmas"]
            self.fBackupSigma = config["like_sign"]["backup_sigma"]
            self.fBackupMean = config["like_sign"]["backup_mean"]
            self.fLikeSignTree = config["like_sign"]["name"]
            if "skip_bins" in config["like_sign"]:
                self.fSkipBins = config["like_sign"]["skip_bins"]
        elif config["type"] == "inv_mass_fit":
            self.fAnalysisType = AnalysisType.InvMassFit
        elif config["type"] == "truth":
            self.fAnalysisType = AnalysisType.Truth
        else:
            print("Analysis type {0} not recognized!".format(config["type"]))

        if (self.fAnalysisType == AnalysisType.SideBand or self.fAnalysisType == AnalysisType.LikeSign or self.fAnalysisType == AnalysisType.LikeSignFit):
            self.fAxis.append(binSet.fBinCountAnalysisAxis)
        else:
            for axis in binSet.fAxis:
                self.fAxis.append(axis)

        print("Spectrum {0} with analysis type {1} added".format(self.fName, self.fAnalysisType.name))
        self.BuildHistograms()

    def GenerateRootList(self):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)
        if self.fHistogram:
            rlist.Add(self.fHistogram)
        if self.fFDHistogram:
            rlist.Add(self.fFDHistogram)
        if self.fFDCorrHistogram:
            rlist.Add(self.fFDCorrHistogram)
        if self.fNormHistogram:
            rlist.Add(self.fNormHistogram)
        if self.fNormFDCorrHistogram:
            rlist.Add(self.fNormFDCorrHistogram)
        if self.fUncertainty:    
            rlist.Add(self.fUncertainty)
        if self.fMass:    
            rlist.Add(self.fMass)
        if self.fMassWidth:    
            rlist.Add(self.fMassWidth)
        if self.fBackground:    
            rlist.Add(self.fBackground)
        if self.fAnalysisType == AnalysisType.SideBand and self.fSideBandHistograms and self.fSignalHistograms:
            SBlist = ROOT.TList()
            SBlist.SetName("SideBandAnalysis")
            for h in self.fSideBandWindowInvMassHistos.itervalues():
                SBlist.Add(h)
            for h in self.fSignalWindowInvMassHistos.itervalues():
                SBlist.Add(h)
            for h in self.fSideBandHistograms:
                SBlist.Add(h)
            for h in self.fSignalHistograms:
                SBlist.Add(h)
            SBlist.Add(self.fSideBandWindowTotalHistogram)
            SBlist.Add(self.fSignalWindowTotalHistogram)
            rlist.Add(SBlist)
        if self.fAnalysisType == AnalysisType.LikeSign and self.fLikeSignSubtractedBinSet:
            LSlist = self.fLikeSignSubtractedBinSet.GenerateInvMassRootList()
            LSlist.SetName("LikeSignAnalysis")
            for h in self.fLikeSignHistograms:
                LSlist.Add(h)
            for h in self.fUnlikeSignHistograms:
                LSlist.Add(h)
            LSlist.Add(self.fLikeSignTotalHistogram)
            LSlist.Add(self.fUnlikeSignTotalHistogram)
            rlist.Add(LSlist)
        return rlist

    def GenerateFDCorrectedSpectrum(self, events, isWeighted):
        self.fFDCorrHistogram = self.fHistogram.Clone("{0}_FDCorr".format(self.fHistogram.GetName()))
        if self.fFDCorrection.fFDHistogram:
            crossSection = 62.3 #mb CINT1
            branchingRatio = 0.0388 # D0->Kpi
            fdhist = self.fFDCorrection.fFDHistogram
            self.fFDHistogram = fdhist.Clone("{0}_FD".format(self.fHistogram.GetName()))
            if not isWeighted:
                self.fFDHistogram.Scale(events / crossSection * branchingRatio)
            self.fFDCorrHistogram.Add(self.fFDHistogram, -1)

    def GenerateNormalizedSpectrum(self, events, weighted=False):
        if self.fHistogram:
            self.fNormHistogram = self.GenerateNormalizedSpectrumForHistogram(self.fHistogram, events, weighted)
        if self.fFDCorrHistogram:
            self.fNormFDCorrHistogram = self.GenerateNormalizedSpectrumForHistogram(self.fFDCorrHistogram, events, weighted)

    def GenerateNormalizedSpectrumForHistogram(self, hist, events, weighted):
        hname = "{0}_Normalized".format(hist.GetName())
        result = hist.Clone(hname)
        result.SetTitle(hname)
        if len(self.fAxis) == 1:
            if weighted:
                axisTitle = "#frac{{d#sigma}}{{d{var}}}".format(var=self.fAxis[0].GetVariableName())
                if self.fAxis[0].GetVariableUnits():
                    axisTitle += " [mb ({0})^{{-1}}]".format(self.fAxis[0].GetVariableUnits())
                else:
                    axisTitle += " (mb)"
            else:
                axisTitle = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d#it{{N}}}}{{d{var}}}".format(var=self.fAxis[0].GetVariableName())
                if self.fAxis[0].GetVariableUnits():
                    axisTitle += " ({0})^{{-1}}".format(self.fAxis[0].GetVariableUnits())
            result.GetYaxis().SetTitle(axisTitle)
        elif len(self.fAxis) == 2:
            units1 = self.fAxis[0].GetVariableUnits()
            units2 = self.fAxis[1].GetVariableUnits()
            if weighted:
                axisTitle = "#frac{{d^{{2}}#sigma}}{{d{var1} d{var2}}}".format(var1=self.fAxis[0].GetVariableName(), var2=self.fAxis[1].GetVariableName())
                if units1 != units2:
                    axisTitle += " [mb"
                    if units1:
                        axisTitle += " ({0})^{{-1}}".format(units1)
                    if units2:
                        axisTitle += " ({0})^{{-1}}".format(units2)
                    axisTitle += "]"
                else:
                    if units1:
                        axisTitle += " [mb ({0})^{{-2}}]".format(units1)
            else:
                axisTitle = "#frac{{1}}{{#it{{N}}_{{evt}}}} #frac{{d^{{2}}#it{{N}}}}{{d{var1} d{var2}}}".format(var1=self.fAxis[0].GetVariableName(), var2=self.fAxis[1].GetVariableName())
                if units1 != units2:
                    if units1:
                        axisTitle += " ({0})^{{-1}}".format(units1)
                    if units2:
                        axisTitle += " ({0})^{{-1}}".format(units2)
                else:
                    if units1:
                        axisTitle += " ({0})^{{-2}}".format(units1)
            result.GetZaxis().SetTitle(axisTitle)

        if not weighted and events > 0:
            result.Scale(1. / events, "width")
        else:
            result.Scale(1, "width")
        return result

    def BuildHistograms(self):
        if len(self.fAxis) == 1:
            self.fHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, self.fName, "counts")
            self.fUncertainty = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Unc".format(self.fName), "relative statistical uncertainty")
            if self.fAnalysisType == AnalysisType.SideBand:
                self.fSideBandHistograms = []
                self.fSignalHistograms = []
                self.fMass = None
                self.fMassWidth = None
                self.fBackground = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Bkg".format(self.fName), "background |#it{{m}} - <#it{{m}}>| < {0}#sigma".format(int(self.fBinCountSignalSigmas)))
                self.fSideBandWindowTotalHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_SideBandWindowTotal".format(self.fName), "counts")
                self.fSignalWindowTotalHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_SignalWindowTotal".format(self.fName), "counts")
            elif self.fAnalysisType == AnalysisType.LikeSign:
                self.fLikeSignHistograms = []
                self.fUnlikeSignHistograms = []
                self.fMass = None
                self.fMassWidth = None
                self.fBackground = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Bkg".format(self.fName), "background |#it{{m}} - <#it{{m}}>| < {0}#sigma".format(int(self.fBinCountSignalSigmas)))
                self.fLikeSignTotalHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_LikeSignTotal".format(self.fName), "counts")
                self.fUnlikeSignTotalHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_UnlikeSignTotal".format(self.fName), "counts")
            elif self.fAnalysisType == AnalysisType.InvMassFit or self.fAnalysisType == AnalysisType.LikeSignFit:
                self.fMass = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Mass".format(self.fName), "D^{0} mass (GeV/#it{c}^{2})")
                self.fMassWidth = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_MassWidth".format(self.fName), "D^{0} mass width (GeV/#it{c}^{2})")
                self.fBackground = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Bkg".format(self.fName), "background |#it{m} - <#it{m}>| < 2#sigma")
        elif len(self.fAxis) == 2:
            self.fHistogram = DMesonJetUtils.BuildHistogram(self.fAxis, self.fName, "counts")
            self.fUncertainty = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Unc".format(self.fName), "relative statistical uncertainty")
            if self.fAnalysisType == AnalysisType.InvMassFit:
                self.fMass = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Mass".format(self.fName), "mass")
                self.fBackground = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Bkg".format(self.fName), "background |#it{m} - <#it{m}>| < 3#sigma")
        else:
            print("Cannot handle spectra with more than two axis!")
            return

    def Fill(self, dmeson, jet, w):
        values = []
        for axis in self.fAxis:
            if axis.fName == "jet_pt":
                values.append(jet.fPt)
            elif axis.fName == "jet_eta":
                values.append(jet.fEta)
            elif axis.fName == "jet_phi":
                values.append(jet.fPhi)
            elif axis.fName == "d_z":
                values.append(jet.fZ)
            elif axis.fName == "d_pt":
                values.append(dmeson.fPt)
            elif axis.fName == "d_eta":
                values.append(dmeson.fEta)
            elif axis.fName == "d_phi":
                values.append(dmeson.fPhi)
        if len(values) == 1:
            self.fHistogram.Fill(values[0], w)
        elif len(values) == 2:
            self.fHistogram.Fill(values[0], values[1], w)
        else:
            print("Cannot handle histograms with more than two axis!")

class DMesonJetCuts:
    def __init__(self, cutList):
        self.fCuts = cutList
        self.InitializeCuts()

    def Dummy(self, obj, min, max):
        return True

    def CutJet(self, dmeson, jet, varname, min, max):
        v = getattr(jet, varname)
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutJetPt(self, dmeson, jet, dummy, min, max):
        v = jet.fPt
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutJetEta(self, dmeson, jet, dummy, min, max):
        v = jet.fEta
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutJetPhi(self, dmeson, jet, dummy, min, max):
        v = jet.fPhi
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutJetZ(self, dmeson, jet, dummy, min, max):
        v = jet.fZ
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutD(self, dmeson, jet, varname, min, max): 
        v = getattr(dmeson, varname)
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutDPt(self, dmeson, jet, dummy, min, max):
        v = dmeson.fPt
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutDEta(self, dmeson, jet, dummy, min, max):
        v = dmeson.fEta
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutDPhi(self, dmeson, jet, dummy, min, max):
        v = dmeson.fPhi
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def InitializeCuts(self):
        self.fFastCuts = []
        for cut in self.fCuts:
            if "min" in cut: min = cut["min"]
            else: min = -sys.float_info.max
            if "max" in cut: max = cut["max"]
            else: max = +sys.float_info.max
            if cut["object"] == "d":
                if cut["variable"] == "fPt":
                    self.fFastCuts.append((self.CutDPt, None, min, max))
                elif cut["variable"] == "fEta":
                    self.fFastCuts.append((self.CutDEta, None, min, max))
                elif cut["variable"] == "fPhi":
                    self.fFastCuts.append((self.CutDPhi, None, min, max))
                else:
                    self.fFastCuts.append((self.CutD, cut["variable"], min, max))
            elif cut["object"] == "jet":
                if cut["variable"] == "fPt":
                    self.fFastCuts.append((self.CutJetPt, None, min, max))
                elif cut["variable"] == "fEta":
                    self.fFastCuts.append((self.CutJetEta, None, min, max))
                elif cut["variable"] == "fPhi":
                    self.fFastCuts.append((self.CutJetPhi, None, min, max))
                elif cut["variable"] == "fZ":
                    self.fFastCuts.append((self.CutJetZ, None, min, max))
                else:
                    self.fFastCuts.append((self.CutJet, cut["variable"], min, max))

    def ApplyCuts(self, dmeson, jet):
        for funct,var,min,max in self.fFastCuts:
            if not funct(dmeson, jet, var, min, max): return False
        return True

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
        self.fCuts = DMesonJetCuts(cutList)
        self.fFitOptions = fitOptions
        self.fWeightEfficiency = weight
        self.fNeedInvMass = need_inv_mass
        if bin_count_axis:
            self.fBinCountAnalysisAxis = Axis(bin_count_axis.keys()[0], bin_count_axis.values()[0], "", True)
        else:
            self.fBinCountAnalysisAxis = None
        limits = dict()
        self.AddBinsRecursive(limitSetList, limits)

    def AmIJetDependent(self):
        for axis in self.fAxis:
            if "jet" in axis.fName or axis.fName == "d_z":
                return True
        if self.fBinCountAnalysisAxis and ("jet" in self.fBinCountAnalysisAxis.fName or self.fBinCountAnalysisAxis.fName == "d_z"):
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
            for a in self.fAxis:
                a.fChargedJet = False
        for s in self.fSpectraConfigs:
            if not dmeson in s["active_mesons"]:
                continue
            if "efficiency" in s and s["efficiency"] and not "MCTruth" in dmeson:
                effWeight = DMesonJetProjectors.EfficiencyWeightCalculator("{0}/{1}".format(inputPath, s["efficiency"]["file_name"]), s["efficiency"]["list_name"], s["efficiency"]["object_name"])
            else:
                effWeight = DMesonJetProjectors.SimpleWeight()
            if "FD" in s and s["FD"]:
                FD = DMesonJetFDCorrection.DMesonJetFDCorrection(s["FD"], inputPath, dmeson, jtype, jradius)
            else:
                FD = DMesonJetFDCorrection.DMesonJetFDCorrection(None)
            spectrum = Spectrum(s, dmeson, jtype, jradius, jtitle, self, effWeight, FD)
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
        return rlist

    def AddBinsRecursive(self, limitSetList, limits):
        if len(limitSetList) > 0:
            (limitSetName, limitSet) = limitSetList[0]
            for min,max in zip(limitSet[:-1], limitSet[1:]):
                limits[limitSetName] = min,max
                self.AddBinsRecursive(limitSetList[1:], limits)
        else:
            bin = BinLimits(limits)
            bin.fBinCountAnalysisAxis = self.fBinCountAnalysisAxis
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
        self.fBinCountAnalysisHisto = None
        self.fCounts = 0
        self.fSumw2 = 0

    def FillInvariantMass(self, dmeson, jet, w):
        self.fCounts += w
        self.fSumw2 += w*w
        if self.fInvMassHisto:
            self.fInvMassHisto.Fill(dmeson.fInvMass, w)
        if self.fBinCountAnalysisHisto:
            if self.fBinCountAnalysisAxis.fName == "jet_pt":
                obsVal = jet.fPt
            else:
                obsVal = 0
            self.fBinCountAnalysisHisto.Fill(dmeson.fInvMass, obsVal, w)
        
    def AddFromAxis(self, axis, binIndex):
        if binIndex == 0:
            min = axis.fBins[0]
            max = axis.fBins[-1]
        else:
            min = axis.fBins[binIndex-1]
            max = axis.fBins[binIndex]                
        self.fLimits[axis.fName] = min, max
      
    def SetMassFitter(self, fitter):
        self.fMassFitter = fitter
    
    def IsSameOf(self, bin):
        for name,(min, max) in self.fLimits.iteritems():
            if bin.fLimits.has_key(name):
                if not bin.fLimits[name] == (min, max):
                    return False
            else:
                return False
        return True

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
        for name,(min,max) in self.fLimits.iteritems():
            if not min < max:
                continue
            if name == "d_pt" and (dmeson.fPt < min or dmeson.fPt >= max):
                return False
            elif name == "jet_pt" and (jet.fPt < min or jet.fPt >= max):
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
        for varName,(min,max) in self.fLimits.iteritems():
            if varName == "d_pt":
                name += "DPt_{0}_{1}_".format(int(min*100), int(max*100))
            elif varName == "jet_pt":
                name += "JetPt_{0}_{1}_".format(int(min*100), int(max*100))
            elif varName == "d_eta":
                name += "DEta_{0}_{1}_".format(int(min*10), int(max*10))
            elif varName == "jet_eta":
                name += "JetEta_{0}_{1}_".format(int(min*10), int(max*10))
            elif varName == "d_z":
                name += "DZ_{0}_{1}_".format(int(min*100), int(max*100))
        
        #remove last "_"
        if name:
            name = name[:-1]
        return name
        
    def GetTitle(self):
        title = ""
        
        for varName,(min,max) in self.fLimits.iteritems():
            if varName == "d_pt":
                title += "{0:.1f} < #it{{p}}_{{T,D}} < {1:.1f} GeV/#it{{c}}, ".format(min, max)
            elif varName == "jet_pt":
                title += "{0:.0f} < #it{{p}}_{{T,jet}} < {1:.0f} GeV/#it{{c}}, ".format(min, max)
            elif varName == "d_eta":
                title += "{0:.1f} < #it{{#eta}}_{{D}} < {1:.1f}, ".format(min, max)
            elif varName == "jet_eta":
                title += "{0:.0f} < #it{{#eta}}_{{jet}} < {1:.0f}, ".format(min, max)
            elif varName == "d_z":
                title += "{0:.1f} < #it{{z}}_{{||, D}} < {1:.1f}, ".format(min, max)

        #remove last ", "
        if title:
            title = title[:-2]
        return title
    
    def Print(self):
        print(self.GetTitle())
    
    def CreateInvMassHisto(self, trigger, DMesonDef, xAxis, yAxis, nMassBins, minMass, maxMass):
        if trigger:
            hname = "InvMass_{0}_{1}_{2}".format(trigger, DMesonDef, self.GetName())
            htitle = "{0} - {1} Invariant Mass: {2};{3};{4}".format(trigger, DMesonDef, self.GetTitle(), xAxis, yAxis)
            if self.fBinCountAnalysisAxis:
                hnameSB = "InvMassBinCounting_{0}_{1}_{2}".format(trigger, DMesonDef, self.GetName())
                htitleSB = "{0} - {1} Invariant Mass: {2};{3};{4};{5}".format(trigger, DMesonDef, self.GetTitle(), xAxis, self.fBinCountAnalysisAxis.GetTitle(), yAxis)
        else:
            hname = "InvMass_{0}_{1}".format(DMesonDef, self.GetName())
            htitle = "{0} Invariant Mass: {1};{2};{3}".format(DMesonDef, self.GetTitle(), xAxis, yAxis)
            if self.fBinCountAnalysisAxis:
                hnameSB = "InvMassBinCounting_{0}_{1}".format(DMesonDef, self.GetName())
                htitleSB = "{0} Invariant Mass: {1};{2};{3};{4}".format(DMesonDef, self.GetTitle(), xAxis, self.fBinCountAnalysisAxis.GetTitle(), yAxis)
        
        if not "MCTruth" in DMesonDef:
            self.fInvMassHisto = ROOT.TH1D(hname, htitle, nMassBins, minMass, maxMass)
            self.fInvMassHisto.Sumw2()
            self.fInvMassHisto.SetMarkerSize(0.9)
            self.fInvMassHisto.SetMarkerStyle(ROOT.kFullCircle)
            self.fInvMassHisto.SetMarkerColor(ROOT.kBlue+2)
            self.fInvMassHisto.SetLineColor(ROOT.kBlue+2)
        
        if self.fBinCountAnalysisAxis:
            self.fBinCountAnalysisHisto = ROOT.TH2D(hnameSB, htitleSB, nMassBins, minMass, maxMass, len(self.fBinCountAnalysisAxis.fBins)-1, array.array('d', self.fBinCountAnalysisAxis.fBins))
            self.fBinCountAnalysisHisto.Sumw2()
