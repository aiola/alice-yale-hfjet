#!/usr/bin/env python
# python base classes and utilities for D Meson jet analysis

import ROOT
import DMesonJetUtils
from DMesonJetBase import AnalysisType

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
            crossSection = 62.3  # mb CINT1
            branchingRatio = 0.0393  # D0->Kpi
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
            elif self.fAnalysisType == AnalysisType.Truth and self.fDMeson in self.fBinSet.fNeedInvMass:
                self.fMass = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_Mass".format(self.fName), "D^{0} mass (GeV/#it{c}^{2})")
                self.fMassWidth = DMesonJetUtils.BuildHistogram(self.fAxis, "{0}_MassWidth".format(self.fName), "D^{0} mass width (GeV/#it{c}^{2})")
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
