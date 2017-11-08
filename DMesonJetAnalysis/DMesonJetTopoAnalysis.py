#!/usr/bin/env python
# python program to perform a D meson jet analysis

import math
import copy
import collections
import os

import ROOT

import DMesonJetProjectors
import DMesonJetCompare
import DMesonJetUtils
from DMesonJetBase import AnalysisType
import BinSet
import Axis
import collections
import numpy

globalList = []

class DMesonJetVariable:
    def __init__(self, name, varname, vartitle, bins):
        self.fName = name
        self.fVariableName = varname
        self.fBins = bins
        self.fXaxisTitle = vartitle
        self.fHistogram = ROOT.TH1F(self.fName, "{};{};counts".format(self.fName, self.fXaxisTitle), len(bins) - 1, bins)
        self.fHistogram.Sumw2()
        self.fIntegral = None
        self.fSignificance = None
        self.fFraction = None
        self.fMaxSignificance = None
        self.fCutValueAtMaxSignificance = None
        self.fMaxSigOverBkg = None
        self.fCutValueAtMaxSigOverBkg = None

    def FillStd(self, dmeson, w):
        v = getattr(dmeson, self.fVariableName)
        self.fHistogram.Fill(v, w)

    def FillAbs(self, dmeson, w):
        v = math.fabs(getattr(dmeson, self.fVariableName))
        self.fHistogram.Fill(v, w)

    def Analyze(self):
        self.fCutEfficiency = self.fHistogram.Clone("{}_CutEfficiency".format(self.fName))
        self.fCutEfficiency.Reset()
        self.fCutCounts = self.fHistogram.Clone("{}_CutCounts".format(self.fName))
        self.fCutCounts.Reset()
        if self.CalculateCutEfficiency == self.CalculateLeftCutEfficiency:
            self.fCutEfficiency.GetYaxis().SetTitle("Efficiency for cut < {}".format(self.fCutEfficiency.GetXaxis().GetTitle()))
            self.fCutCounts.GetYaxis().SetTitle("Counts for cut < {}".format(self.fCutEfficiency.GetXaxis().GetTitle()))
        else:
            self.fCutEfficiency.GetYaxis().SetTitle("Efficiency for cut > {}".format(self.fCutEfficiency.GetXaxis().GetTitle()))
            self.fCutCounts.GetYaxis().SetTitle("Counts for cut > {}".format(self.fCutEfficiency.GetXaxis().GetTitle()))
        self.fIntegral = self.fHistogram.Integral(0, -1)
        if self.fIntegral == 0: return
        self.CalculateCutEfficiency()

    def CalculateLeftCutEfficiency(self):
        partialIntErr2 = 0.
        partialInt = 0.
        for ibin in xrange(0, self.fHistogram.GetNbinsX() + 2):
            partialIntErr2 += self.fHistogram.GetBinError(ibin) ** 2
            partialInt += self.fHistogram.GetBinContent(ibin)
            eff = partialInt / self.fIntegral
            effErr = math.sqrt(partialIntErr2) / self.fIntegral
            self.fCutCounts.SetBinContent(ibin, partialInt)
            self.fCutCounts.SetBinError(ibin, math.sqrt(partialIntErr2))
            self.fCutEfficiency.SetBinContent(ibin, eff)
            self.fCutEfficiency.SetBinError(ibin, effErr)

    def CalculateRightCutEfficiency(self):
        partialIntErr2 = 0.
        partialInt = 0.
        for ibin in reversed(xrange(0, self.fHistogram.GetNbinsX() + 2)):
            partialIntErr2 += self.fHistogram.GetBinError(ibin) ** 2
            partialInt += self.fHistogram.GetBinContent(ibin)
            eff = partialInt / self.fIntegral
            effErr = math.sqrt(partialIntErr2) / self.fIntegral
            self.fCutCounts.SetBinContent(ibin, partialInt)
            self.fCutCounts.SetBinError(ibin, math.sqrt(partialIntErr2))
            self.fCutEfficiency.SetBinContent(ibin, eff)
            self.fCutEfficiency.SetBinError(ibin, effErr)

    def CalculateFractions(self, signal, anaName):
        if "Bkg" in anaName:
            self.CalculateSignificance(signal)
            self.CalculateFraction(signal, True)
            self.fSignificance.GetYaxis().SetTitle("S / #sqrt{S+B}")
            self.fFraction.GetYaxis().SetTitle("S / B")
        else:
            self.CalculateFraction(signal, False)
            self.fFraction.GetYaxis().SetTitle("Non-Prompt / Prmpt")

    def CalculateSignificance(self, signal):
        self.fSignificance = self.fHistogram.Clone("{}_CutSignificance".format(self.fName))
        self.fSignificance.Reset()
        self.fMaxSignificance = -1
        for ibin in xrange(0, self.fHistogram.GetNbinsX() + 2):
            totSig = signal.fCutEfficiency.GetBinContent(ibin) * signal.fIntegral
            totBkg = self.fCutEfficiency.GetBinContent(ibin) * self.fIntegral
            totSigErr = signal.fCutEfficiency.GetBinError(ibin) * signal.fIntegral
            totBkgErr = self.fCutEfficiency.GetBinError(ibin) * self.fIntegral
            if totSig == 0: continue
            significance = totSig / math.sqrt(totSig + totBkg)
            if significance > self.fMaxSignificance:
                self.fMaxSignificance = significance
                if self.CalculateCutEfficiency == self.CalculateLeftCutEfficiency:
                    self.fCutValueAtMaxSignificance = self.fHistogram.GetXaxis().GetBinUpEdge(ibin)
                else:
                    self.fCutValueAtMaxSignificance = self.fHistogram.GetXaxis().GetBinLowEdge(ibin)
            significanceErr = math.sqrt(((totSig * totBkgErr) ** 2 + ((2 * totBkg + totSig) * totSigErr) ** 2) / ((totSig + totBkg) ** 3)) / 2
            self.fSignificance.SetBinContent(ibin, significance)
            self.fSignificance.SetBinError(ibin, significanceErr)

    def CalculateFraction(self, signal, rev):
        self.fFraction = self.fHistogram.Clone("{}_CutFraction".format(self.fName))
        self.fFraction.Reset()
        if rev: self.fMaxSigOverBkg = -1
        for ibin in xrange(0, self.fHistogram.GetNbinsX() + 2):
            totSig = signal.fCutEfficiency.GetBinContent(ibin) * signal.fIntegral
            totBkg = self.fCutEfficiency.GetBinContent(ibin) * self.fIntegral
            totSigErr = signal.fCutEfficiency.GetBinError(ibin) * signal.fIntegral
            totBkgErr = self.fCutEfficiency.GetBinError(ibin) * self.fIntegral
            if totBkg == 0 or totSig == 0: continue
            if rev:
                frac = totSig / totBkg
                if frac > self.fMaxSigOverBkg:
                    self.fMaxSigOverBkg = frac
                    if self.CalculateCutEfficiency == self.CalculateLeftCutEfficiency:
                        self.fCutValueAtMaxSigOverBkg = self.fHistogram.GetXaxis().GetBinUpEdge(ibin)
                    else:
                        self.fCutValueAtMaxSigOverBkg = self.fHistogram.GetXaxis().GetBinLowEdge(ibin)
            else:
                frac = totBkg / totSig
            fracErr = frac * math.sqrt((totSigErr / totSig) ** 2 + (totBkgErr / totBkg) ** 2)
            self.fFraction.SetBinContent(ibin, frac)
            self.fFraction.SetBinError(ibin, fracErr)

    @classmethod
    def DCA(cls):
        bins = numpy.linspace(0.0, 0.3, 61, True)
        obj = cls("DCA", "fDCA", "DCA (cm)", bins)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def CosThetaStar(cls):
        bins = numpy.linspace(0.0, 1.0, 51, True)
        obj = cls("CosThetaStar", "fCosThetaStar", "|cos(#theta*)|", bins)
        obj.Fill = obj.FillAbs
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def d0d0(cls):
        bins = numpy.linspace(-500e-6, 500e-6, 251, True)
        obj = cls("d0d0", "fd0d0", "#it{d}_{0,#pi}#it{d}_{0,K} (cm^{2})", bins)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def MaxNormd0(cls):
        bins = numpy.linspace(0.0, 10.0, 21, True)
        obj = cls("MaxNormd0", "fMaxNormd0", "max(|#it{d}_{0,#pi}|, |#it{d}_{0,K}|) / #sigma(#it{d}_{0})", bins)
        obj.Fill = obj.FillAbs
        obj.CalculateCutEfficiency = obj.CalculateLeftCutEfficiency
        return obj

    @classmethod
    def CosPointing(cls):
        bins = numpy.linspace(0.0, 1.0, 51, True)
        obj = cls("CosPointing", "fCosPointing", "cos(#theta_{p})", bins)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateRightCutEfficiency
        return obj

    @classmethod
    def PtK(cls):
        bins = numpy.concatenate((numpy.linspace(0.0, 5.0, 50), numpy.linspace(5.0, 20.0, 30), numpy.linspace(20.0, 40.0, 5, True)))
        obj = cls("PtK", "fPtK", "#it{p}_{T,K}", bins)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateRightCutEfficiency
        return obj

    @classmethod
    def PtPi(cls):
        bins = numpy.concatenate((numpy.linspace(0.0, 5.0, 50), numpy.linspace(5.0, 20.0, 30), numpy.linspace(20.0, 40.0, 5, True)))
        obj = cls("PtPi", "fPtPi", "#it{p}_{T,#pi}", bins)
        obj.Fill = obj.FillStd
        obj.CalculateCutEfficiency = obj.CalculateRightCutEfficiency
        return obj

DMesonJetVariable.fgVariableList = [DMesonJetVariable.DCA, DMesonJetVariable.CosThetaStar, DMesonJetVariable.d0d0, DMesonJetVariable.CosPointing]

class DMesonJetTopoContainer:
    def __init__(self, jet_pt_bins, d_pt_bins, sigma_mass, jet_type, jet_radius):
        self.fJetPtBins = jet_pt_bins
        self.fDPtBins = d_pt_bins
        self.fVariables = collections.OrderedDict()
        self.fJetType = jet_type
        self.fJetRadius = jet_radius
        self.fPDGMass = 1.86484
        self.fSigmaMass = sigma_mass
        for jetPtBin in self.fJetPtBins:
            self.fVariables[jetPtBin] = collections.OrderedDict()
            for dPtBin in self.fDPtBins[jetPtBin]:
                self.fVariables[jetPtBin][dPtBin] = self.GenerateVariables()

    def GetMassRange(self, jetPtBin, dPtBin):
        sigma = self.fSigmaMass[jetPtBin][dPtBin]
        return self.fPDGMass - 3 * sigma, self.fPDGMass + 3 * sigma

    def AddVariable(self, variables, v):
        variables[v.fName] = v

    def GenerateVariables(self):
        variables = collections.OrderedDict()

        for var in DMesonJetVariable.fgVariableList:
            self.AddVariable(variables, var())

        return variables

    def Fill(self, event, eventWeight):
        dmeson = event.DmesonJet

        if self.fJetType or self.fJetRadius:
            jetName = "Jet_AKT{0}{1}_pt_scheme".format(self.fJetType, self.fJetRadius)
            jet = getattr(event, jetName)
        else:
            print("Error no jet found!")
            exit(1)

        for (jetPtMin, jetPtMax) in self.fJetPtBins:
            if jet.fPt < jetPtMin or jet.fPt >= jetPtMax: continue
            for (dPtMin, dPtMax) in self.fDPtBins[(jetPtMin, jetPtMax)]:
                if dmeson.fPt < dPtMin or dmeson.fPt >= dPtMax: continue
                (minMass, maxMass) = self.GetMassRange((jetPtMin, jetPtMax), (dPtMin, dPtMax))
                if dmeson.fInvMass < minMass or dmeson.fInvMass > maxMass: continue
                for var in self.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)].itervalues():
                    var.Fill(dmeson, eventWeight)

class DMesonJetTopoAnalysis:
    def __init__(self, name, title, trigger, dmeson, jet_pt_bins, d_pt_bins, sigma_mass, projector, norm):
        self.fName = name
        self.fTitle = title
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fProjector = projector
        self.fJetPtBins = jet_pt_bins
        self.fDPtBins = d_pt_bins
        self.fSigmaMass = sigma_mass
        self.fNormalization = norm

    def SaveRootFile(self, file):
        rlist = ROOT.TList()
        rlist.SetName(self.fName)

        for (jetPtMin, jetPtMax), jetPtBin in self.fVariables.iteritems():
            for (dPtMin, dPtMax), dPtBin in jetPtBin.iteritems():
                hListName = "Histograms_JetPt{}_{}_DPt_{}_{}".format(jetPtMin, jetPtMax, dPtMin, dPtMax)
                hList = ROOT.TList()
                hList.SetName(hListName)
                rlist.Add(hList)

                eListName = "CutEfficiency_JetPt{}_{}_DPt_{}_{}".format(jetPtMin, jetPtMax, dPtMin, dPtMax)
                eList = ROOT.TList()
                eList.SetName(eListName)
                rlist.Add(eList)

                cListName = "CutCounts_JetPt{}_{}_DPt_{}_{}".format(jetPtMin, jetPtMax, dPtMin, dPtMax)
                cList = ROOT.TList()
                cList.SetName(cListName)
                rlist.Add(cList)

                fList = None
                sList = None

                if "Bkg" in self.fName:
                    fListName = "SigOverBkg_JetPt{}_{}_DPt_{}_{}".format(jetPtMin, jetPtMax, dPtMin, dPtMax)
                    fList = ROOT.TList()
                    fList.SetName(fListName)
                    rlist.Add(fList)

                    sListName = "Significance_JetPt{}_{}_DPt_{}_{}".format(jetPtMin, jetPtMax, dPtMin, dPtMax)
                    sList = ROOT.TList()
                    sList.SetName(sListName)
                    rlist.Add(sList)

                if "NonPrompt" in self.fName:
                    fListName = "NonPromptFraction_JetPt{}_{}_DPt_{}_{}".format(jetPtMin, jetPtMax, dPtMin, dPtMax)
                    fList = ROOT.TList()
                    fList.SetName(fListName)
                    rlist.Add(fList)

                for v in dPtBin.itervalues():
                    hList.Add(v.fHistogram)
                    eList.Add(v.fCutEfficiency)
                    cList.Add(v.fCutCounts)
                    if fList and v.fFraction: fList.Add(v.fFraction)
                    if sList and v.fSignificance: fList.Add(v.fSignificance)

        if rlist.GetEntries() > 0:
            file.cd()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)

    def DoProjections(self):
        data_container = DMesonJetTopoContainer(self.fJetPtBins, self.fDPtBins, self.fSigmaMass, "Charged", "R040")
        self.fProjector.StartProjection(self.fTrigger, self.fDMeson, None, data_container, self.fNormalization)
        self.fVariables = data_container.fVariables

    def Analyze(self):
        for (jetPtMin, jetPtMax), jetPtBin in self.fVariables.iteritems():
            for (dPtMin, dPtMax), dPtBin in jetPtBin.iteritems():
                for vname, v in dPtBin.iteritems():
                    v.Analyze()

class DMesonJetTopoAnalysisManager:
    def __init__(self, dmeson, jet_pt_bins):
        self.fDMeson = dmeson
        self.fAnalysisList = collections.OrderedDict()
        self.fDPtBins = dict()
        self.fSigmaMass = dict()
        self.fJetPtBins = []
        for jet_pt_bin in jet_pt_bins:
            self.fJetPtBins.append((jet_pt_bin["min"], jet_pt_bin["max"]))
            d_pt_bin_limits = jet_pt_bin["d_pt_bins"]
            self.fDPtBins[(jet_pt_bin["min"], jet_pt_bin["max"])] = []
            self.fSigmaMass[(jet_pt_bin["min"], jet_pt_bin["max"])] = dict()
            for i, (minpt, maxpt) in enumerate(zip(d_pt_bin_limits[:-1], d_pt_bin_limits[1:])):
                self.fDPtBins[(jet_pt_bin["min"], jet_pt_bin["max"])].append((minpt, maxpt))
                self.fSigmaMass[(jet_pt_bin["min"], jet_pt_bin["max"])][(minpt, maxpt)] = jet_pt_bin["sigma_mass"][i]

    def AddAnalysis(self, name, title, trigger, projector, dmesonSuffix, norm):
        self.fAnalysisList[name] = DMesonJetTopoAnalysis(name, title, trigger, "{}_{}".format(self.fDMeson, dmesonSuffix), self.fJetPtBins, self.fDPtBins, self.fSigmaMass, projector, norm)

    def Analyze(self):
        for ana in self.fAnalysisList.itervalues():
            ana.Analyze()

    def DoProjections(self):
        for ana in self.fAnalysisList.itervalues():
            ana.DoProjections()

    def StartAnalysis(self):
        self.DoProjections()
        self.Analyze()
        self.CalculateFractions()
        self.Plot()

    def SaveRootFile(self, path):
        fname = "{}/TopoAnalysis.root".format(path)
        file = ROOT.TFile(fname, "recreate")
        if not file or file.IsZombie():
            print("Could not open file '{}'".format(fname))
            return
        for ana in self.fAnalysisList.itervalues():
            ana.SaveRootFile(file)
        print("Results stored in '{}'".format(fname))

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            if c: c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))

    def CalculateFractions(self):
        for (jetPtMin, jetPtMax) in self.fJetPtBins:
            for (dPtMin, dPtMax) in self.fDPtBins[(jetPtMin, jetPtMax)]:
                for variableName in [varFunct().fName for varFunct in DMesonJetVariable.fgVariableList]:
                    if not "Signal" in self.fAnalysisList.values()[0].fName and not "Sig" in self.fAnalysisList.values()[0].fName:
                        print("**ERROR** The first analysis in the list must be the signal. No fractions calculated (see CalculateFractions).")
                        return
                    signalAna = self.fAnalysisList.values()[0]
                    for ana in self.fAnalysisList.values()[1:]:
                        ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].CalculateFractions(signalAna.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName], ana.fName)

    def Plot(self):
        self.fCanvases = []
        self.fCompareObjects = []
        self.fKeepObjects = []
        summary = ""
        for (jetPtMin, jetPtMax) in self.fJetPtBins:
            for (dPtMin, dPtMax) in self.fDPtBins[(jetPtMin, jetPtMax)]:
                for variableName in [varFunct().fName for varFunct in DMesonJetVariable.fgVariableList]:
                    summary += "******\nJet Pt bin [{}. {}]\nD Pt bin [{}. {}]\nVariable name {}\n******\n".format(jetPtMin, jetPtMax, dPtMin, dPtMax, variableName)
                    hVariable = []
                    hEfficiency = []
                    hCounts = []
                    cutValueMaxSignificance = None
                    cutValueMaxSigOverBkg = None
                    for ana in self.fAnalysisList.itervalues():
                        if not ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fCutValueAtMaxSignificance is None:
                            cutValueMaxSignificance = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fCutValueAtMaxSignificance
                            significance = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fMaxSignificance
                            summary += "Name: {}\nCut value: {}\nSignificance: {}\n\n".format(ana.fTitle, cutValueMaxSignificance, significance)

                        if not ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fCutValueAtMaxSigOverBkg is None:
                            cutValueMaxSigOverBkg = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fCutValueAtMaxSigOverBkg
                            sigOverBkg = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fMaxSigOverBkg
                            summary += "Name: {}\nCut value: {}\nS/B: {}\n\n".format(ana.fTitle, cutValueMaxSigOverBkg, sigOverBkg)

                        h = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fHistogram
                        hVar = h.Clone("{}_copy".format(h.GetName()))
                        hVar.SetTitle(ana.fTitle)
                        hVariable.append(hVar)

                        h = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fCutEfficiency
                        hEff = h.Clone("{}_copy".format(h.GetName()))
                        hEff.SetTitle(ana.fTitle)
                        hEfficiency.append(hEff)

                        h = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fCutCounts
                        hCts = h.Clone("{}_copy".format(h.GetName()))
                        hCts.SetTitle(ana.fTitle)
                        hCounts.append(hCts)

                        hFraction = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fFraction
                        if hFraction:
                            cname = "{}_{}_JetPt{}_{}_DPt{}_{}".format(hFraction.GetName(), ana.fName, jetPtMin, jetPtMax, dPtMin, dPtMax)
                            c = ROOT.TCanvas(cname, cname)
                            c.SetGridx()
                            c.SetGridy()
                            c.cd()
                            hFraction.SetMarkerStyle(ROOT.kFullCircle)
                            hFraction.SetMarkerSize(0.9)
                            hFraction.SetMarkerColor(ROOT.kBlue + 2)
                            hFraction.SetLineColor(ROOT.kBlue + 2)
                            hFraction.Draw()
                            c.Update()
                            if not cutValueMaxSignificance is None:
                                line = ROOT.TLine(cutValueMaxSignificance, c.GetUymin(), cutValueMaxSignificance, c.GetUymax())
                                line.SetLineColor(ROOT.kMagenta + 2)
                                line.SetLineStyle(2)
                                line.SetLineWidth(2)
                                line.Draw()
                                self.fKeepObjects.append(line)
                            if not cutValueMaxSigOverBkg is None:
                                line = ROOT.TLine(cutValueMaxSigOverBkg, c.GetUymin(), cutValueMaxSigOverBkg, c.GetUymax())
                                line.SetLineColor(ROOT.kGreen + 2)
                                line.SetLineStyle(1)
                                line.SetLineWidth(2)
                                line.Draw()
                                self.fKeepObjects.append(line)
                            self.fCanvases.append(c)

                        hSignificance = ana.fVariables[(jetPtMin, jetPtMax)][(dPtMin, dPtMax)][variableName].fSignificance
                        if hSignificance:
                            cname = "{}_{}_JetPt{}_{}_DPt{}_{}".format(hSignificance.GetName(), ana.fName, jetPtMin, jetPtMax, dPtMin, dPtMax)
                            c = ROOT.TCanvas(cname, cname)
                            c.SetGridx()
                            c.SetGridy()
                            c.cd()
                            hSignificance.SetMarkerStyle(ROOT.kFullCircle)
                            hSignificance.SetMarkerSize(0.9)
                            hSignificance.SetMarkerColor(ROOT.kBlue + 2)
                            hSignificance.SetLineColor(ROOT.kBlue + 2)
                            hSignificance.Draw()
                            c.Update()
                            if not cutValueMaxSignificance is None:
                                line = ROOT.TLine(cutValueMaxSignificance, c.GetUymin(), cutValueMaxSignificance, c.GetUymax())
                                line.SetLineColor(ROOT.kMagenta + 2)
                                line.SetLineStyle(2)
                                line.SetLineWidth(2)
                                line.Draw()
                                self.fKeepObjects.append(line)
                            if not cutValueMaxSigOverBkg is None:
                                line = ROOT.TLine(cutValueMaxSigOverBkg, c.GetUymin(), cutValueMaxSigOverBkg, c.GetUymax())
                                line.SetLineColor(ROOT.kGreen + 2)
                                line.SetLineStyle(1)
                                line.SetLineWidth(2)
                                line.Draw()
                                self.fKeepObjects.append(line)
                            self.fCanvases.append(c)

                    cname = "{}_JetPt{}_{}_DPt{}_{}".format(hEfficiency[0].GetName(), jetPtMin, jetPtMax, dPtMin, dPtMax)
                    comp = DMesonJetCompare.DMesonJetCompare(cname)
                    comp.fDoSpectraPlot = "lineary"
                    comp.fDoRatioPlot = False
                    comp.CompareSpectra(hEfficiency[0], hEfficiency[1:])
                    if not cutValueMaxSignificance is None:
                        line = ROOT.TLine(cutValueMaxSignificance, comp.fMainHistogram.GetMinimum(), cutValueMaxSignificance, comp.fMainHistogram.GetMaximum())
                        line.SetLineColor(ROOT.kMagenta + 2)
                        line.SetLineStyle(2)
                        line.SetLineWidth(2)
                        line.Draw()
                        self.fKeepObjects.append(line)
                    if not cutValueMaxSigOverBkg is None:
                        line = ROOT.TLine(cutValueMaxSigOverBkg, comp.fMainHistogram.GetMinimum(), cutValueMaxSigOverBkg, comp.fMainHistogram.GetMaximum())
                        line.SetLineColor(ROOT.kGreen + 2)
                        line.SetLineStyle(1)
                        line.SetLineWidth(2)
                        line.Draw()
                        self.fKeepObjects.append(line)
                    comp.fCanvasSpectra.SetGridx()
                    comp.fCanvasSpectra.SetGridy()
                    self.fCompareObjects.append(comp)
                    self.fCanvases.append(comp.fCanvasSpectra)

                    cname = "{}_JetPt{}_{}_DPt{}_{}".format(hVariable[0].GetName(), jetPtMin, jetPtMax, dPtMin, dPtMax)
                    comp = DMesonJetCompare.DMesonJetCompare(cname)
                    comp.fDoSpectraPlot = "logy"
                    comp.fDoRatioPlot = False
                    comp.CompareSpectra(hVariable[0], hVariable[1:])
                    if not cutValueMaxSignificance is None:
                        line = ROOT.TLine(cutValueMaxSignificance, comp.fMainHistogram.GetMinimum(), cutValueMaxSignificance, comp.fMainHistogram.GetMaximum())
                        line.SetLineColor(ROOT.kMagenta + 2)
                        line.SetLineStyle(2)
                        line.SetLineWidth(2)
                        line.Draw()
                        self.fKeepObjects.append(line)
                    if not cutValueMaxSigOverBkg is None:
                        line = ROOT.TLine(cutValueMaxSigOverBkg, comp.fMainHistogram.GetMinimum(), cutValueMaxSigOverBkg, comp.fMainHistogram.GetMaximum())
                        line.SetLineColor(ROOT.kGreen + 2)
                        line.SetLineStyle(1)
                        line.SetLineWidth(2)
                        line.Draw()
                        self.fKeepObjects.append(line)
                    comp.fCanvasSpectra.SetGridx()
                    comp.fCanvasSpectra.SetGridy()
                    self.fCompareObjects.append(comp)
                    self.fCanvases.append(comp.fCanvasSpectra)

                    cname = "{}_JetPt{}_{}_DPt{}_{}".format(hCounts[0].GetName(), jetPtMin, jetPtMax, dPtMin, dPtMax)
                    comp = DMesonJetCompare.DMesonJetCompare(cname)
                    comp.fDoSpectraPlot = "logy"
                    comp.fDoRatioPlot = False
                    comp.CompareSpectra(hCounts[0], hCounts[1:])
                    if not cutValueMaxSignificance is None:
                        line = ROOT.TLine(cutValueMaxSignificance, comp.fMainHistogram.GetMinimum(), cutValueMaxSignificance, comp.fMainHistogram.GetMaximum())
                        line.SetLineColor(ROOT.kMagenta + 2)
                        line.SetLineStyle(2)
                        line.SetLineWidth(2)
                        line.Draw()
                        self.fKeepObjects.append(line)
                    if not cutValueMaxSigOverBkg is None:
                        line = ROOT.TLine(cutValueMaxSigOverBkg, comp.fMainHistogram.GetMinimum(), cutValueMaxSigOverBkg, comp.fMainHistogram.GetMaximum())
                        line.SetLineColor(ROOT.kGreen + 2)
                        line.SetLineStyle(1)
                        line.SetLineWidth(2)
                        line.Draw()
                        self.fKeepObjects.append(line)
                    comp.fCanvasSpectra.SetGridx()
                    comp.fCanvasSpectra.SetGridy()
                    self.fCompareObjects.append(comp)
                    self.fCanvases.append(comp.fCanvasSpectra)

        print(summary)

