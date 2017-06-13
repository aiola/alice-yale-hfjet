#!/usr/local/bin/python
# python program to perform a D meson jet analysis

import math
import array
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

globalList = []

class DMesonJetAnalysisEngine:
    def __init__(self, collision, trigger, dmeson, binSet, nMassBins, minMass, maxMass, jets, projector):
        self.fCollision = collision
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fBinMultiSets = dict()
        self.fNMassBins = nMassBins
        self.fMinMass = minMass
        self.fMaxMass = maxMass
        self.fCanvases = []

        for jetDef in self.fJetDefinitions:
            if "active" in jetDef and not jetDef["active"]:
                continue
            binset_copy = copy.deepcopy(binSet)
            binset_copy.Initialize(self.fDMeson, jetDef["type"], jetDef["radius"], jetDef["title"], self.fProjector.fInputPath)
            self.fBinMultiSets[jetDef["type"], jetDef["radius"]] = binset_copy

        binset_copy = copy.deepcopy(binSet)
        binset_copy.Initialize(self.fDMeson, None, None, None, self.fProjector.fInputPath)
        self.fBinMultiSets[None, None] = binset_copy

    def CompareSpectra(self):
        for binMultiSet in self.fBinMultiSets.itervalues():
            self.CompareSpectraForAxis("jet_pt", binMultiSet)
            self.CompareSpectraForAxis("d_pt", binMultiSet)
            self.CompareSpectraForAxis("d_z", binMultiSet)

    def CompareSpectraForAxis(self, axisName, binMultiSet):
        cname = '_'.join(obj for obj in [self.fTrigger, self.fDMeson, binMultiSet.fJetType, binMultiSet.fJetRadius, axisName, "SpectraComparison"] if obj)
        spectraToCompare = []
        spectra = binMultiSet.FindAllSpectra()
        for s in spectra:
            if not s.fCompare:
                continue
            if not s.fNormHistogram:
                continue
            if len(s.fAxis) != 1:
                continue
            if axisName != s.fAxis[0].fName:
                continue
            h = s.fNormHistogram.Clone("{0}_{1}_copy".format(s.fNormHistogram.GetName(), cname))
            if s.fTitle:
                h.SetTitle(s.fTitle)
            globalList.append(h)
            spectraToCompare.append(h)
        if len(spectraToCompare) < 2:
            return

        comp = DMesonJetCompare.DMesonJetCompare(cname)
        comp.fOptRatio = "hist"
        results = comp.CompareSpectra(spectraToCompare[0], spectraToCompare[1:])
        for obj in results:
            if obj and isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)
                obj.cd()
                pave = ROOT.TPaveText(0.12, 0.70, 0.40, 0.85, "NB NDC")
                pave.SetTextAlign(11)
                pave.SetFillStyle(0)
                pave.SetBorderSize(0)
                pave.SetTextFont(43)
                pave.SetTextSize(15)
                pave.AddText(self.fCollision)
                if s.fJetTitle: pave.AddText(s.fJetTitle)
                pave.AddText(DMesonJetUtils.ConvertDMesonName(self.fDMeson))
                pave.Draw()
                globalList.append(pave)
            globalList.append(obj)

    def SaveRootFile(self, file):
        rlist = ROOT.TList()
        if self.fTrigger:
            rlist.SetName("{}_{}".format(self.fTrigger, self.fDMeson))
        else:
            rlist.SetName(self.fDMeson)

        for (jtype, jradius), binMultiSet in self.fBinMultiSets.iteritems():
            if jtype or jradius:
                jetName = "_".join(obj for obj in [jtype, jradius] if obj)
                jlist = ROOT.TList()
                jlist.SetName(jetName)
            else:
                jetName = None
                jlist = rlist
            for invmasslist in binMultiSet.GenerateInvMassRootLists():
                if invmasslist.GetEntries() > 0:
                    jlist.Add(invmasslist)
            spectra = binMultiSet.FindAllSpectra()
            for s in spectra:
                slist = s.GenerateRootList()
                if slist.GetEntries() > 0:
                    jlist.Add(slist)
            if jlist is not rlist and jlist.GetEntries() > 0:
                rlist.Add(jlist)

        if rlist.GetEntries() > 0:
            file.cd()
            rlist.Write(rlist.GetName(), ROOT.TObject.kSingleKey)

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))

    def CreateMassFitter(self, name, invMassHist):
        if "D0" in self.fDMeson:
            minMass = self.fMinMass
            maxMass = self.fMaxMass
            minFitRange = self.fMinMass
            maxFitRange = self.fMaxMass
            DMesonType = ROOT.MassFitter.kDzeroKpi

            # some defualt params
            pdgMass = ROOT.TDatabasePDG.Instance().GetParticle(421).Mass()
            massWidth = 0.012

            minMass = invMassHist.GetXaxis().GetBinLowEdge(1)
            maxMass = invMassHist.GetXaxis().GetBinUpEdge(invMassHist.GetXaxis().GetNbins())
            totIntegral = invMassHist.Integral(1, invMassHist.GetXaxis().GetNbins(), "width")
            integral3sigma = invMassHist.Integral(invMassHist.GetXaxis().FindBin(pdgMass - 3 * massWidth), invMassHist.GetXaxis().FindBin(pdgMass + 3 * massWidth), "width")

            if invMassHist.GetBinContent(invMassHist.GetXaxis().GetNbins()) > 0 and invMassHist.GetBinContent(1) > 0:
                expoParBkg1 = math.log(invMassHist.GetBinContent(invMassHist.GetXaxis().GetNbins()) / invMassHist.GetBinContent(1)) / (maxMass - minMass)
            else:
                expoParBkg1 = -1.0
            if expoParBkg1 == 0: expoParBkg1 = -1.0
            expoParBkg0 = totIntegral / (math.exp(expoParBkg1 * minMass) - math.exp(expoParBkg1 * maxMass)) * (-expoParBkg1)

            sig = integral3sigma - (math.exp(expoParBkg1 * (pdgMass - 3 * massWidth)) - math.exp(expoParBkg1 * (pdgMass + 3 * massWidth))) / (-expoParBkg1) * expoParBkg0
            GaussConst = sig

            fitter = ROOT.MassFitter(name, ROOT.MassFitter.kDzeroKpi, minMass, maxMass)
            fitter.SetFitRange(minMass, maxMass)
            fitter.SetHistogram(invMassHist)
            fitter.GetFitFunction().SetParameter(0, expoParBkg0)
            fitter.GetFitFunction().SetParameter(1, expoParBkg1)
            fitter.GetFitFunction().SetParameter(2, GaussConst)
            fitter.GetFitFunction().SetParameter(3, pdgMass)  # start fitting using PDG mass
            fitter.GetFitFunction().SetParLimits(3, pdgMass * 0.95, pdgMass * 1.05)  # limiting mass parameter +/- 5% of PDG value
            fitter.GetFitFunction().SetParameter(4, massWidth)
            fitter.GetFitFunction().SetParLimits(4, 0, 1)  # limiting width to being positive

        elif self.fDMeson == "DStar":
            print("Not implemented for DStar!")

        if "SignalOnly" in self.fDMeson:
            fitter.DisableBkg()
        elif "BackgroundOnly" in self.fDMeson:
            fitter.DisableSig()
        globalList.append(fitter)

        return fitter

    def DoProjections(self):
        self.fProjector.StartProjection(self.fTrigger, self.fDMeson,
                                        self.fBinMultiSets, self.fNMassBins, self.fMinMass, self.fMaxMass)

        self.fEvents = self.fProjector.fTotalEvents
        self.fIsWeighted = not (self.fProjector.fWeight == 1)

    def Start(self, ana):
        self.fEngines = ana.fAnalysisEngine
        if not "MCTruth" in self.fDMeson and not "WrongPID" in self.fDMeson:
            self.FitInvMassPlots()
        if not "BackgroundOnly" in self.fDMeson and not "WrongPID" in self.fDMeson:
            print("Skipping spectra generation for {0}".format(self.fDMeson))
            self.GenerateSpectra()

        self.PlotInvMassPlots()
        if not "BackgroundOnly" in self.fDMeson and not "WrongPID" in self.fDMeson:
            self.PlotSpectra()

    def PlotSpectra(self):
        for binMultiSet in self.fBinMultiSets.itervalues():
            spectra = binMultiSet.FindAllSpectra()
            for s in spectra:
                if len(s.fAxis) == 1:
                    self.PlotSpectrum1D(s)
                elif len(s.fAxis) == 2:
                    self.PlotSpectrum2D(s)
                elif len(s.fAxis) == 3:
                    self.PlotSpectrum3D(s)
                else:
                    print("Not able to plot spectra with dim > 3!")

    def PlotSpectrum1D(self, s):
        # Spectrum
        if self.fTrigger:
            cname = "{0}_{1}_canvas".format(self.fTrigger, s.fNormHistogram.GetName())
        else:
            cname = "{0}_canvas".format(s.fNormHistogram.GetName())
        c = ROOT.TCanvas(cname, cname)
        c.SetLogy()
        self.fCanvases.append(c)
        c.cd()
        h = s.fNormHistogram.DrawCopy()

        h.SetMarkerColor(ROOT.kBlue + 2)
        h.SetMarkerStyle(ROOT.kFullCircle)
        h.SetMarkerSize(0.9)
        h.SetLineColor(ROOT.kBlue + 2)

        pave = ROOT.TPaveText(0.60, 0.88, 0.9, 0.55, "NB NDC")
        pave.SetFillStyle(0)
        pave.SetBorderSize(0)
        pave.SetTextFont(43)
        pave.SetTextSize(15)
        pave.AddText(self.fCollision)
        if s.fJetTitle: pave.AddText(s.fJetTitle)
        pave.AddText(DMesonJetUtils.ConvertDMesonName(self.fDMeson))
        pave.AddText(s.fTitle)
        pave.Draw()

        globalList.append(c)
        globalList.append(h)
        globalList.append(pave)

        # Uncertainty
        if self.fTrigger:
            cname = "{0}_{1}_canvas".format(self.fTrigger, s.fUncertainty.GetName())
        else:
            cname = "{0}_canvas".format(s.fUncertainty.GetName())
        c = ROOT.TCanvas(cname, cname)
        self.fCanvases.append(c)
        c.cd()
        h = s.fUncertainty.DrawCopy("hist")

        h.SetLineColor(ROOT.kBlue + 2)
        h.SetLineWidth(2)
        h.SetFillColorAlpha(ROOT.kBlue + 2, 0.25)
        h.SetFillStyle(1001)

        pave = ROOT.TPaveText(0.12, 0.88, 0.4, 0.55, "NB NDC")
        pave.SetFillStyle(0)
        pave.SetBorderSize(0)
        pave.SetTextFont(43)
        pave.SetTextSize(15)
        pave.AddText(self.fCollision)
        if s.fJetTitle: pave.AddText(s.fJetTitle)
        pave.AddText(DMesonJetUtils.ConvertDMesonName(self.fDMeson))
        pave.AddText(s.fTitle)
        pave.Draw()

        globalList.append(c)
        globalList.append(h)
        globalList.append(pave)

        # Mass
        if s.fMass:
            if self.fTrigger:
                cname = "{0}_{1}_canvas".format(self.fTrigger, s.fMass.GetName())
            else:
                cname = "{0}_canvas".format(s.fMass.GetName())
            c = ROOT.TCanvas(cname, cname)
            self.fCanvases.append(c)
            c.cd()
            h = s.fMass.DrawCopy()

            h.SetMarkerColor(ROOT.kBlue + 2)
            h.SetMarkerStyle(ROOT.kFullCircle)
            h.SetMarkerSize(0.9)
            h.SetLineColor(ROOT.kBlue + 2)
            h.GetYaxis().SetRangeUser(1.84, 1.91)

            pave = ROOT.TPaveText(0.12, 0.88, 0.8, 0.68, "NB NDC")
            pave.SetFillStyle(0)
            pave.SetBorderSize(0)
            pave.SetTextFont(43)
            pave.SetTextSize(15)
            pave.AddText(self.fCollision)
            if s.fJetTitle: pave.AddText(s.fJetTitle)
            pave.AddText(DMesonJetUtils.ConvertDMesonName(self.fDMeson))
            pave.AddText(s.fTitle)
            pave.Draw()

            line = ROOT.TLine(h.GetXaxis().GetBinLowEdge(1), 1.86484, h.GetXaxis().GetBinUpEdge(h.GetXaxis().GetNbins()), 1.86484)
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()

            globalList.append(c)
            globalList.append(h)
            globalList.append(pave)
            globalList.append(line)

        # Mass width
        if s.fMassWidth:
            if self.fTrigger:
                cname = "{0}_{1}_canvas".format(self.fTrigger, s.fMassWidth.GetName())
            else:
                cname = "{0}_canvas".format(s.fMassWidth.GetName())
            c = ROOT.TCanvas(cname, cname)
            self.fCanvases.append(c)
            c.cd()
            h = s.fMassWidth.DrawCopy()

            h.SetMarkerColor(ROOT.kBlue + 2)
            h.SetMarkerStyle(ROOT.kFullCircle)
            h.SetMarkerSize(0.9)
            h.SetLineColor(ROOT.kBlue + 2)

            pave = ROOT.TPaveText(0.12, 0.88, 0.8, 0.68, "NB NDC")
            pave.SetFillStyle(0)
            pave.SetBorderSize(0)
            pave.SetTextFont(43)
            pave.SetTextSize(15)
            pave.AddText(self.fCollision)
            if s.fJetTitle: pave.AddText(s.fJetTitle)
            pave.AddText(DMesonJetUtils.ConvertDMesonName(self.fDMeson))
            pave.AddText(s.fTitle)
            pave.Draw()

            globalList.append(c)
            globalList.append(h)
            globalList.append(pave)
            globalList.append(line)

        if s.fLikeSignTotalHistogram:
            if self.fTrigger:
                cname = "{0}_{1}_TotalBkgVsSig".format(self.fTrigger, s.fName)
            else:
                cname = "{0}_TotalBkgVsSig".format(s.fName)
            self.PlotBackgroundVsSignalSpectra(cname, s.fUnlikeSignTotalHistogram, "Unlike Sign", s.fLikeSignTotalHistogram, "Like Sign")
        if s.fUnlikeSignHistograms and s.fLikeSignHistograms:
            if self.fTrigger:
                cname = "{0}_{1}_BkgVsSig".format(self.fTrigger, s.fName)
            else:
                cname = "{0}_BkgVsSig".format(s.fName)
            self.PlotMultiCanvasBkgVsSigSpectra(cname, s.fUnlikeSignHistograms, "Unlike Sign", s.fLikeSignHistograms, "Like Sign")

        if s.fSideBandWindowTotalHistogram:
            if self.fTrigger:
                cname = "{0}_{1}_TotalBkgVsSig".format(self.fTrigger, s.fName)
            else:
                cname = "{0}_TotalBkgVsSig".format(s.fName)
            self.PlotBackgroundVsSignalSpectra(cname, s.fSignalWindowTotalHistogram, "Sig. Window", s.fSideBandWindowTotalHistogram, "SB Window")
        if s.fSignalHistograms and s.fSideBandHistograms:
            if self.fTrigger:
                cname = "{0}_{1}_BkgVsSig".format(self.fTrigger, s.fName)
            else:
                cname = "{0}_BkgVsSig".format(s.fName)
            self.PlotMultiCanvasBkgVsSigSpectra(cname, s.fSignalHistograms, "Sig. Window", s.fSideBandHistograms, "SB Window")

    def PlotMultiCanvasBkgVsSigSpectra(self, cname, sigHistograms, sigTitle, bkgHistograms, bkgTitle):
        ncanvases = len(sigHistograms)
        c = DMesonJetUtils.GenerateMultiCanvas(cname, ncanvases)
        self.fCanvases.append(c)
        globalList.append(c)
        for i, (sig, bkg) in enumerate(zip(sigHistograms, bkgHistograms)):
            pad = c.cd(i + 1)
            pad.SetLeftMargin(0.12)
            pad.SetRightMargin(0.05)
            pad.SetTopMargin(0.08)
            pad.SetBottomMargin(0.13)
            (hSig, hBkg, hSub) = self.PlotBackgroundVsSignalSpectra(None, sig, None, bkg, None)
            hSig.GetXaxis().SetTitleFont(43)
            hSig.GetXaxis().SetTitleOffset(2.3)
            hSig.GetXaxis().SetTitleSize(19)
            hSig.GetXaxis().SetLabelFont(43)
            hSig.GetXaxis().SetLabelOffset(0.009)
            hSig.GetXaxis().SetLabelSize(18)
            hSig.GetYaxis().SetTitleFont(43)
            hSig.GetYaxis().SetTitleOffset(2.3)
            hSig.GetYaxis().SetTitleSize(19)
            hSig.GetYaxis().SetLabelFont(43)
            hSig.GetYaxis().SetLabelOffset(0.009)
            hSig.GetYaxis().SetLabelSize(18)
            htitle = ROOT.TPaveText(0.12, 0.99, 0.95, 0.93, "NB NDC")
            htitle.SetBorderSize(0)
            htitle.SetFillStyle(0)
            htitle.SetTextFont(43)
            htitle.SetTextSize(18)
            htitle.AddText(sig.GetTitle())
            htitle.Draw()
            globalList.append(htitle)
        c.cd(1)
        leg = ROOT.TLegend(0.27, 0.72, 0.72, 0.87, "", "NB NDC")
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(16)
        leg.SetMargin(0.2)
        leg.AddEntry(hSig, sigTitle, "pe")
        leg.AddEntry(hBkg, bkgTitle, "pe")
        leg.AddEntry(hSub, "{0} - {1}".format(sigTitle, bkgTitle), "pe")
        leg.Draw()
        globalList.append(leg)

    def PlotBackgroundVsSignalSpectra(self, cname, hSigOrig, hTitleSig, hBkgOrig, hTitleBkg):
        if cname:
            c = ROOT.TCanvas(cname, cname, 700, 700)
            c.SetTopMargin(0.05)
            self.fCanvases.append(c)
            globalList.append(c)

        hSig = hSigOrig.DrawCopy()
        hSig.SetMarkerColor(ROOT.kBlue + 2)
        hSig.SetMarkerStyle(ROOT.kOpenCircle)
        hSig.SetMarkerSize(0.9)
        hSig.SetLineColor(ROOT.kBlue + 2)
        hSig.GetYaxis().SetTitleFont(43)
        hSig.GetYaxis().SetTitleSize(22)
        hSig.GetYaxis().SetTitleOffset(1.1)
        hSig.GetYaxis().SetTitle("yield")
        hSig.GetYaxis().SetLabelFont(43)
        hSig.GetYaxis().SetLabelSize(20)
        hSig.GetXaxis().SetTitleFont(43)
        hSig.GetXaxis().SetTitleSize(22)
        hSig.GetXaxis().SetLabelFont(43)
        hSig.GetXaxis().SetLabelSize(20)
        hSig.SetMaximum(hSig.GetMaximum() * 1.3)
        hSig.SetMinimum(0)

        hBkg = hBkgOrig.DrawCopy("same")
        hBkg.SetMarkerColor(ROOT.kRed + 2)
        hBkg.SetMarkerStyle(ROOT.kOpenSquare)
        hBkg.SetMarkerSize(0.9)
        hBkg.SetLineColor(ROOT.kRed + 2)

        hSub = hSigOrig.DrawCopy("same")
        hSub.Add(hBkgOrig, -1)
        hSub.SetMarkerColor(ROOT.kGreen + 2)
        hSub.SetMarkerStyle(ROOT.kStar)
        hSub.SetMarkerSize(0.9)
        hSub.SetLineColor(ROOT.kGreen + 2)

        if hTitleSig and hTitleBkg:
            leg = ROOT.TLegend(0.42, 0.84, 0.87, 0.93, "", "NB NDC")
            leg.SetBorderSize(0)
            leg.SetFillStyle(0)
            leg.SetTextFont(43)
            leg.SetTextSize(20)
            leg.SetMargin(0.2)
            leg.AddEntry(hSig, hTitleSig, "pe")
            leg.AddEntry(hBkg, hTitleBkg, "pe")
            leg.AddEntry(hSub, "{0} - {1}".format(hTitleSig, hTitleBkg), "pe")
            leg.Draw()
            globalList.append(leg)

        globalList.append(hSig)
        globalList.append(hBkg)
        globalList.append(hSub)

        return hSig, hBkg, hSub

    def PlotSpectrum2D(self, s):
        if self.fTrigger:
            cname = "{0}_{1}_canvas".format(self.fTrigger, s.fNormHistogram.GetName())
        else:
            cname = "{0}_canvas".format(s.fNormHistogram.GetName())

        c = ROOT.TCanvas(cname, cname)
        c.SetRightMargin(0.18)
        c.SetLogz()
        self.fCanvases.append(c)
        c.cd()
        h = s.fNormHistogram.DrawCopy("colz")
        h.GetZaxis().SetTitleOffset(1.4)

        globalList.append(c)
        globalList.append(h)

        if self.fTrigger:
            cname = "{0}_{1}_canvas".format(self.fTrigger, s.fUncertainty.GetName())
        else:
            cname = "{0}_canvas".format(s.fUncertainty.GetName())

        c = ROOT.TCanvas(cname, cname)
        c.SetRightMargin(0.18)
        self.fCanvases.append(c)
        c.cd()
        h = s.fUncertainty.DrawCopy("colz")
        h.GetZaxis().SetTitleOffset(1.4)

        globalList.append(c)
        globalList.append(h)

    def PlotSpectrum3D(self, s):
        print("PlotSpectrum3D not implemented!")

    def GenerateSpectra(self):
        for binMultiSet in self.fBinMultiSets.itervalues():
            for binSet in binMultiSet.fBinSets.itervalues():
                for s in binSet.fSpectra.itervalues():
                    if len(s.fAxis) == 1:
                        self.GenerateSpectrum1D(s)
                    elif len(s.fAxis) == 2:
                        self.GenerateSpectrum2D(s)
                    elif len(s.fAxis) == 3:
                        self.GenerateSpectrum3D(s)
                    else:
                        print("Not able to generate spectra with dim > 3!")
                    s.GenerateNormalizedSpectrum(self.fEvents, self.fIsWeighted)

    def GenerateSpectrum1DInvMassFit(self, s):
        if s.fAnalysisType == AnalysisType.InvMassFit:
            for bin in s.fBinSet.fBins:
                self.ProcessInvMassFitSpectrumBin(s, bin, s.fBinSet)
        elif s.fAnalysisType == AnalysisType.LikeSignFit:
            if not s.fLikeSignSubtractedBinSet:
                print("No LS subtracted invariant mass bin sets was found!")
                return
            for bin in s.fLikeSignSubtractedBinSet.fBins:
                self.ProcessInvMassFitSpectrumBin(s, bin, s.fLikeSignSubtractedBinSet)

    def ProcessInvMassFitSpectrumBin(self, s, bin, binSet):
        w = s.fEfficiencyWeight.GetEfficiencyWeightTH1ForPt(bin.GetBinCenter("d_pt"))

        xbin = s.fHistogram.GetXaxis().FindBin(bin.GetBinCenter(s.fAxis[0].fName))
        if bin.fMassFitter is None:
            print("The bin printed below does not have a mass fitter!")
            bin.Print()
            return
        if not bin.fMassFitter.FitSuccessfull():
            print("The bin printed below does not have a successful invariant mass fit!")
            bin.Print()
            return
        signal = bin.fMassFitter.GetSignal() * w
        signal_unc = bin.fMassFitter.GetSignalError() * w

        s.fHistogram.SetBinContent(xbin, signal)
        s.fHistogram.SetBinError(xbin, signal_unc)
        if signal > 0:
            s.fUncertainty.SetBinContent(xbin, signal_unc / signal)

        if bin.fMassFitter:
            if s.fBackground:
                s.fBackground.SetBinContent(xbin, bin.fMassFitter.GetBackground(2))
                s.fBackground.SetBinError(xbin, bin.fMassFitter.GetBackgroundError(2))
            if s.fMass:
                s.fMass.SetBinContent(xbin, bin.fMassFitter.GetSignalMean())
                s.fMass.SetBinError(xbin, bin.fMassFitter.GetSignalMeanError())
            if s.fMassWidth:
                s.fMassWidth.SetBinContent(xbin, bin.fMassFitter.GetSignalWidth())
                s.fMassWidth.SetBinError(xbin, bin.fMassFitter.GetSignalWidthError())


    def GenerateInvMassWidonws(self, invMassHisto, binSBL_1, binSBL_2, binSBR_1, binSBR_2, binSig_1, binSig_2):
        print("Generating invariant mass windows for {0}".format(invMassHisto.GetName()))
        sideBandWindowInvMassHisto = invMassHisto.Clone(invMassHisto.GetName().replace("InvMass", "InvMassSBWindow"))
        signalWindowInvMassHisto = invMassHisto.Clone(invMassHisto.GetName().replace("InvMass", "InvMassSigWindow"))
        for xbin in xrange(0, invMassHisto.GetNbinsX() + 2):
            if not ((xbin >= binSBL_1 and xbin < binSBL_2) or (xbin >= binSBR_1 and xbin < binSBR_2)):
                sideBandWindowInvMassHisto.SetBinContent(xbin, 0)
                sideBandWindowInvMassHisto.SetBinError(xbin, 0)
            if not (xbin >= binSig_1 and xbin < binSig_2):
                signalWindowInvMassHisto.SetBinContent(xbin, 0)
                signalWindowInvMassHisto.SetBinError(xbin, 0)

        return signalWindowInvMassHisto, sideBandWindowInvMassHisto

    def GenerateSpectrum1DSideBandMethod(self, s):
        if s.fSkipBins:
            print("I will skip the following bins: {0}".format(s.fSkipBins))
        else:
            print("I will not skip any bin")
        for ibin, bin in enumerate(s.fBinSet.fBins):
            if s.fSkipBins and ibin in s.fSkipBins:
                print("Skipping bin {0} as requested".format(bin.GetTitle()))
                continue
            if not bin.fMassFitter or not bin.fMassFitter.FitSuccessfull():
                print("Skipping bin {0} because fit was unsuccessful".format(bin.GetTitle()))
                continue
            print("Bin: {0}".format(bin.GetTitle()))
            sigma = bin.fMassFitter.GetSignalWidth()
            mean = bin.fMassFitter.GetSignalMean()
            print("Sigma={0:.3f}, Mean={1:.3f}".format(sigma, mean))

            binSBL_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(mean - s.fSideBandMaxSigmas * sigma)
            binSBL_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(mean - s.fSideBandMinSigmas * sigma)
            binSBR_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(mean + s.fSideBandMinSigmas * sigma)
            binSBR_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(mean + s.fSideBandMaxSigmas * sigma)
            binSig_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(mean - s.fBinCountSignalSigmas * sigma)
            binSig_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(mean + s.fBinCountSignalSigmas * sigma)
            if binSBL_1 < 1:
                binSBL_1 = 1
            if binSBR_2 > bin.fBinCountAnalysisHisto.GetXaxis().GetNbins():
                binSBR_2 = bin.fBinCountAnalysisHisto.GetXaxis().GetNbins()
            effSigma1 = (mean - bin.fBinCountAnalysisHisto.GetXaxis().GetBinLowEdge(binSig_1)) / sigma
            effSigma2 = (bin.fBinCountAnalysisHisto.GetXaxis().GetBinUpEdge(binSig_2) - mean) / sigma
            print("The left effective sigma is {0:.3f}; the right effective sigma is {1:.3f}".format(effSigma1, effSigma2))

            binSBL_1_invMass = bin.fInvMassHisto.GetXaxis().FindBin(mean - s.fSideBandMaxSigmas * sigma)
            binSBL_2_invMass = bin.fInvMassHisto.GetXaxis().FindBin(mean - s.fSideBandMinSigmas * sigma)
            binSBR_1_invMass = bin.fInvMassHisto.GetXaxis().FindBin(mean + s.fSideBandMinSigmas * sigma)
            binSBR_2_invMass = bin.fInvMassHisto.GetXaxis().FindBin(mean + s.fSideBandMaxSigmas * sigma)
            binSig_1_invMass = bin.fInvMassHisto.GetXaxis().FindBin(mean - s.fBinCountSignalSigmas * sigma)
            binSig_2_invMass = bin.fInvMassHisto.GetXaxis().FindBin(mean + s.fBinCountSignalSigmas * sigma)
            (signalWindowInvMassHisto, sideBandWindowInvMassHisto) = self.GenerateInvMassWidonws(bin.fInvMassHisto, binSBL_1_invMass, binSBL_2_invMass, binSBR_1_invMass, binSBR_2_invMass, binSig_1_invMass, binSig_2_invMass)

            s.fSideBandWindowInvMassHistos[sideBandWindowInvMassHisto.GetName()] = sideBandWindowInvMassHisto
            s.fSignalWindowInvMassHistos[signalWindowInvMassHisto.GetName()] = signalWindowInvMassHisto

            if bin.fBinCountAnalysisHisto.GetDimension() == 2:
                sbL = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SideBandWindowL_{1}".format(s.fName, bin.GetName()), binSBL_1, binSBL_2, "e")
                sbR = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SideBandWindowR_{1}".format(s.fName, bin.GetName()), binSBR_1, binSBR_2, "e")
                sig = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SignalWindow_{1}".format(s.fName, bin.GetName()), binSig_1, binSig_2, "e")
            elif bin.fBinCountAnalysisHisto.GetDimension() == 3:
                bin.fBinCountAnalysisHisto.GetXaxis().SetRange(binSBL_1, binSBL_2)
                sbL = bin.fBinCountAnalysisHisto.Project3D("zye")
                sbL.SetName("{0}_SideBandWindowL_{1}".format(s.fName, bin.GetName()))
                bin.fBinCountAnalysisHisto.GetXaxis().SetRange(binSBR_1, binSBR_2)
                sbR = bin.fBinCountAnalysisHisto.Project3D("zye")
                sbR.SetName("{0}_SideBandWindowR_{1}".format(s.fName, bin.GetName()))
                bin.fBinCountAnalysisHisto.GetXaxis().SetRange(binSig_1, binSig_2)
                sig = bin.fBinCountAnalysisHisto.Project3D("zye")
                sig.SetName("{0}_SignalWindow_{1}".format(s.fName, bin.GetName()))
                bin.fBinCountAnalysisHisto.GetXaxis().SetRange(1, bin.fBinCountAnalysisHisto.GetNbinsX())

            sbTotal = sbL.Clone("{0}_SideBandWindow_{1}".format(s.fName, bin.GetName()))
            sbTotal.Add(sbR)

            # peakAreaBkgError = ROOT.Double(0.)
            # peakAreaBkg = bin.fMassFitter.GetBackgroundAndError(peakAreaBkgError, s.fBinCountSignalSigmas)
            intSigErr = ROOT.Double(0.)
            if sig.GetDimension() == 1:
                intSig = sig.IntegralAndError(0, -1, intSigErr)
            elif sig.GetDimension() == 2:
                intSig = sig.IntegralAndError(0, -1, 0, -1, intSigErr)
            peakAreaBkgError = math.sqrt(intSigErr ** 2 + (bin.fMassFitter.GetSignalError(effSigma1) / 2) ** 2 + (bin.fMassFitter.GetSignalError(effSigma2) / 2) ** 2)
            peakAreaBkg = intSig - bin.fMassFitter.GetSignal(effSigma1) / 2 - bin.fMassFitter.GetSignal(effSigma2) / 2
            if sbL.GetDimension() == 1:
                print("The background in side bands is: {0} + {1} = {2}".format(sbL.Integral(0, -1), sbR.Integral(0, -1), sbTotal.Integral(0, -1)))
            elif sbL.GetDimension() == 2:
                print("The background in side bands is: {0} + {1} = {2}".format(sbL.Integral(0, -1, 0, -1), sbR.Integral(0, -1, 0, -1), sbTotal.Integral(0, -1, 0, -1)))
            print("The estimated background in the signal window is {0} +/- {1}".format(peakAreaBkg, peakAreaBkgError))
            print("The total signal+background is {0}, which is the same from the invariant mass plot {1} or summing signal and background {2}".format(intSig, bin.fInvMassHisto.Integral(binSig_1, binSig_2), bin.fMassFitter.GetSignal(s.fBinCountSignalSigmas) + bin.fMassFitter.GetBackground(s.fBinCountSignalSigmas)))

            sbTotalIntegralError = ROOT.Double(0)
            if sbTotal.GetDimension() == 1:
                sbTotalIntegral = sbTotal.IntegralAndError(0, -1, sbTotalIntegralError)
            elif sbTotal.GetDimension() == 2:
                sbTotalIntegral = sbTotal.IntegralAndError(0, -1, 0, -1, sbTotalIntegralError)
            if sbTotalIntegral > 0:
                if sbTotal.GetDimension() == 1:
                    for xbin in xrange(0, sbTotal.GetNbinsX() + 2):
                        if sbTotal.GetBinContent(xbin) == 0:
                            continue
                        # Error propagation
                        error2_1 = peakAreaBkgError ** 2 * sbTotal.GetBinContent(xbin) ** 2
                        error2_2 = peakAreaBkg ** 2 / sbTotalIntegral ** 2 * sbTotal.GetBinError(xbin) ** 2 * (sbTotalIntegral - sbTotal.GetBinContent(xbin)) ** 2
                        error2_3 = peakAreaBkg ** 2 / sbTotalIntegral ** 2 * sbTotal.GetBinContent(xbin) ** 2 * (sbTotalIntegralError ** 2 - sbTotal.GetBinError(xbin) ** 2)
                        # print("bin {0}: error1 = {1}, error2 = {2}, error3 = {3}".format(xbin, math.sqrt(error2_1)/sbTotalIntegral, math.sqrt(error2_2)/sbTotalIntegral, math.sqrt(error2_3)/sbTotalIntegral))
                        error = math.sqrt(error2_1 + error2_2 + error2_3) / sbTotalIntegral
                        cont = sbTotal.GetBinContent(xbin) * peakAreaBkg / sbTotalIntegral
                        sbTotal.SetBinError(xbin, error)
                        sbTotal.SetBinContent(xbin, cont)
                elif sbTotal.GetDimension() == 2:
                    for xbin in xrange(0, sbTotal.GetNbinsX() + 2):
                        for ybin in xrange(0, sbTotal.GetNbinsY() + 2):
                            if sbTotal.GetBinContent(xbin, ybin) == 0: continue
                            # Error propagation
                            error2_1 = peakAreaBkgError ** 2 * sbTotal.GetBinContent(xbin, ybin) ** 2
                            error2_2 = peakAreaBkg ** 2 / sbTotalIntegral ** 2 * sbTotal.GetBinError(xbin, ybin) ** 2 * (sbTotalIntegral - sbTotal.GetBinContent(xbin, ybin)) ** 2
                            error2_3 = peakAreaBkg ** 2 / sbTotalIntegral ** 2 * sbTotal.GetBinContent(xbin, ybin) ** 2 * (sbTotalIntegralError ** 2 - sbTotal.GetBinError(xbin, ybin) ** 2)
                            # print("bin {0}: error1 = {1}, error2 = {2}, error3 = {3}".format(xbin, math.sqrt(error2_1)/sbTotalIntegral, math.sqrt(error2_2)/sbTotalIntegral, math.sqrt(error2_3)/sbTotalIntegral))
                            error = math.sqrt(error2_1 + error2_2 + error2_3) / sbTotalIntegral
                            cont = sbTotal.GetBinContent(xbin, ybin) * peakAreaBkg / sbTotalIntegral
                            sbTotal.SetBinError(xbin, ybin, error)
                            sbTotal.SetBinContent(xbin, ybin, cont)

            if sbTotal.GetDimension() == 1:
                integralError = ROOT.Double(0.)
                integral = sbTotal.IntegralAndError(0, -1, integralError)
            elif sbTotal.GetDimension() == 2:
                integralError = ROOT.Double(0.)
                integral = sbTotal.IntegralAndError(0, -1, 0, -1, integralError)
            print("The total normalized side-band background is {0} +/- {1}".format(integral, integralError))

            w = s.fEfficiencyWeight.GetEfficiencyWeightTH1ForPt(bin.GetBinCenter("d_pt"))
            if w != 1:
                sbTotal.Scale(w)
                sig.Scale(w)

            scaleGaussianLimit = 1.0 - math.erfc(effSigma1 / math.sqrt(2.0)) / 2 - math.erfc(effSigma2 / math.sqrt(2.0)) / 2
            print("Scaling for the Gaussian limit ({0} sigmas): {1}".format(s.fBinCountSignalSigmas, scaleGaussianLimit))
            sbTotal.Scale(1.0 / scaleGaussianLimit)
            sig.Scale(1.0 / scaleGaussianLimit)

            sbTotal.SetTitle(bin.GetTitle())
            s.fSideBandHistograms.append(sbTotal)
            s.fSideBandLeftHistograms.append(sbL)
            s.fSideBandRightHistograms.append(sbR)
            s.fSideBandWindowTotalHistogram.Add(sbTotal)

            sig.SetTitle(bin.GetTitle())
            s.fSignalHistograms.append(sig)
            s.fSignalWindowTotalHistogram.Add(sig)

        SBinvMassName = "{0}_SideBand_{1}".format(s.fBinSet.fName, s.fName)
        self.PlotInvMassPlotsBinSet(SBinvMassName, s.fBinSet.fBins, None, s)

        s.fHistogram.Add(s.fSignalWindowTotalHistogram)
        s.fHistogram.Add(s.fSideBandWindowTotalHistogram, -1)
        if s.fBackground: s.fBackground.Add(s.fSideBandWindowTotalHistogram)

        if s.fHistogram.GetDimension() == 1:
            for xbin in xrange(0, s.fHistogram.GetNbinsX() + 2):
                if s.fHistogram.GetBinContent(xbin) > 0:
                    s.fUncertainty.SetBinContent(xbin, s.fHistogram.GetBinError(xbin) / s.fHistogram.GetBinContent(xbin))
                else:
                    s.fUncertainty.SetBinContent(xbin, 0)
        if s.fHistogram.GetDimension() == 2:
            for xbin in xrange(0, s.fHistogram.GetNbinsX() + 2):
                for ybin in xrange(0, s.fHistogram.GetNbinsY() + 2):
                    if s.fHistogram.GetBinContent(xbin, ybin) > 0:
                        s.fUncertainty.SetBinContent(xbin, ybin, s.fHistogram.GetBinError(xbin, ybin) / s.fHistogram.GetBinContent(xbin, ybin))
                    else:
                        s.fUncertainty.SetBinContent(xbin, ybin, 0)

    def GenerateSpectrum1DTruth(self, s):
        # The truth spectrum is already done, only need to apply the efficiency
        for ibin in xrange(0, s.fHistogram.GetNbinsX() + 2):
            if not s.fHistogram.GetBinContent(ibin) > 0:
                continue
            w = s.fEfficiencyWeight.GetEfficiencyWeightTH1ForPt(s.fHistogram.GetXaxis().GetBinCenter(ibin))
            s.fHistogram.SetBinContent(ibin, s.fHistogram.GetBinContent(ibin) * w)
            s.fHistogram.SetBinError(ibin, s.fHistogram.GetBinError(ibin) * w)
            s.fUncertainty.SetBinContent(ibin, s.fHistogram.GetBinError(ibin) / s.fHistogram.GetBinContent(ibin))
        if self.fDMeson in s.fBinSet.fNeedInvMass:
            for xbin, bin in enumerate(s.fBinSet.fBins):
                if bin.fMassFitter:
                    if s.fMass:
                        s.fMass.SetBinContent(xbin + 1, bin.fMassFitter.GetSignalMean())
                        s.fMass.SetBinError(xbin + 1, bin.fMassFitter.GetSignalMeanError())
                    if s.fMassWidth:
                        s.fMassWidth.SetBinContent(xbin + 1, bin.fMassFitter.GetSignalWidth())
                        s.fMassWidth.SetBinError(xbin + 1, bin.fMassFitter.GetSignalWidthError())

    def GenerateSpectrum1D(self, s):
        print("Generating spectrum {0}".format(s.fName))
        if s.fAnalysisType == AnalysisType.SideBand:
            self.GenerateSpectrum1DSideBandMethod(s)
        elif s.fAnalysisType == AnalysisType.LikeSign or s.fAnalysisType == AnalysisType.LikeSignFit:
            print("Like sign analysis not implemented!")
        elif s.fAnalysisType == AnalysisType.InvMassFit:
            self.GenerateSpectrum1DInvMassFit(s)
        elif s.fAnalysisType == AnalysisType.Truth:
            self.GenerateSpectrum1DTruth(s)
        else:
            print("Analysis type {0} not recognized!".format(s.fAnalysisType))

    def GenerateSpectrum2D(self, s):
        print("Generating spectrum {0}".format(s.fName))
        if s.fAnalysisType == AnalysisType.InvMassFit:
            self.GenerateSpectrum2DInvMassFit(s)
        if s.fAnalysisType == AnalysisType.SideBand:
            self.GenerateSpectrum1DSideBandMethod(s)
        elif s.fAnalysisType == AnalysisType.Truth:
            self.GenerateSpectrum2DTruth(s)
        else:
            print("Analysis type {0} not implemented for 2D spectra!".format(s.fAnalysisType))

    def GenerateSpectrum2DInvMassFit(self, s):
        for binSetName in s.fBins:
            for bin in self.fBinMultiSets[s.fJetType, s.fJetRadius].fBinSets[binSetName].fBins:
                if bin.fMassFitter is None:
                    print("The bin printed below does not have a mass fitter!")
                    bin.Print()
                    continue
                xbin = s.fHistogram.GetXaxis().FindBin(bin.GetBinCenter(s.fAxis[0].fName))
                ybin = s.fHistogram.GetYaxis().FindBin(bin.GetBinCenter(s.fAxis[1].fName))
                signal = bin.fMassFitter.GetSignal()
                signal_unc = bin.fMassFitter.GetSignalError()

                s.fHistogram.SetBinContent(xbin, ybin, signal)
                s.fHistogram.SetBinError(xbin, ybin, signal_unc)
                s.fBackground.SetBinContent(xbin, ybin, bin.fMassFitter.GetBackground())
                s.fBackground.SetBinError(xbin, ybin, bin.fMassFitter.GetBackgroundError())
                s.fUncertainty.SetBinContent(xbin, ybin, signal_unc / signal)
                s.fMass.SetBinContent(xbin, ybin, bin.fMassFitter.GetSignalMean())
                s.fMass.SetBinError(xbin, ybin, bin.fMassFitter.GetSignalMeanError())
                s.fMassWidth.SetBinContent(xbin, ybin, bin.fMassFitter.GetSignalWidth())
                s.fMassWidth.SetBinError(xbin, ybin, bin.fMassFitter.GetSignalWidthError())

    def GenerateSpectrum2DTruth(self, s):
        w = 1
        for xbin in xrange(0, s.fHistogram.GetNbinsX() + 2):
            if s.fAxis[0].fName == "d_pt":
                w = s.fEfficiencyWeight.GetEfficiencyWeightTH1ForPt(s.fHistogram.GetXaxis().GetBinCenter(xbin))
            for ybin in xrange(0, s.fHistogram.GetNbinsY() + 2):
                if s.fAxis[1].fName == "d_pt":
                    w = s.fEfficiencyWeight.GetEfficiencyWeightTH1ForPt(s.fHistogram.GetYaxis().GetBinCenter(ybin))
                if not s.fHistogram.GetBinContent(xbin, ybin) > 0:
                    continue
                s.fHistogram.SetBinContent(xbin, ybin, s.fHistogram.GetBinContent(xbin, ybin) * w)
                s.fHistogram.SetBinError(xbin, ybin, s.fHistogram.GetBinError(xbin, ybin) * w)
                s.fUncertainty.SetBinContent(xbin, ybin, s.fHistogram.GetBinError(xbin, ybin) / s.fHistogram.GetBinContent(xbin, ybin))

    def GenerateSpectrum3D(self, s):
        print("GenerateSpectrum3D not implemented!")

    def DrawFitResults(self, bin):
        if bin.fMassFitter is None or not bin.fMassFitter.FitSuccessfull():
            return

        bin.fMassFitter.Draw("same");

        chi2Text = bin.fMassFitter.GetChisquareString().Data()

        paveSig = ROOT.TPaveText(0.165, 0.795, 0.490, 0.92, "NB NDC")
        globalList.append(paveSig)
        paveSig.SetBorderSize(0)
        paveSig.SetFillStyle(0)
        paveSig.SetTextFont(43)
        paveSig.SetTextSize(14)
        paveSig.SetTextAlign(13)
        paveSig.AddText("{0}, {1}".format(bin.fMassFitter.GetSignalString().Data(),
                                          bin.fMassFitter.GetBackgroundString().Data()))
        paveSig.AddText("{0}, {1}, {2}".format(bin.fMassFitter.GetSignalOverBackgroundString().Data(),
                                          bin.fMassFitter.GetSignalOverSqrtSignalBackgroundString().Data(),
                                          chi2Text))
        paveSig.Draw()

        paveFit = ROOT.TPaveText(0.48, 0.51, 0.97, 0.77, "NB NDC")
        globalList.append(paveFit)
        paveFit.SetBorderSize(0)
        paveFit.SetFillStyle(0)
        paveFit.SetTextFont(43)
        paveFit.SetTextSize(14)
        paveFit.SetTextAlign(23)

        paveFit.AddText(bin.fMassFitter.GetSignalMeanString().Data())
        paveFit.AddText(bin.fMassFitter.GetSignalWidthString().Data())
        paveFit.AddText(bin.fMassFitter.GetBkgPar1String().Data())
        paveFit.AddText(bin.fMassFitter.GetTotalEntriesString().Data())
        paveFit.Draw()

    def FitInvMassPlots(self):
        for binMultiSet in self.fBinMultiSets.itervalues():
            for binSet in binMultiSet.fBinSets.itervalues():
                if self.fDMeson in binSet.fNeedInvMass:
                    self.FitInvMassPlotsBinSet(binSet.fName, binSet.fBins, binSet.fFitOptions)

    def FitInvMassPlotsBinSet(self, name, bins, fitOptions):
        print("Fitting {0}".format(name))

        for i, bin in enumerate(bins):
            if not bin.fInvMassHisto:
                continue
            if self.fTrigger:
                fitterName = "InvMass_{0}_{1}_{2}_fitter".format(self.fTrigger, self.fDMeson, bin.GetName())
            else:
                fitterName = "InvMass_{0}_{1}_fitter".format(self.fDMeson, bin.GetName())

            fitter = self.CreateMassFitter(fitterName, bin.fInvMassHisto)
            bin.SetMassFitter(fitter)

            print("Fitting bin {0}".format(bin.GetTitle()))

            fitter.Fit(bin.fInvMassHisto, fitOptions)

    def PlotInvMassPlots(self):
        for binMultiSet in self.fBinMultiSets.itervalues():
            for binSet in binMultiSet.fBinSets.itervalues():
                if self.fDMeson in binSet.fNeedInvMass:
                    self.PlotInvMassPlotsBinSet(binSet.fName, binSet.fBins)

    def PlotInvMassLikeSign(self, bin):
        hls = bin.fInvMassHisto.DrawCopy("hist same")
        hls.SetLineColor(ROOT.kMagenta)
        hls.SetLineStyle(2)
        hls.SetLineWidth(2)
        hls.SetFillStyle(0)
        globalList.append(hls)

        return hls

    def PlotInvMassSideBands(self, bin, spectrum):
        if not spectrum:
            return None
        invMassSBHistoName = bin.fInvMassHisto.GetName().replace("InvMass", "InvMassSBWindow")
        invMassSigHistoName = bin.fInvMassHisto.GetName().replace("InvMass", "InvMassSigWindow")
        if not invMassSBHistoName in spectrum.fSideBandWindowInvMassHistos.keys() or not invMassSigHistoName in spectrum.fSignalWindowInvMassHistos.keys():
            return None
        hsb = spectrum.fSideBandWindowInvMassHistos[invMassSBHistoName].DrawCopy("hist")
        hsb.SetFillColorAlpha(ROOT.kGreen + 2, 0.4)
        hsb.SetFillStyle(1001)
        hsb.SetLineColorAlpha(ROOT.kGreen + 2, 0.4)
        hsig = spectrum.fSignalWindowInvMassHistos[invMassSigHistoName].DrawCopy("hist same")
        hsig.SetFillColorAlpha(ROOT.kRed + 2, 0.4)
        hsig.SetFillStyle(1001)
        hsig.SetLineColorAlpha(ROOT.kRed + 2, 0.4)
        globalList.append(hsb)
        globalList.append(hsig)

        return hsb

    def PlotInvMassPlotsBinSet(self, name, bins, LS_bins=None, spectrum=None):
        if self.fTrigger:
            cname = "{0}_{1}".format(self.fTrigger, name)
        else:
            cname = name
        nbins = len(bins)
        if spectrum and spectrum.fSkipBins:
            nbins -= len(spectrum.fSkipBins)
        c = DMesonJetUtils.GenerateMultiCanvas(cname, nbins)
        self.fCanvases.append(c)
        globalList.append(c)
        icanvas = 1
        for i, bin in enumerate(bins):
            if spectrum and spectrum.fSkipBins and i in spectrum.fSkipBins:
                continue
            if not bin.fInvMassHisto:
                continue
            pad = c.cd(icanvas)
            pad.SetLeftMargin(0.12)
            pad.SetRightMargin(0.05)
            pad.SetTopMargin(0.08)
            pad.SetBottomMargin(0.13)
            if spectrum:
                SB = self.PlotInvMassSideBands(bin, spectrum)
            else:
                SB = None
            if SB:
                h = bin.fInvMassHisto.DrawCopy("same")
                globalList.append(h)
                SB.SetMaximum(h.GetMaximum() * 1.8)
                h = SB
            else:
                h = bin.fInvMassHisto.DrawCopy()
                globalList.append(h)
                h.SetMaximum(h.GetMaximum() * 1.8)
            h.GetXaxis().SetTitleFont(43)
            h.GetXaxis().SetTitleOffset(2.3)
            h.GetXaxis().SetTitleSize(19)
            h.GetXaxis().SetLabelFont(43)
            h.GetXaxis().SetLabelOffset(0.009)
            h.GetXaxis().SetLabelSize(18)
            h.GetYaxis().SetTitleFont(43)
            h.GetYaxis().SetTitleOffset(2.3)
            h.GetYaxis().SetTitleSize(19)
            h.GetYaxis().SetLabelFont(43)
            h.GetYaxis().SetLabelOffset(0.009)
            h.GetYaxis().SetLabelSize(18)
            htitle = ROOT.TPaveText(0.12, 0.99, 0.95, 0.93, "NB NDC")
            htitle.SetBorderSize(0)
            htitle.SetFillStyle(0)
            htitle.SetTextFont(43)
            htitle.SetTextSize(18)
            htitle.AddText(bin.GetTitle())
            htitle.Draw()
            globalList.append(htitle)
            self.DrawFitResults(bin)
            if LS_bins:
                self.PlotInvMassLikeSign(LS_bins[i])
            icanvas += 1

class DMesonJetAnalysis:
    def __init__(self, name):
        self.fName = name
        self.fAnalysisEngine = []
        self.fCanvases = []
        self.fJets = None

    def SetProjector(self, projector):
        self.fProjector = projector

    def StartAnalysis(self, collision, config):
        self.fCollision = collision
        self.fJets = config["jets"]
        binMultiSet = BinSet.BinMultiSet()
        for binLists in config["binLists"]:
            if len(binLists["active_mesons"]) == 0:
                continue
            limitSetList = []
            axis = []
            for name, binList in binLists["bins"].iteritems():
                limitSetList.append((name, binList))
                axis.append(Axis.Axis(name, binList, "", True))
            if "cuts" in binLists:
                cuts = binLists["cuts"]
            else:
                cuts = []
            if "bin_count_analysis" in binLists:
                bin_count_analysis = binLists["bin_count_analysis"]
            else:
                bin_count_analysis = None
            if "efficiency" in binLists and binLists["efficiency"]:
                effWeight = DMesonJetProjectors.EfficiencyWeightCalculator("{0}/{1}".format(self.fProjector.fInputPath, binLists["efficiency"]["file_name"]), binLists["efficiency"]["list_name"], binLists["efficiency"]["object_name"])
                fitOptions = "0 WL S"
            else:
                effWeight = DMesonJetProjectors.SimpleWeight()
                if self.fProjector.fMergingType == "simple_sum":
                    fitOptions = "0 L S"
                else:
                    fitOptions = "0 WL S"
            if not "spectra" in binLists:
                binLists["spectra"] = []
            binMultiSet.AddBinSet(BinSet.BinSet(binLists["name"], binLists["title"], binLists["need_inv_mass"], limitSetList, binLists["spectra"], axis, cuts, bin_count_analysis, effWeight, fitOptions))

        for trigger in config["trigger"]:
            for d_meson in config["d_meson"]:
                if trigger:
                    print("Projecting trigger {0}, D meson {1}".format(trigger, d_meson))
                else:
                    print("Projecting D meson {0}".format(d_meson))
                eng = DMesonJetAnalysisEngine(collision, trigger, d_meson,
                                              binMultiSet, config["n_mass_bins"], config["min_mass"], config["max_mass"],
                                              self.fJets, self.fProjector)
                self.fAnalysisEngine.append(eng)
                eng.DoProjections()

        for eng in self.fAnalysisEngine:
            if not "LikeSign" in eng.fDMeson:
                eng.Start(self)
                eng.CompareSpectra()

        for jetDef in self.fJets:
            self.CompareDMesons(config["binLists"], jetDef)

        jetDef = dict()
        jetDef["type"] = None
        jetDef["radius"] = None
        jetDef["title"] = None
        self.CompareDMesons(config["binLists"], jetDef)

    def CompareDMesons(self, binLists, jetDef):
        if len(self.fAnalysisEngine) <= 1:
            return
        allSpectrumNames = set()
        for binList in binLists:
            if len(binList["active_mesons"]) == 0:
                continue
            if len(binList["bins"]) != 1:
                continue
            for s in binList["spectra"]:
                allSpectrumNames.add(s["name"])
        for spectrum_name in allSpectrumNames:
            spectraToCompare = []
            for eng in self.fAnalysisEngine:
                for binList in binLists:
                    if not binList["name"] in eng.fBinMultiSets[jetDef["type"], jetDef["radius"]].fBinSets:
                        continue
                    if len(binList["active_mesons"]) == 0:
                        continue
                    if len(binList["bins"]) != 1:
                        continue
                    if not eng.fDMeson in binList["active_mesons"]:
                        continue
                    for s in binList["spectra"]:
                        if not s["name"] == spectrum_name:
                            continue
                        if not eng.fDMeson in s["active_mesons"]:
                            continue
                        if "suffix" in s:
                            suffix = s["suffix"]
                        else:
                            suffix = None
                        sname = '_'.join(obj for obj in [eng.fDMeson, jetDef["type"], jetDef["radius"], s["name"], suffix] if obj)
                        binSet = eng.fBinMultiSets[jetDef["type"], jetDef["radius"]].fBinSets[binList["name"]]
                        h = binSet.fSpectra[sname].fHistogram
                        if not h:
                            continue
                        if "MCTruth" in eng.fDMeson:
                            continue
                        h_copy = h.Clone("{0}_copy".format(h.GetName()))
                        if s["title"]:
                            if binList["title"]:
                                h_copy.SetTitle("{0}, {1}".format(s["title"], binList["title"]))
                            else:
                                h_copy.SetTitle("{0}".format(s["title"]))
                        else:
                            if binList["title"]:
                                h_copy.SetTitle("{0}".format(binList["title"]))
                            else:
                                h_copy.SetTitle(eng.fDMeson)
                        globalList.append(h_copy)
                        spectraToCompare.append(h_copy)
            if len(spectraToCompare) > 1:
                cname = '_'.join(obj for obj in [jetDef["type"], jetDef["radius"], spectrum_name, "SpectraComparison"] if obj)
                comp = DMesonJetCompare.DMesonJetCompare(cname)
                results = comp.CompareSpectra(spectraToCompare[0], spectraToCompare[1:])
                for obj in results:
                    if isinstance(obj, ROOT.TCanvas):
                        self.fCanvases.append(obj)
                        obj.cd()
                        pave = ROOT.TPaveText(0.12, 0.70, 0.40, 0.85, "NB NDC")
                        pave.SetTextAlign(11)
                        pave.SetFillStyle(0)
                        pave.SetBorderSize(0)
                        pave.SetTextFont(43)
                        pave.SetTextSize(15)
                        pave.AddText(self.fCollision)
                        if jetDef["title"]: pave.AddText(jetDef["title"])
                        pave.Draw()
                        globalList.append(pave)
                    globalList.append(obj)

    def SaveRootFile(self, path):
        file = ROOT.TFile.Open("{0}/{1}.root".format(path, self.fName), "recreate")
        for eng in self.fAnalysisEngine:
            eng.SaveRootFile(file)
        file.Close()

    def SavePlots(self, path, format):
        fullPath = "{0}/{1}/{2}".format(path, self.fName, format)
        if not os.path.isdir(fullPath):
            os.makedirs(fullPath)
        for eng in self.fAnalysisEngine:
            eng.SavePlots(fullPath, format)

        for c in self.fCanvases:
            c.SaveAs("{0}/{1}.{2}".format(fullPath, c.GetName(), format))
