#!/usr/bin/env python
# python program to perform a D meson jet analysis

import math
import array
import copy
import collections
import os

import ROOT

import DMesonJetProjectors
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
            h = s.fNormHistogram.Clone("{0}_copy".format(s.fNormHistogram.GetName()))
            if s.fTitle:
                h.SetTitle(s.fTitle)
            globalList.append(h)
            spectraToCompare.append(h)
        if len(spectraToCompare) < 2:
            return
        results = DMesonJetUtils.CompareSpectra(spectraToCompare[0], spectraToCompare[1:], "{0}_{1}_{2}_{3}_SpectraComparison".format(self.fDMeson, binMultiSet.fJetType, binMultiSet.fJetRadius, axisName), "", "hist")
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

    def CreateMassFitter(self, name):
        if "D0" in self.fDMeson:
            minMass = self.fMinMass
            maxMass = self.fMaxMass
            minFitRange = self.fMinMass
            maxFitRange = self.fMaxMass
            startingSigma = 0.01
            startingSigmaBkg = -1
            DMesonType = ROOT.MassFitter.kDzeroKpi
        elif self.fDMeson == "DStar":
            print("Not implemented for DStar!")

        fitter = ROOT.MassFitter(name, DMesonType, minMass, maxMass)
        fitter.GetFitFunction().SetParameter(1, startingSigmaBkg)
        fitter.GetFitFunction().SetParameter(4, startingSigma)
        fitter.SetFitRange(minFitRange, maxFitRange)
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
        if not "MCTruth" in self.fDMeson:
            self.FitInvMassPlots()
        if not "BackgroundOnly" in self.fDMeson:
            self.GenerateSpectra()

        if not "MCTruth" in self.fDMeson:
            self.PlotInvMassPlots()
        if not "BackgroundOnly" in self.fDMeson:
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
        c = ROOT.TCanvas("{0}_canvas".format(s.fNormHistogram.GetName()), s.fNormHistogram.GetName())
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
        c = ROOT.TCanvas("{0}_canvas".format(s.fUncertainty.GetName()), s.fUncertainty.GetName())
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
            c = ROOT.TCanvas("{0}_canvas".format(s.fMass.GetName()), s.fMass.GetName())
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
            c = ROOT.TCanvas("{0}_canvas".format(s.fMassWidth.GetName()), s.fMassWidth.GetName())
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
            self.PlotBackgroundVsSignalSpectra("{0}_TotalBkgVsSig".format(s.fName), s.fUnlikeSignTotalHistogram, "Unlike Sign", s.fLikeSignTotalHistogram, "Like Sign")
        if s.fUnlikeSignHistograms and s.fLikeSignHistograms:
            self.PlotMultiCanvasBkgVsSigSpectra("{0}_BkgVsSig".format(s.fName), s.fUnlikeSignHistograms, "Unlike Sign", s.fLikeSignHistograms, "Like Sign")

        if s.fSideBandWindowTotalHistogram:
            self.PlotBackgroundVsSignalSpectra("{0}_TotalBkgVsSig".format(s.fName), s.fSignalWindowTotalHistogram, "Sig. Window", s.fSideBandWindowTotalHistogram, "SB Window")
        if s.fSignalHistograms and s.fSideBandHistograms:
            self.PlotMultiCanvasBkgVsSigSpectra("{0}_BkgVsSig".format(s.fName), s.fSignalHistograms, "Sig. Window", s.fSideBandHistograms, "SB Window")

        self.CompareFeedDown(s)

    def CompareFeedDown(self, s):
        if not s.fFDHistogram: return
        before = s.fHistogram.Clone("{0}_compareFD".format(s.fHistogram.GetName()))
        before.SetTitle("Before FD correction")
        globalList.append(before)

        after = s.fFDCorrHistogram.Clone("{0}_compareFD".format(s.fFDCorrHistogram.GetName()))
        after.SetTitle("After FD correction")
        globalList.append(after)

        fd = s.fFDHistogram.Clone("{0}_compareFD".format(s.fFDHistogram.GetName()))
        fd.SetTitle("FD correction")
        globalList.append(fd)

        cname = "{0}_FDCorrection".format(s.fName)
        r = DMesonJetUtils.CompareSpectra(before, [after, fd], cname)
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)

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
        c = ROOT.TCanvas("{0}_canvas".format(s.fNormHistogram.GetName()), s.fNormHistogram.GetName())
        c.SetRightMargin(0.18)
        c.SetLogz()
        self.fCanvases.append(c)
        c.cd()
        h = s.fNormHistogram.DrawCopy("colz")
        h.GetZaxis().SetTitleOffset(1.4)

        globalList.append(c)
        globalList.append(h)

        c = ROOT.TCanvas("{0}_canvas".format(s.fUncertainty.GetName()), s.fUncertainty.GetName())
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
                    s.GenerateFDCorrectedSpectrum(self.fEvents, self.fIsWeighted)
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
                s.fBackground.SetBinContent(xbin, bin.fMassFitter.GetBackground())
                s.fBackground.SetBinError(xbin, bin.fMassFitter.GetBackgroundError())
            if s.fMass:
                s.fMass.SetBinContent(xbin, bin.fMassFitter.GetSignalMean())
                s.fMass.SetBinError(xbin, bin.fMassFitter.GetSignalMeanError())
            if s.fMassWidth:
                s.fMassWidth.SetBinContent(xbin, bin.fMassFitter.GetSignalWidth())
                s.fMassWidth.SetBinError(xbin, bin.fMassFitter.GetSignalWidthError())

    def GenerateSpectrum1DLikeSignMethod(self, s):
        binSetName = s.fBinSet.fBinSetName

        eng_LS = None
        for engTest in self.fEngines:
            if engTest.fTrigger == self.fTrigger and engTest.fDMeson == s.fLikeSignTree:
                eng_LS = engTest
                break
        if not eng_LS:
            print("Could not find engine with trigger '{0}' and tree name '{1}'".format(self.fTrigger, s.fLikeSignTree))
            return

        s.fLikeSignSubtractedBinSet = copy.deepcopy(s.fBinSet)
        s.fLikeSignSubtractedBinSet.fName = "{0}_LikeSignSubtracted_{1}".format(s.fLikeSignSubtractedBinSet.fName, s.fName)
        s.fLikeSignNormalizedBinSet = copy.deepcopy(eng_LS.fBinMultiSets[s.fJetType, s.fJetRadius].fBinSets[binSetName])
        s.fLikeSignNormalizedBinSet.fName = "{0}_Normalized_{1}".format(s.fLikeSignNormalizedBinSet.fName, s.fName)

        for ibin, (LS_sub_bin, LSbin, bin) in enumerate(zip(s.fLikeSignSubtractedBinSet.fBins, s.fLikeSignNormalizedBinSet.fBins, s.fBinSet.fBins)):
            if not bin.fInvMassHisto:
                continue
            if s.fSkipBins and ibin in s.fSkipBins:
                continue
            # Calculate the projections in the peak area for L-S and U-S
            if bin.fMassFitter and bin.fMassFitter.FitSuccessfull():
                sigma = bin.fMassFitter.GetSignalWidth()
                mean = bin.fMassFitter.GetSignalMean()
            else:
                sigma = s.fBackupSigma
                mean = s.fBackupMean
            binSig_1 = bin.fInvMassHisto.GetXaxis().FindBin(mean - s.fBinCountSignalSigmas * sigma)
            binSig_2 = bin.fInvMassHisto.GetXaxis().FindBin(mean + s.fBinCountSignalSigmas * sigma)

            if s.fAnalysisType == AnalysisType.LikeSign:
                if s.fAxis[0].fName == bin.fBinCountAnalysisAxis.fName:
                    sig = bin.fBinCountAnalysisHisto.ProjectionY("{0}_UnlikeSign_{1}".format(s.fName, bin.GetName()), binSig_1, binSig_2, "e")
                    LStotal = LSbin.fBinCountAnalysisHisto.ProjectionY("{0}_LikeSign_{1}".format(s.fName, bin.GetName()), binSig_1, binSig_2, "e")
                else:
                    sig = self.BuildSpectrum1D(s, "{0}_UnlikeSign_{1}".format(s.fName, bin.GetName()), "counts")
                    LStotal = self.BuildSpectrum1D(s, "{0}_LikeSign_{1}".format(s.fName, bin.GetName()), "counts")
                    ibin = sig.GetXaxis().FindBin(bin.GetBinCenter(s.fAxis[0].fName))
                    sigValueErr = ROOT.Double(0.)
                    sigValue = bin.fInvMassHisto.IntegralAndError(binSig_1, binSig_2, sigValueErr)
                    sig.SetBinContent(ibin, sigValue)
                    sig.SetBinError(ibin, sigValueErr)
                    lsValueErr = ROOT.Double(0.)
                    lsValue = LSbin.fInvMassHisto.IntegralAndError(binSig_1, binSig_2, lsValueErr)
                    LStotal.SetBinContent(ibin, lsValue)
                    LStotal.SetBinError(ibin, lsValueErr)
            else:
                sig = None
                LStotal = None

            print("Bin: {0}".format(bin.GetTitle()))
            if sig:
                print("The total signal+background is {0}".format(sig.Integral(0, -1)))
            print("The signal+background from the invariant mass plot {0} or summing signal and background {1}".format(bin.fInvMassHisto.Integral(binSig_1, binSig_2), bin.fMassFitter.GetSignal() + bin.fMassFitter.GetBackground(s.fBinCountSignalSigmas)))
            peakAreaBkgError = ROOT.Double(0.)
            peakAreaBkg = bin.fMassFitter.GetBackgroundAndError(peakAreaBkgError, s.fBinCountSignalSigmas)
            print("The total background in the peak area estimated from the fit is {0} +/- {1}".format(peakAreaBkg, peakAreaBkgError))
            bkgNorm = 0
            bkgNorm_err = 0
            # Two methods to normalize the background
            # Method 1: use the side-bands
            if s.fSideBandMaxSigmas > s.fSideBandMinSigmas:
                binSBL_1 = bin.fInvMassHisto.GetXaxis().FindBin(mean - s.fSideBandMaxSigmas * sigma)
                binSBL_2 = bin.fInvMassHisto.GetXaxis().FindBin(mean - s.fSideBandMinSigmas * sigma)
                binSBR_1 = bin.fInvMassHisto.GetXaxis().FindBin(mean + s.fSideBandMinSigmas * sigma)
                binSBR_2 = bin.fInvMassHisto.GetXaxis().FindBin(mean + s.fSideBandMaxSigmas * sigma)
                if binSBL_1 < 1:
                    binSBL_1 = 1
                if binSBR_2 > bin.fInvMassHisto.GetXaxis().GetNbins():
                    binSBR_2 = bin.fInvMassHisto.GetXaxis().GetNbins()

                SB_L_err = ROOT.Double(0)
                SB_L = bin.fInvMassHisto.IntegralAndError(binSBL_1, binSBL_2, SB_L_err)
                SB_R_err = ROOT.Double(0)
                SB_R = bin.fInvMassHisto.IntegralAndError(binSBR_1, binSBR_2, SB_R_err)
                SB_total = SB_L + SB_R
                SB_total_err = math.sqrt(SB_L_err ** 2 + SB_R_err ** 2)

                LS_SB_L_err = ROOT.Double(0)
                LS_SB_L = LSbin.fInvMassHisto.IntegralAndError(binSBL_1, binSBL_2, SB_L_err)
                LS_SB_R_err = ROOT.Double(0)
                LS_SB_R = LSbin.fInvMassHisto.IntegralAndError(binSBR_1, binSBR_2, SB_R_err)
                LS_SB_total = LS_SB_L + LS_SB_R
                LS_SB_total_err = math.sqrt(LS_SB_L_err ** 2 + LS_SB_R_err ** 2)

                if LS_SB_total > 0:
                    bkgNorm = SB_total / LS_SB_total
                    bkgNorm_err = ((SB_total_err / SB_total) ** 2 + (LS_SB_total_err / LS_SB_total) ** 2) * bkgNorm
                else:
                    bkgNorm = 1
                    bkgNorm_err = 0

                print("The background in side bands U-S is: {0} + {1} = {2}".format(SB_L, SB_R, SB_total))
                print("The background in side bands L-S is: {0} + {1} = {2}".format(LS_SB_L, SB_R, LS_SB_total))
                print("The background normalization is {0} +/- {1}".format(bkgNorm, bkgNorm_err))

                if LStotal:
                    for xbin in range(0, LStotal.GetNbinsX() + 2):
                        if LStotal.GetBinContent(xbin) == 0:
                            continue
                        error = math.sqrt(LStotal.GetBinError(xbin) ** 2 / LStotal.GetBinContent(xbin) ** 2 + bkgNorm_err ** 2 / bkgNorm ** 2) * LStotal.GetBinContent(xbin) * bkgNorm
                        cont = LStotal.GetBinContent(xbin) * bkgNorm
                        LStotal.SetBinError(xbin, error)
                        LStotal.SetBinContent(xbin, cont)

                (signalWindowInvMassHisto, sideBandWindowInvMassHisto) = self.GenerateInvMassWidonws(bin.fInvMassHisto, binSBL_1, binSBL_2, binSBR_1, binSBR_2, binSig_1, binSig_2)
                s.fSideBandWindowInvMassHistos[sideBandWindowInvMassHisto.GetName()] = sideBandWindowInvMassHisto
                s.fSignalWindowInvMassHistos[signalWindowInvMassHisto.GetName()] = signalWindowInvMassHisto

            # Method 2: use the peak area
            elif s.fBinCountNormSignalSigmas > 0:
                if not bin.fMassFitter.FitSuccessfull():
                    continue
                LStotalIntegralError = ROOT.Double(0)
                LStotalIntegral = bin.fInvMassHisto.IntegralAndError(binSig_1, binSig_2, LStotalIntegralError)
                if LStotalIntegral > 0:
                    if LStotal:
                        for xbin in range(0, LStotal.GetNbinsX() + 2):
                            if LStotal.GetBinContent(xbin) == 0:
                                continue
                            # Error propagation
                            error2_1 = peakAreaBkgError ** 2 * LStotal.GetBinContent(xbin) ** 2
                            error2_2 = peakAreaBkg ** 2 / LStotalIntegral ** 2 * LStotal.GetBinError(xbin) ** 2 * (LStotalIntegral - LStotal.GetBinContent(xbin)) ** 2
                            error2_3 = peakAreaBkg ** 2 / LStotalIntegral ** 2 * LStotal.GetBinContent(xbin) ** 2 * (LStotalIntegralError ** 2 - LStotal.GetBinError(xbin) ** 2)
                            # print("bin {0}: error1 = {1}, error2 = {2}, error3 = {3}".format(xbin, math.sqrt(error2_1)/LStotalIntegral, math.sqrt(error2_2)/LStotalIntegral, math.sqrt(error2_3)/LStotalIntegral))
                            error = math.sqrt(error2_1 + error2_2 + error2_3) / LStotalIntegral
                            cont = LStotal.GetBinContent(xbin) * peakAreaBkg / LStotalIntegral
                            LStotal.SetBinError(xbin, error)
                            LStotal.SetBinContent(xbin, cont)
                    bkgNorm = peakAreaBkg / LStotalIntegral
                    # this is only a rough estimate
                    bkgNorm_err = ((peakAreaBkgError / peakAreaBkg) ** 2 + (LStotalIntegralError / LStotalIntegral) ** 2) * bkgNorm
                else:
                    bkgNorm = 1
                    bkgNorm_err = 0
            # no background normalization
            else:
                bkgNorm = 1
                bkgNorm_err = 0

            print("The background normalization is {0} +/- {1}".format(bkgNorm, bkgNorm_err))

            if LStotal:
                integralError = ROOT.Double(0.)
                integral = LStotal.IntegralAndError(0, -1, integralError)
                print("The total normalized like-sign background in the peak area is {0} +/- {1}".format(integral, integralError))

            LSbin.fInvMassHisto.Scale(bkgNorm)

            w = s.fEfficiencyWeight.GetEfficiencyWeightTH1ForPt(bin.GetBinCenter("d_pt"))
            if w != 1:
                if LStotal:
                    LStotal.Scale(w)
                if sig:
                    sig.Scale(w)

            if LStotal:
                LStotal.SetTitle(bin.GetTitle())
                s.fLikeSignHistograms.append(LStotal)
                s.fLikeSignTotalHistogram.Add(LStotal)

            if sig:
                sig.SetTitle(bin.GetTitle())
                s.fUnlikeSignHistograms.append(sig)
                s.fUnlikeSignTotalHistogram.Add(sig)

            LS_sub_bin.fInvMassHisto.Add(LSbin.fInvMassHisto, -1)
            LS_sub_bin.fMassFitter = None

        self.PlotInvMassPlotsBinSet(s.fLikeSignNormalizedBinSet.fName, s.fBinSet.fBins, s.fLikeSignNormalizedBinSet.fBins, s)

        if s.fAnalysisType == AnalysisType.LikeSign:
            s.fHistogram.Add(s.fUnlikeSignTotalHistogram)
            s.fHistogram.Add(s.fLikeSignTotalHistogram, -1)
            if s.fBackground:
                s.fBackground.Add(s.fLikeSignTotalHistogram)

            for xbin in range(0, s.fHistogram.GetNbinsX() + 2):
                if s.fHistogram.GetBinContent(xbin) > 0:
                    s.fUncertainty.SetBinContent(xbin, s.fHistogram.GetBinError(xbin) / s.fHistogram.GetBinContent(xbin))
                else:
                    s.fUncertainty.SetBinContent(xbin, 0)
        elif s.fAnalysisType == AnalysisType.LikeSignFit:
            self.FitInvMassPlotsBinSet(s.fLikeSignSubtractedBinSet.fName, s.fLikeSignSubtractedBinSet.fBins, s.fLikeSignSubtractedBinSet.fFitOptions, 0.7)
            self.PlotInvMassPlotsBinSet(s.fLikeSignSubtractedBinSet.fName, s.fLikeSignSubtractedBinSet.fBins)
            self.GenerateSpectrum1DInvMassFit(s)

    def GenerateInvMassWidonws(self, invMassHisto, binSBL_1, binSBL_2, binSBR_1, binSBR_2, binSig_1, binSig_2):
        print("Generating invariant mass windows for {0}".format(invMassHisto.GetName()))
        sideBandWindowInvMassHisto = invMassHisto.Clone(invMassHisto.GetName().replace("InvMass", "InvMassSBWindow"))
        signalWindowInvMassHisto = invMassHisto.Clone(invMassHisto.GetName().replace("InvMass", "InvMassSigWindow"))
        for xbin in range(0, invMassHisto.GetNbinsX() + 2):
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

            sigma = bin.fMassFitter.GetSignalWidth()
            mean = bin.fMassFitter.GetSignalMean()

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

            (signalWindowInvMassHisto, sideBandWindowInvMassHisto) = self.GenerateInvMassWidonws(bin.fInvMassHisto, binSBL_1, binSBL_2, binSBR_1, binSBR_2, binSig_1, binSig_2)

            s.fSideBandWindowInvMassHistos[sideBandWindowInvMassHisto.GetName()] = sideBandWindowInvMassHisto
            s.fSignalWindowInvMassHistos[signalWindowInvMassHisto.GetName()] = signalWindowInvMassHisto

            sbL = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SideBandWindowL_{1}".format(s.fName, bin.GetName()), binSBL_1, binSBL_2, "e")
            sbR = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SideBandWindowR_{1}".format(s.fName, bin.GetName()), binSBR_1, binSBR_2, "e")
            sig = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SignalWindow_{1}".format(s.fName, bin.GetName()), binSig_1, binSig_2, "e")
            sbTotal = sbL.Clone("{0}_SideBandWindow_{1}".format(s.fName, bin.GetName()))
            sbTotal.Add(sbR)

            peakAreaBkgError = ROOT.Double(0.)
            peakAreaBkg = bin.fMassFitter.GetBackgroundAndError(peakAreaBkgError, s.fBinCountSignalSigmas)
            print("Bin: {0}".format(bin.GetTitle()))
            print("The background in side bands is: {0} + {1} = {2}".format(sbL.Integral(0, -1), sbR.Integral(0, -1), sbTotal.Integral(0, -1)))
            print("The estimated background in the signal window is {0} +/- {1}".format(peakAreaBkg, peakAreaBkgError))
            print("The total signal+background is {0}, which is the same from the invariant mass plot {1} or summing signal and background {2}".format(sig.Integral(0, -1), bin.fInvMassHisto.Integral(binSig_1, binSig_2), bin.fMassFitter.GetSignal() + bin.fMassFitter.GetBackground(s.fBinCountSignalSigmas)))

            sbTotalIntegralError = ROOT.Double(0)
            sbTotalIntegral = sbTotal.IntegralAndError(0, -1, sbTotalIntegralError)
            if sbTotalIntegral > 0:
                for xbin in range(0, sbTotal.GetNbinsX() + 2):
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

            integralError = ROOT.Double(0.)
            integral = sbTotal.IntegralAndError(0, -1, integralError)
            print("The total normalized side-band background is {0} +/- {1}".format(integral, integralError))

            w = s.fEfficiencyWeight.GetEfficiencyWeightTH1ForPt(bin.GetBinCenter("d_pt"))
            if w != 1:
                sbTotal.Scale(w)
                sig.Scale(w)

            sbTotal.SetTitle(bin.GetTitle())
            s.fSideBandHistograms.append(sbTotal)
            s.fSideBandWindowTotalHistogram.Add(sbTotal)

            sig.SetTitle(bin.GetTitle())
            s.fSignalHistograms.append(sig)
            s.fSignalWindowTotalHistogram.Add(sig)

        SBinvMassName = "{0}_SideBand_{1}".format(s.fBinSet.fName, s.fName)
        self.PlotInvMassPlotsBinSet(SBinvMassName, s.fBinSet.fBins, None, s)

        s.fHistogram.Add(s.fSignalWindowTotalHistogram)
        s.fHistogram.Add(s.fSideBandWindowTotalHistogram, -1)
        if s.fBackground:
            s.fBackground.Add(s.fSideBandWindowTotalHistogram)

        for xbin in range(0, s.fHistogram.GetNbinsX() + 2):
            if s.fHistogram.GetBinContent(xbin) > 0:
                s.fUncertainty.SetBinContent(xbin, s.fHistogram.GetBinError(xbin) / s.fHistogram.GetBinContent(xbin))
            else:
                s.fUncertainty.SetBinContent(xbin, 0)

    def GenerateSpectrum1DTruth(self, s):
        # The truth spectrum is already done, only need to apply the efficiency
        for ibin in range(0, s.fHistogram.GetNbinsX() + 2):
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
            self.GenerateSpectrum1DLikeSignMethod(s)
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
        for xbin in range(0, s.fHistogram.GetNbinsX() + 2):
            if s.fAxis[0].fName == "d_pt":
                w = s.fEfficiencyWeight.GetEfficiencyWeightTH1ForPt(s.fHistogram.GetXaxis().GetBinCenter(xbin))
            for ybin in range(0, s.fHistogram.GetNbinsY() + 2):
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
                self.FitInvMassPlotsBinSet(binSet.fName, binSet.fBins, binSet.fFitOptions)

    def FitInvMassPlotsBinSet(self, name, bins, fitOptions, initialSigOverBkg=0.1):
        print("Fitting {0}".format(name))
        pdgMass = ROOT.TDatabasePDG.Instance().GetParticle(421).Mass()
        massWidth = 0.015

        for i, bin in enumerate(bins):
            if not bin.fInvMassHisto:
                continue
            if self.fTrigger:
                fitterName = "InvMass_{0}_{1}_{2}_fitter".format(self.fTrigger, self.fDMeson, bin.GetName())
            else:
                fitterName = "InvMass_{0}_{1}_fitter".format(self.fDMeson, bin.GetName())

            fitter = self.CreateMassFitter(fitterName)
            bin.SetMassFitter(fitter)
            integral = bin.fInvMassHisto.Integral(1, bin.fInvMassHisto.GetXaxis().GetNbins())
            if integral > 10:
                fitter.GetFitFunction().SetParameter(0, integral)  # total integral
                fitter.GetFitFunction().SetParLimits(0, integral - 3 * math.sqrt(integral), integral + 3 * math.sqrt(integral))  # total integral
                if initialSigOverBkg <= 1 and initialSigOverBkg >= 0:
                    fitter.GetFitFunction().SetParameter(2, integral * initialSigOverBkg)  # signal integral
                    fitter.GetFitFunction().SetParLimits(2, 0, integral + 3 * math.sqrt(integral))  # signal integral has to be contained in the total integral
                else:
                    fitter.GetFitFunction().SetParameter(2, integral)  # signal integral
                    fitter.GetFitFunction().SetParLimits(2, integral * (initialSigOverBkg - 1), integral + 3 * math.sqrt(integral))  # signal integral has to be contained in the total integral
            else:
                fitter.GetFitFunction().SetParameter(0, 10)  # total integral
                fitter.GetFitFunction().SetParameter(2, 10)  # signal integral

            fitter.GetFitFunction().SetParameter(3, pdgMass)  # start fitting using PDG mass
            fitter.GetFitFunction().SetParLimits(3, pdgMass * 0.9975, pdgMass * 1.0025)  # start fitting using PDG mass
            fitter.GetFitFunction().SetParameter(4, massWidth)  # start fitting using mass peak width =
            print("Fitting bin {0}".format(bin.GetTitle()))

            fitter.Fit(bin.fInvMassHisto, fitOptions)

    def PlotInvMassPlots(self):
        for binMultiSet in self.fBinMultiSets.itervalues():
            for binSet in binMultiSet.fBinSets.itervalues():
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
                results = DMesonJetUtils.CompareSpectra(spectraToCompare[0], spectraToCompare[1:], cname)
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
