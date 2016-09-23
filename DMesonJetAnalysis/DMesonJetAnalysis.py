#!/usr/bin/env python
#python program to perform a D meson jet analysis

import ROOT
import math
import DMesonJetProjectors
from DMesonJetBase import *
import array
import copy
import DMesonJetUtils
import collections

globalList = []

class DMesonJetAnalysisEngine:
    def __init__(self, figTitle, collision, trigger, dmeson, binSet, nMassBins, minMass, maxMass, jets, spectra, projector):
        self.fFigureTitle = figTitle
        self.fCollision = collision
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fBinMultiSet = binSet
        self.fNMassBins = nMassBins
        self.fMinMass = minMass
        self.fMaxMass = maxMass
        self.fSpectra = collections.OrderedDict()
        self.fCanvases = []
        for s in spectra:
            if "active" in s and not s["active"]:
                continue
            name = "{0}_{1}".format(self.fDMeson, s["name"])
            self.fSpectra[s["name"]] = Spectrum(s, name, self.fBinMultiSet.fBinSets)
            
    def CompareSpectra(self):
        spectraToCompare = []
        axisBaseline = None
        for name,s in self.fSpectra.iteritems():
            if not s.fNormHistogram:
                continue
            if len(s.fAxis) != 1:
                continue
            if axisBaseline:
                if axisBaseline.fName != s.fAxis[0].fName:
                    continue
            else:
                axisBaseline = s.fAxis[0]
            h = s.fNormHistogram.Clone("{0}_copy".format(s.fNormHistogram.GetName()))
            if s.fTitle:
                h.SetTitle(s.fTitle)
            globalList.append(h)
            spectraToCompare.append(h)
        if len(spectraToCompare) < 2:
            return
        results = DMesonJetUtils.CompareSpectra(spectraToCompare[0], spectraToCompare[1:], "{0}_SpectraComparison".format(self.fDMeson), "", "hist")
        for obj in results:
            if obj and isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)
            globalList.append(obj)        

    def SaveRootFile(self, file):
        file.cd()
        for rlist in self.fBinMultiSet.GenerateInvMassRootLists():
            rlist.Write("{0}_{1}".format(self.fDMeson, rlist.GetName()), ROOT.TObject.kSingleKey)
        for s in self.fSpectra.itervalues():
            if s.fNormHistogram:
                s.fNormHistogram.Write()
            slist = s.GenerateRootList()
            if slist.GetEntries() > 0:
                slist.Write(slist.GetName(), ROOT.TObject.kSingleKey)

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
            massFitTypeSig = ROOT.MassFitter.kGaus
            massFitTypeBkg = ROOT.MassFitter.kExpo
        elif self.fDMeson == "DStar":
            print("Not implemented for DStar!")

        fitter = ROOT.MassFitter(name, massFitTypeSig, massFitTypeBkg, minMass, maxMass)
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
        self.fProjector.GetInvMassHisograms(self.fTrigger, self.fDMeson, self.fJetDefinitions, 
                                            self.fBinMultiSet, self.fNMassBins, self.fMinMass, self.fMaxMass)

        self.fEvents = self.fProjector.fTotalEvents

    def Start(self, engines):
        self.fEngines = engines
        if not "MCTruth" in self.fDMeson:
            self.FitInvMassPlots()
        if not "BackgroundOnly" in self.fDMeson:
            self.GenerateSpectra()

        if not "MCTruth" in self.fDMeson:
            self.PlotInvMassPlots()            
        if not "BackgroundOnly" in self.fDMeson:
            self.PlotSpectra()

    def PlotSpectra(self):
        for s in self.fSpectra.itervalues():
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

        h.SetMarkerColor(ROOT.kBlue+2)
        h.SetMarkerStyle(ROOT.kFullCircle)
        h.SetMarkerSize(0.9)
        h.SetLineColor(ROOT.kBlue+2)

        pave = ROOT.TPaveText(0.60, 0.88, 0.9, 0.55, "NB NDC")
        pave.SetFillStyle(0)
        pave.SetBorderSize(0)
        pave.SetTextFont(43)
        pave.SetTextSize(15)
        pave.AddText(self.fFigureTitle)
        pave.AddText(self.fCollision)
        pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4")
        pave.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
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

        h.SetLineColor(ROOT.kBlue+2)
        h.SetLineWidth(2)
        h.SetFillColorAlpha(ROOT.kBlue+2, 0.25)
        h.SetFillStyle(1001)

        pave = ROOT.TPaveText(0.10, 0.88, 0.4, 0.55, "NB NDC")
        pave.SetFillStyle(0)
        pave.SetBorderSize(0)
        pave.SetTextFont(43)
        pave.SetTextSize(15)
        pave.AddText(self.fFigureTitle)
        pave.AddText(self.fCollision)
        pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4")
        pave.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
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
    
            h.SetMarkerColor(ROOT.kBlue+2)
            h.SetMarkerStyle(ROOT.kFullCircle)
            h.SetMarkerSize(0.9)
            h.SetLineColor(ROOT.kBlue+2)
            h.GetYaxis().SetRangeUser(1.84, 1.91)
    
            pave = ROOT.TPaveText(0.10, 0.88, 0.8, 0.68, "NB NDC")
            pave.SetFillStyle(0)
            pave.SetBorderSize(0)
            pave.SetTextFont(43)
            pave.SetTextSize(15)
            pave.AddText("{0} {1}".format(self.fFigureTitle, self.fCollision))
            pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4 with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
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
    
            h.SetMarkerColor(ROOT.kBlue+2)
            h.SetMarkerStyle(ROOT.kFullCircle)
            h.SetMarkerSize(0.9)
            h.SetLineColor(ROOT.kBlue+2)
            #h.GetYaxis().SetRangeUser(1.84, 1.91)
    
            pave = ROOT.TPaveText(0.10, 0.88, 0.8, 0.68, "NB NDC")
            pave.SetFillStyle(0)
            pave.SetBorderSize(0)
            pave.SetTextFont(43)
            pave.SetTextSize(15)
            pave.AddText("{0} {1}".format(self.fFigureTitle, self.fCollision))
            pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4 with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
            pave.AddText(s.fTitle)
            pave.Draw()
    
            globalList.append(c)
            globalList.append(h)
            globalList.append(pave)
            globalList.append(line)

        if s.fLikeSignTotalHistogram:
            self.PlotBackgroundVsSignalSpectra("{0}_BkgVsSig".format(s.fUnlikeSignTotalHistogram.GetName()), s.fUnlikeSignTotalHistogram, "Unlike Sign", s.fLikeSignTotalHistogram, "Like Sign")
            for hSig, hBkg in zip(s.fUnlikeSignHistograms, s.fLikeSignHistograms):
                self.PlotBackgroundVsSignalSpectra("{0}_BkgVsSig".format(hSig.GetName()), hSig, "Unlike Sign", hBkg, "Like Sign")

        if s.fSideBandWindowTotalHistogram:
            self.PlotBackgroundVsSignalSpectra("{0}_BkgVsSig".format(s.fSignalWindowTotalHistogram.GetName()), s.fSideBandWindowTotalHistogram, "Signal Window", s.fSignalWindowTotalHistogram, "Side Band Window")
            for hSig, hBkg in zip(s.fSignalHistograms, s.fSideBandHistograms):
                self.PlotBackgroundVsSignalSpectra("{0}_BkgVsSig".format(hSig.GetName()), hSig, "Signal Window", hBkg, "Side Band Window")

    def PlotBackgroundVsSignalSpectra(self, cname, hSigOrig, hTitleSig, hBkgOrig, hTitleBkg):
        c = ROOT.TCanvas(cname, cname, 700, 700)
        c.SetTopMargin(0.05)
        self.fCanvases.append(c)

        hSig = hSigOrig.DrawCopy()
        hSig.SetMarkerColor(ROOT.kBlue+2)
        hSig.SetMarkerStyle(ROOT.kOpenCircle)
        hSig.SetMarkerSize(0.9)
        hSig.SetLineColor(ROOT.kBlue+2)
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
        hSig.SetMaximum(hSig.GetMaximum()*1.3)
        hSig.SetMinimum(0)

        hBkg = hBkgOrig.DrawCopy("same")
        hBkg.SetMarkerColor(ROOT.kRed+2)
        hBkg.SetMarkerStyle(ROOT.kOpenSquare)
        hBkg.SetMarkerSize(0.9)
        hBkg.SetLineColor(ROOT.kRed+2)

        hSub = hSigOrig.DrawCopy("same")
        hSub.Add(hBkgOrig, -1)
        hSub.SetMarkerColor(ROOT.kGreen+2)
        hSub.SetMarkerStyle(ROOT.kStar)
        hSub.SetMarkerSize(0.9)
        hSub.SetLineColor(ROOT.kGreen+2)

        leg = ROOT.TLegend(0.12, 0.84, 0.57, 0.93, "", "NB NDC")
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
        globalList.append(c)
        globalList.append(hSig)
        globalList.append(hBkg)
        globalList.append(hSub)

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

    def PlotSpectrum3D(self, s):
        print("PlotSpectrum3D not implemented!")

    def GenerateSpectra(self):
        for s in self.fSpectra.itervalues():
            if len(s.fAxis) == 1:
                self.GenerateSpectrum1D(s)
            elif len(s.fAxis) == 2:
                self.GenerateSpectrum2D(s)
            elif len(s.fAxis) == 3:
                self.GenerateSpectrum3D(s)
            else:
                print("Not able to generate spectra with dim > 3!")
            s.GenerateNormalizedSpectrum(self.fEvents)

    def BuildSpectrum1D(self, s, name, yaxis):
        hist = ROOT.TH1D(name, name, len(s.fAxis[0].fBins)-1, array.array('d',s.fAxis[0].fBins))
        hist.GetXaxis().SetTitle(s.fAxis[0].GetTitle())
        hist.GetYaxis().SetTitle(yaxis)
        hist.Sumw2()
        return hist

    def GenerateSpectrum1DInvMassFit(self, s):
        s.fHistogram = self.BuildSpectrum1D(s, s.fName, "counts")
        s.fUncertainty = self.BuildSpectrum1D(s, "{0}_Unc".format(s.fName), "relative statistical uncertainty")
        s.fMass = self.BuildSpectrum1D(s, "{0}_Mass".format(s.fName), "D^{0} mass (GeV/#it{c}^{2})")
        s.fMassWidth = self.BuildSpectrum1D(s, "{0}_MassWidth".format(s.fName), "D^{0} mass width (GeV/#it{c}^{2})")
        s.fBackground = self.BuildSpectrum1D(s, "{0}_Bkg".format(s.fName), "background |#it{m} - <#it{m}>| < 2#sigma")
        for binSetName in s.fBins:
            for bin in self.fBinMultiSet.fBinSets[binSetName].fBins:
                if bin.fMassFitter is None:
                    print("The bin printed below does not have a mass fitter!")
                    bin.Print()
                    continue
                if self.fBinMultiSet.fBinSets[binSetName].fApplyEfficiencyToSpectrum:
                    w = self.fBinMultiSet.fBinSets[binSetName].fWeightEfficiency.GetEfficiencyWeightTH1ForPt(bin.GetBinCenter("d_pt"))
                else:
                    w = 1
                xbin = s.fHistogram.GetXaxis().FindBin(bin.GetBinCenter(s.fAxis[0].fName))
                if "SignalOnly" in self.fDMeson:
                    signal = bin.fInvMassHisto.Integral()*w
                    signal_unc = math.sqrt(bin.fInvMassHisto.Integral())*w
                else:
                    signal = bin.fMassFitter.GetSignal()*w
                    signal_unc = bin.fMassFitter.GetSignalError()*w

                s.fHistogram.SetBinContent(xbin, signal)
                s.fHistogram.SetBinError(xbin, signal_unc)
                s.fBackground.SetBinContent(xbin, bin.fMassFitter.GetBackground()*w)
                s.fBackground.SetBinError(xbin, bin.fMassFitter.GetBackgroundError()*w)
                s.fUncertainty.SetBinContent(xbin, signal_unc/signal) 
                s.fMass.SetBinContent(xbin, bin.fMassFitter.GetSignalMean())
                s.fMass.SetBinError(xbin, bin.fMassFitter.GetSignalMeanError())
                s.fMassWidth.SetBinContent(xbin, bin.fMassFitter.GetSignalWidth())
                s.fMassWidth.SetBinError(xbin, bin.fMassFitter.GetSignalWidthError())

    def GenerateSpectrum1DLikeSignMethod(self, s):
        s.fLikeSignHistograms = []
        s.fUnlikeSignHistograms = []
        s.fMass = None
        s.fMassWidth = None
        s.fHistogram = self.BuildSpectrum1D(s, s.fName, "counts")
        s.fUncertainty = self.BuildSpectrum1D(s, "{0}_Unc".format(s.fName), "relative statistical uncertainty")
        s.fBackground = self.BuildSpectrum1D(s, "{0}_Bkg".format(s.fName), "background |#it{{m}} - <#it{{m}}>| < {0}#sigma".format(int(s.fBinCountSignalSigmas)))
        s.fLikeSignTotalHistogram = self.BuildSpectrum1D(s, "{0}_LikeSignTotal".format(s.fName), "counts")
        s.fUnlikeSignTotalHistogram = self.BuildSpectrum1D(s, "{0}_UnlikeSignTotal".format(s.fName), "counts")
        print("Generating spectrum {0}".format(s.fName))
        binSetName = s.fBins[0]

        if "SignalOnly" in self.fDMeson:
            for bin in self.fBinMultiSet.fBinSets[binSetName].fBins:
                sig = bin.fBinCountAnalysisHisto.ProjectionY("{0}_UnlikeSign_{1}".format(s.fName, bin.GetName()), 0, -1, "e")
                if self.fBinMultiSet.fBinSets[binSetName].fApplyEfficiencyToSpectrum:
                    w = self.fBinMultiSet.fBinSets[binSetName].fWeightEfficiency.GetEfficiencyWeightTH1ForPt(bin.GetBinCenter("d_pt"))
                    sig.Scale(w)
                sig.SetTitle(bin.GetTitle())
                s.fUnlikeSignHistograms.append(sig)
                s.fUnlikeSignTotalHistogram.Add(sig)
        else:
            eng_LS = None
            for engTest in self.fEngines:
                if engTest.fTrigger == self.fTrigger and engTest.fDMeson == s.fLikeSignTree:
                    eng_LS = engTest
                    break
            if not eng_LS:
                print("Could not find engine with trigger '{0}' and tree name '{1}'".format(self.fTrigger, s.fLikeSignTree))
                return

            s.fLikeSignSubtractedBinSet = copy.deepcopy(self.fBinMultiSet.fBinSets[binSetName])
            s.fLikeSignSubtractedBinSet.fName = "{0}_LikeSignSubtracted_{1}".format(s.fLikeSignSubtractedBinSet.fName, s.fName)

            for (LS_sub_bin, LSbin, bin) in zip(s.fLikeSignSubtractedBinSet.fBins, eng_LS.fBinMultiSet.fBinSets[binSetName].fBins, self.fBinMultiSet.fBinSets[binSetName].fBins):
                # Calculate the projections in the peak area for L-S and U-S
                binSig_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - s.fBinCountSignalSigmas*bin.fMassFitter.GetSignalWidth())
                binSig_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + s.fBinCountSignalSigmas*bin.fMassFitter.GetSignalWidth())
                sig = bin.fBinCountAnalysisHisto.ProjectionY("{0}_UnlikeSign_{1}".format(s.fName, bin.GetName()), binSig_1, binSig_2, "e")
                LStotal = LSbin.fBinCountAnalysisHisto.ProjectionY("{0}_LikeSign_{1}".format(s.fName, bin.GetName()), binSig_1, binSig_2, "e")

                print("Bin: {0}".format(bin.GetTitle()))
                print("The total signal+background is {0}, which is the same from the invariant mass plot {1} or summing signal and background {2}".format(sig.Integral(0,-1), bin.fInvMassHisto.Integral(binSig_1, binSig_2), bin.fMassFitter.GetSignal()+bin.fMassFitter.GetBackground(s.fBinCountSignalSigmas)))
                peakAreaBkgError = ROOT.Double(0.)
                peakAreaBkg = bin.fMassFitter.GetBackgroundAndError(peakAreaBkgError, s.fBinCountSignalSigmas)
                print("The total background in the peak area estimated from the fit is {0} +/- {1}".format(peakAreaBkg, peakAreaBkgError))
                bkgNorm = 0
                bkgNorm_err = 0
                # Two methods to normalize the background
                # Method 1: use the side-bands
                if s.fSideBandMaxSigmas > s.fSideBandMinSigmas:
                    binSBL_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - s.fSideBandMaxSigmas*bin.fMassFitter.GetSignalWidth())
                    binSBL_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - s.fSideBandMinSigmas*bin.fMassFitter.GetSignalWidth())
                    binSBR_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + s.fSideBandMinSigmas*bin.fMassFitter.GetSignalWidth())
                    binSBR_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + s.fSideBandMaxSigmas*bin.fMassFitter.GetSignalWidth())
                    if binSBL_1 < 1:
                        binSBL_1 = 1
                    if binSBR_2 > bin.fBinCountAnalysisHisto.GetXaxis().GetNbins():
                        binSBR_2 = bin.fBinCountAnalysisHisto.GetXaxis().GetNbins()
    
                    SB_L_err = ROOT.Double(0)
                    SB_L = bin.fInvMassHisto.IntegralAndError(binSBL_1, binSBL_2, SB_L_err)
                    SB_R_err = ROOT.Double(0)
                    SB_R = bin.fInvMassHisto.IntegralAndError(binSBR_1, binSBR_2, SB_R_err)
                    SB_total = SB_L + SB_R
                    SB_total_err = math.sqrt(SB_L_err**2 + SB_R_err**2)
    
                    LS_SB_L_err = ROOT.Double(0)
                    LS_SB_L = LSbin.fInvMassHisto.IntegralAndError(binSBL_1, binSBL_2, SB_L_err)
                    LS_SB_R_err = ROOT.Double(0)
                    LS_SB_R = LSbin.fInvMassHisto.IntegralAndError(binSBR_1, binSBR_2, SB_R_err)
                    LS_SB_total = LS_SB_L + LS_SB_R
                    LS_SB_total_err = math.sqrt(LS_SB_L_err**2 + LS_SB_R_err**2)
    
                    if LS_SB_total > 0:
                        bkgNorm = SB_total / LS_SB_total
                        bkgNorm_err = ((SB_total_err/SB_total)**2 + (LS_SB_total_err/LS_SB_total)**2) * bkgNorm
                    else:
                        bkgNorm = 0
                        bkgNorm_err = 0
    
                    print("The background in side bands U-S is: {0} + {1} = {2}".format(SB_L, SB_R, SB_total))
                    print("The background in side bands L-S is: {0} + {1} = {2}".format(LS_SB_L, SB_R, LS_SB_total))
                    print("The background normalization is {0} +/- {1}".format(bkgNorm, bkgNorm_err))

                    for xbin in range(0, LStotal.GetNbinsX()+2):
                        if LStotal.GetBinContent(xbin) == 0:
                            continue
                        error = math.sqrt(LStotal.GetBinError(xbin)**2/LStotal.GetBinContent(xbin)**2+bkgNorm_err**2/bkgNorm**2)*LStotal.GetBinContent(xbin)*bkgNorm
                        cont = LStotal.GetBinContent(xbin)*bkgNorm
                        LStotal.SetBinError(xbin, error)
                        LStotal.SetBinContent(xbin, cont)
                # Method 2: use the peak area
                else:
                    LStotalIntegralError = ROOT.Double(0)
                    LStotalIntegral = LStotal.IntegralAndError(0, -1, LStotalIntegralError)
                    if LStotalIntegral > 0:
                        for xbin in range(0, LStotal.GetNbinsX()+2):
                            if LStotal.GetBinContent(xbin) == 0:
                                continue
                            # Error propagation
                            error2_1 = peakAreaBkgError**2 * LStotal.GetBinContent(xbin)**2
                            error2_2 = peakAreaBkg**2 / LStotalIntegral**2 * LStotal.GetBinError(xbin)**2 * (LStotalIntegral-LStotal.GetBinContent(xbin))**2
                            error2_3 = peakAreaBkg**2 / LStotalIntegral**2 * LStotal.GetBinContent(xbin)**2 * (LStotalIntegralError**2-LStotal.GetBinError(xbin)**2)
                            #print("bin {0}: error1 = {1}, error2 = {2}, error3 = {3}".format(xbin, math.sqrt(error2_1)/LStotalIntegral, math.sqrt(error2_2)/LStotalIntegral, math.sqrt(error2_3)/LStotalIntegral))
                            error = math.sqrt(error2_1 + error2_2 + error2_3) / LStotalIntegral
                            cont = LStotal.GetBinContent(xbin)*peakAreaBkg/LStotalIntegral
                            LStotal.SetBinError(xbin, error)
                            LStotal.SetBinContent(xbin, cont)
                        bkgNorm = peakAreaBkg / LStotalIntegral
                        # this is only a rough estimate
                        bkgNorm_err = ((peakAreaBkgError/peakAreaBkg)**2 + (LStotalIntegralError/LStotalIntegral)**2) * bkgNorm

                print("The background normalization is {0} +/- {1}".format(bkgNorm, bkgNorm_err))
                integralError = ROOT.Double(0.)
                integral = LStotal.IntegralAndError(0, -1, integralError)
                print("The total normalized like-sign background in the peak area is {0} +/- {1}".format(integral, integralError))

                LSbin.fInvMassHisto.Scale(bkgNorm)

                if self.fBinMultiSet.fBinSets[binSetName].fApplyEfficiencyToSpectrum:
                    w = self.fBinMultiSet.fBinSets[binSetName].fWeightEfficiency.GetEfficiencyWeightTH1ForPt(bin.GetBinCenter("d_pt"))
                    LStotal.Scale(w)
                    sig.Scale(w)

                LStotal.SetTitle(bin.GetTitle())
                s.fLikeSignHistograms.append(LStotal)
                s.fLikeSignTotalHistogram.Add(LStotal)

                sig.SetTitle(bin.GetTitle())
                s.fUnlikeSignHistograms.append(sig)
                s.fUnlikeSignTotalHistogram.Add(sig)
                
                LS_sub_bin.fInvMassHisto.Add(LSbin.fInvMassHisto, -1)
                LS_sub_bin.fMassFitter = None
                
            self.FitInvMassPlotsBinSet(s.fLikeSignSubtractedBinSet.fName, s.fLikeSignSubtractedBinSet.fBins, s.fLikeSignSubtractedBinSet.fFitOptions)
            self.PlotInvMassPlotsBinSet(s.fLikeSignSubtractedBinSet.fName, s.fLikeSignSubtractedBinSet.fBins)

        s.fHistogram.Add(s.fUnlikeSignTotalHistogram)
        s.fHistogram.Add(s.fLikeSignTotalHistogram, -1)
        s.fBackground.Add(s.fLikeSignTotalHistogram)

        for xbin in range(0, s.fHistogram.GetNbinsX()+2):
            if s.fHistogram.GetBinContent(xbin) > 0:
                s.fUncertainty.SetBinContent(xbin, s.fHistogram.GetBinError(xbin) / s.fHistogram.GetBinContent(xbin))
            else:
                s.fUncertainty.SetBinContent(xbin, 0)

    def GenerateSpectrum1DSideBandMethod(self, s):
        s.fSideBandHistograms = []
        s.fSignalHistograms = []
        s.fMass = None
        s.fMassWidth = None
        s.fHistogram = self.BuildSpectrum1D(s, s.fName, "counts")
        s.fUncertainty = self.BuildSpectrum1D(s, "{0}_Unc".format(s.fName), "relative statistical uncertainty")
        s.fBackground = self.BuildSpectrum1D(s, "{0}_Bkg".format(s.fName), "background |#it{{m}} - <#it{{m}}>| < {0}#sigma".format(int(s.fBinCountSignalSigmas)))
        s.fSideBandWindowTotalHistogram = self.BuildSpectrum1D(s, "{0}_SideBandWindowTotal".format(s.fName), "counts")
        s.fSignalWindowTotalHistogram = self.BuildSpectrum1D(s, "{0}_SignalWindowTotal".format(s.fName), "counts")
        print("Generating spectrum {0}".format(s.fName))
        for binSetName in s.fBins:
            for bin in self.fBinMultiSet.fBinSets[binSetName].fBins:
                if self.fBinMultiSet.fBinSets[binSetName].fApplyEfficiencyToSpectrum:
                    w = self.fBinMultiSet.fBinSets[binSetName].fWeightEfficiency.GetEfficiencyWeightTH1ForPt(bin.GetBinCenter("d_pt"))
                else:
                    w = 1
                if "SignalOnly" in self.fDMeson:
                    sig = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SignalWindow_{1}".format(s.fName, bin.GetName()), 0, -1, "e")
                    sig.SetTitle(bin.GetTitle())
                    if w != 1:
                        sig.Scale(w)
                    s.fSignalHistograms.append(sig)
                    s.fSignalWindowTotalHistogram.Add(sig)
                else:
                    binSBL_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - s.fSideBandMaxSigmas*bin.fMassFitter.GetSignalWidth())
                    binSBL_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - s.fSideBandMinSigmas*bin.fMassFitter.GetSignalWidth())
                    binSBR_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + s.fSideBandMinSigmas*bin.fMassFitter.GetSignalWidth())
                    binSBR_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + s.fSideBandMaxSigmas*bin.fMassFitter.GetSignalWidth())
                    binSig_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - s.fBinCountSignalSigmas*bin.fMassFitter.GetSignalWidth())
                    binSig_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + s.fBinCountSignalSigmas*bin.fMassFitter.GetSignalWidth())

                    if binSBL_1 < 1:
                        binSBL_1 = 1
                    if binSBR_2 > bin.fBinCountAnalysisHisto.GetXaxis().GetNbins():
                        binSBR_2 = bin.fBinCountAnalysisHisto.GetXaxis().GetNbins()
                    sbL = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SideBandWindowL_{1}".format(s.fName, bin.GetName()), binSBL_1, binSBL_2, "e")
                    sbR = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SideBandWindowR_{1}".format(s.fName, bin.GetName()), binSBR_1, binSBR_2, "e")
                    sig = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SignalWindow_{1}".format(s.fName, bin.GetName()), binSig_1, binSig_2, "e")
                    sbTotal = sbL.Clone("{0}_SideBandWindow_{1}".format(s.fName, bin.GetName()))
                    sbTotal.Add(sbR)

                    sideBandWindowInvMassHisto = bin.fInvMassHisto.Clone(bin.fInvMassHisto.GetName().replace("InvMass", "InvMassSBWindow"))
                    signalWindowInvMassHisto = bin.fInvMassHisto.Clone(bin.fInvMassHisto.GetName().replace("InvMass", "InvMassSigWindow"))

                    for xbin in range(0,bin.fInvMassHisto.GetNbinsX()+2):
                        if not ((xbin >= binSBL_1 and xbin < binSBL_2) or (xbin >= binSBR_1 and xbin < binSBR_2)):
                            sideBandWindowInvMassHisto.SetBinContent(xbin, 0)
                            sideBandWindowInvMassHisto.SetBinError(xbin, 0)
                        if not (xbin >= binSig_1 and xbin < binSig_2):
                            signalWindowInvMassHisto.SetBinContent(xbin, 0)
                            signalWindowInvMassHisto.SetBinError(xbin, 0)
                    s.fSideBandWindowInvMassHistos[sideBandWindowInvMassHisto.GetName()] = sideBandWindowInvMassHisto
                    s.fSignalWindowInvMassHistos[signalWindowInvMassHisto.GetName()] = signalWindowInvMassHisto

                    peakAreaBkgError = ROOT.Double(0.)
                    peakAreaBkg = bin.fMassFitter.GetBackgroundAndError(peakAreaBkgError, s.fBinCountSignalSigmas)
                    print("Bin: {0}".format(bin.GetTitle()))
                    print("The background in side bands is: {0} + {1} = {2}".format(sbL.Integral(0,-1), sbR.Integral(0,-1), sbTotal.Integral(0,-1)))
                    print("The estimated background in the signal window is {0} +/- {1}".format(peakAreaBkg, peakAreaBkgError))
                    print("The total signal+background is {0}, which is the same from the invariant mass plot {1} or summing signal and background {2}".format(sig.Integral(0,-1), bin.fInvMassHisto.Integral(binSig_1, binSig_2), bin.fMassFitter.GetSignal()+bin.fMassFitter.GetBackground(s.fBinCountSignalSigmas)))

                    sbTotalIntegralError = ROOT.Double(0)
                    sbTotalIntegral = sbTotal.IntegralAndError(0, -1, sbTotalIntegralError)
                    if sbTotalIntegral > 0:
                        for xbin in range(0, sbTotal.GetNbinsX()+2):
                            if sbTotal.GetBinContent(xbin) == 0:
                                continue
                            # Error propagation
                            error2_1 = peakAreaBkgError**2 * sbTotal.GetBinContent(xbin)**2
                            error2_2 = peakAreaBkg**2 / sbTotalIntegral**2 * sbTotal.GetBinError(xbin)**2 * (sbTotalIntegral-sbTotal.GetBinContent(xbin))**2
                            error2_3 = peakAreaBkg**2 / sbTotalIntegral**2 * sbTotal.GetBinContent(xbin)**2 * (sbTotalIntegralError**2-sbTotal.GetBinError(xbin)**2)
                            #print("bin {0}: error1 = {1}, error2 = {2}, error3 = {3}".format(xbin, math.sqrt(error2_1)/sbTotalIntegral, math.sqrt(error2_2)/sbTotalIntegral, math.sqrt(error2_3)/sbTotalIntegral))
                            error = math.sqrt(error2_1 + error2_2 + error2_3) / sbTotalIntegral
                            cont = sbTotal.GetBinContent(xbin)*peakAreaBkg/sbTotalIntegral
                            sbTotal.SetBinError(xbin, error)
                            sbTotal.SetBinContent(xbin, cont)

                    integralError = ROOT.Double(0.)
                    integral = sbTotal.IntegralAndError(0, -1, integralError)
                    print("The total normalized side-band background is {0} +/- {1}".format(integral, integralError))

                    if w != 1:
                        sig.Scale(w)
                        sbTotal.Scale(w)

                    sbTotal.SetTitle(bin.GetTitle())
                    s.fSideBandHistograms.append(sbTotal)
                    s.fSideBandWindowTotalHistogram.Add(sbTotal)

                    sig.SetTitle(bin.GetTitle())
                    s.fSignalHistograms.append(sig)
                    s.fSignalWindowTotalHistogram.Add(sig)

        s.fHistogram.Add(s.fSignalWindowTotalHistogram)
        s.fHistogram.Add(s.fSideBandWindowTotalHistogram, -1)
        s.fBackground.Add(s.fSideBandWindowTotalHistogram)

        for xbin in range(0, s.fHistogram.GetNbinsX()+2):
            if s.fHistogram.GetBinContent(xbin) > 0:
                s.fUncertainty.SetBinContent(xbin, s.fHistogram.GetBinError(xbin) / s.fHistogram.GetBinContent(xbin))
            else:
                s.fUncertainty.SetBinContent(xbin, 0)

    def GenerateSpectrum1D(self, s):
        if s.fAnalysisType == AnalysisType.SideBand:
            self.GenerateSpectrum1DSideBandMethod(s)
        elif s.fAnalysisType == AnalysisType.LikeSign:
            self.GenerateSpectrum1DLikeSignMethod(s)
        else:
            self.GenerateSpectrum1DInvMassFit(s)

    def BuildSpectrum2D(self, s, name, zaxis):
        hist = ROOT.TH2D(name, name, len(s.fAxis[0].fBins)-1, array.array('d',s.fAxis[0].fBins), len(s.fAxis[1].fBins)-1, array.array('d',s.fAxis[1].fBins))
        hist.GetXaxis().SetTitle(s.fAxis[0].GetTitle())
        hist.GetYaxis().SetTitle(s.fAxis[1].GetTitle())
        hist.GetZaxis().SetTitle(zaxis)
        hist.Sumw2()
        return hist

    def GenerateSpectrum2D(self, s):
        s.fHistogram = self.BuildSpectrum2D(s, s.fName, "counts")
        s.fUncertainty = self.BuildSpectrum2D(s, "{0}_Unc".format(s.fName), "relative statistical uncertainty")
        s.fMass = self.BuildSpectrum2D(s, "{0}_Mass".format(s.fName), "mass")
        s.fBackground = self.BuildSpectrum2D(s, "{0}_Bkg".format(s.fName), "background |#it{m} - <#it{m}>| < 3#sigma")
        for binSetName in s.fBins:
            for bin in self.fBinMultiSet.fBinSets[binSetName].fBins:
                if bin.fMassFitter is None:
                    print("The bin printed below does not have a mass fitter!")
                    bin.Print()
                    continue
                xbin = s.fHistogram.GetXaxis().FindBin(bin.GetBinCenter(s.fAxis[0].fName))
                ybin = s.fHistogram.GetYaxis().FindBin(bin.GetBinCenter(s.fAxis[1].fName))
                if "SignalOnly" in self.fDMeson:
                    signal = bin.fInvMassHisto.Integral()
                    signal_unc = math.sqrt(bin.fInvMassHisto.Integral())
                else:
                    signal = bin.fMassFitter.GetSignal()
                    signal_unc = bin.fMassFitter.GetSignalError()

                s.fHistogram.SetBinContent(xbin, ybin, signal)
                s.fHistogram.SetBinError(xbin, ybin, signal_unc)
                s.fBackground.SetBinContent(xbin, ybin, bin.fMassFitter.GetBackground())
                s.fBackground.SetBinError(xbin, ybin, bin.fMassFitter.GetBackgroundError())
                s.fUncertainty.SetBinContent(xbin, ybin, signal_unc/signal)
                s.fMass.SetBinContent(xbin, ybin, bin.fMassFitter.GetSignalMean())
                s.fMass.SetBinError(xbin, ybin, bin.fMassFitter.GetSignalMeanError())
                s.fMassWidth.SetBinContent(xbin, ybin, bin.fMassFitter.GetSignalWidth())
                s.fMassWidth.SetBinError(xbin, ybin, bin.fMassFitter.GetSignalWidthError())
                
    def GenerateSpectrum3D(self, s):
        print("GenerateSpectrum3D not implemented!")
                
    def DrawFitResults(self,bin):
        if bin.fMassFitter is None:
            return
        
        bin.fMassFitter.Draw("same");
        
        fitStatus = int(bin.fMassFitter.GetFitStatus())
        if fitStatus == 0:
            chi2Text = bin.fMassFitter.GetChisquareString().Data()
        else:
            chi2Text = "Fit failed"

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
        for binSet in self.fBinMultiSet.fBinSets.itervalues():
            self.FitInvMassPlotsBinSet(binSet.fName,binSet.fBins,binSet.fFitOptions)

    def FitInvMassPlotsBinSet(self,name,bins,fitOptions):
        pdgMass = ROOT.TDatabasePDG.Instance().GetParticle(421).Mass()

        for i,bin in enumerate(bins):
            if not bin.fInvMassHisto:
                continue
            if self.fTrigger:
                fitterName = "InvMass_{0}_{1}_{2}_fitter".format(self.fTrigger, self.fDMeson, bin.GetName())
            else:
                fitterName = "InvMass_{0}_{1}_fitter".format(self.fDMeson, bin.GetName())

            fitter = self.CreateMassFitter(fitterName)
            bin.SetMassFitter(fitter)
            integral = bin.fInvMassHisto.Integral(1, bin.fInvMassHisto.GetXaxis().GetNbins())
            fitter.GetFitFunction().FixParameter(0, integral) # total integral is fixed
            fitter.GetFitFunction().SetParameter(2, integral / 10) # signal integral (start with small signal)
            fitter.GetFitFunction().SetParLimits(2, 0, integral) # signal integral has to be contained in the total integral
            fitter.GetFitFunction().SetParameter(3, pdgMass) # start fitting using PDG mass
            print("Fitting bin {0}".format(bin.GetTitle()))

            fitter.Fit(bin.fInvMassHisto, fitOptions)

    def PlotInvMassPlots(self):
        for binSet in self.fBinMultiSet.fBinSets.itervalues():
            self.PlotInvMassPlotsBinSet(binSet.fName,binSet.fBins)

    def PlotInvMassLikeSign(self, name, bin):
        spectrum = None
        for s in self.fSpectra.itervalues():
            for binSetName in s.fBins:
                if binSetName == name and s.fAnalysisType == AnalysisType.LikeSign:
                    spectrum = s
                    break
        if not spectrum:
            return None

        eng = None
        for engTest in self.fEngines:
            if engTest.fTrigger == self.fTrigger and engTest.fDMeson == s.fLikeSignTree:
                eng = engTest
                break
        if not eng:
            return None

        binLS = None
        for binTest in eng.fBinMultiSet.fBinSets[name].fBins:
            if binTest.IsSameOf(bin):
                binLS = binTest
                break

        hls = binLS.fInvMassHisto.DrawCopy("hist same")
        hls.SetLineColor(ROOT.kMagenta)
        hls.SetLineStyle(2)
        hls.SetLineWidth(2)
        hls.SetFillStyle(0)
        globalList.append(hls)

        return hls

    def PlotInvMassSideBands(self, name, bin):
        spectrum = None
        for s in self.fSpectra.itervalues():
            for binSetName in s.fBins:
                if binSetName == name and s.fAnalysisType == AnalysisType.SideBand:
                    spectrum = s
                    break
        if not spectrum:
            return None

        invMassSBHistoName = bin.fInvMassHisto.GetName().replace("InvMass", "InvMassSBWindow")
        invMassSigHistoName = bin.fInvMassHisto.GetName().replace("InvMass", "InvMassSigWindow")
        if not invMassSBHistoName in spectrum.fSideBandWindowInvMassHistos.keys() or not invMassSigHistoName in spectrum.fSignalWindowInvMassHistos.keys():
            return None
        hsb = spectrum.fSideBandWindowInvMassHistos[invMassSBHistoName].DrawCopy("hist")
        hsb.SetFillColorAlpha(ROOT.kGreen+2, 0.4)
        hsb.SetFillStyle(1001)
        hsb.SetLineColorAlpha(ROOT.kGreen+2, 0.4)
        hsig = spectrum.fSignalWindowInvMassHistos[invMassSigHistoName].DrawCopy("hist same")
        hsig.SetFillColorAlpha(ROOT.kRed+2, 0.4)
        hsig.SetFillStyle(1001)
        hsig.SetLineColorAlpha(ROOT.kRed+2, 0.4)
        globalList.append(hsb)
        globalList.append(hsig)

        return hsb

    def PlotInvMassPlotsBinSet(self, name, bins):
        cname = "{0}_{1}".format(self.fDMeson, name)
        c = DMesonJetUtils.GenerateMultiCanvas(cname, len(bins))
        self.fCanvases.append(c)
        globalList.append(c)
        for i,bin in enumerate(bins):
            if not bin.fInvMassHisto:
                continue
            pad = c.cd(i+1)
            pad.SetLeftMargin(0.12)
            pad.SetRightMargin(0.05)
            pad.SetTopMargin(0.08)
            pad.SetBottomMargin(0.13)
            if not "SignalOnly" in self.fDMeson:
                SB = self.PlotInvMassSideBands(name, bin)
            else:
                SB = None
            if SB:
                h = bin.fInvMassHisto.DrawCopy("same")
                globalList.append(h)
                SB.SetMaximum(h.GetMaximum()*1.8)
                h = SB
            else:
                h = bin.fInvMassHisto.DrawCopy()
                globalList.append(h)
                h.SetMaximum(h.GetMaximum()*1.8)
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
            if not "SignalOnly" in self.fDMeson:
                self.PlotInvMassLikeSign(name, bin)

class DMesonJetAnalysis:
    def __init__(self, name):
        self.fName = name
        self.fAnalysisEngine = []
        self.fCanvases = []

    def SetProjector(self, projector):
        self.fProjector = projector
        
    def StartAnalysis(self, figTitle, collision, config):
        binMultiSet = DMesonJetProjectors.BinMultiSet()
        for binLists in config["binLists"]:
            if "active" in binLists and not binLists["active"]:
                continue
            limitSetList = []
            for name, binList in binLists["bins"].iteritems():
                limitSetList.append((name, binList))
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
                effToSpectrum = binLists["efficiency"]["apply_to_final_spectrum"]
                if effToSpectrum:
                    fitOptions = "0 L S"
                else:
                    fitOptions = "0 WL S"
            else:
                effWeight = DMesonJetProjectors.SimpleWeight()
                fitOptions = "0 L S"
                effToSpectrum = False
            binMultiSet.AddBinSet(BinSet(binLists["name"], limitSetList, cuts, bin_count_analysis, effWeight, effToSpectrum, fitOptions))

        for trigger in config["trigger"]:
            for d_meson in config["d_meson"]:
                binset_copy = copy.deepcopy(binMultiSet)
                eng = DMesonJetAnalysisEngine(figTitle, collision, trigger, d_meson, 
                                              binset_copy, config["n_mass_bins"], config["min_mass"], config["max_mass"],
                                              config["jets"], config["spectra"], self.fProjector)
                self.fAnalysisEngine.append(eng)
                eng.DoProjections()

        for eng in self.fAnalysisEngine:
            if not "LikeSign" in eng.fDMeson:
                eng.Start(self.fAnalysisEngine)
                eng.CompareSpectra()

        self.CompareSpectra(config["spectra"])

    def CompareSpectra(self, spectra):
        if len(self.fAnalysisEngine) <= 1:
            return
        for s in spectra:
            if "active" in s and not s["active"]:
                continue
            if "axis" in s and len(s["axis"]) != 1:
                continue
            baseline = self.fAnalysisEngine[0].fSpectra[s["name"]].fHistogram.Clone("{0}_copy".format(self.fAnalysisEngine[0].fSpectra[s["name"]].fHistogram.GetName()))
            baseline.SetTitle(self.fAnalysisEngine[0].fDMeson)
            globalList.append(baseline)
            spectraToCompare = []
            for eng in self.fAnalysisEngine[1:]:
                if not eng.fSpectra[s["name"]].fHistogram:
                    continue
                h = eng.fSpectra[s["name"]].fHistogram.Clone("{0}_copy".format(eng.fSpectra[s["name"]].fHistogram.GetName()))
                h.SetTitle(eng.fDMeson)
                globalList.append(h)
                spectraToCompare.append(h)
            if len(spectraToCompare) > 0:
                results = DMesonJetUtils.CompareSpectra(baseline, spectraToCompare, "{0}_SpectraComparison".format(s["name"]))
                for obj in results:
                    if isinstance(obj, ROOT.TCanvas):
                        self.fCanvases.append(obj)
                    globalList.append(obj)

    def SaveRootFile(self, path):
        file = ROOT.TFile.Open("{0}/{1}.root".format(path, self.fName), "recreate")
        file.cd()
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
