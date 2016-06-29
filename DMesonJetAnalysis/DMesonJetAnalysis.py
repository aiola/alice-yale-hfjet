#!/usr/bin/env python
#python program to perform a D meson jet analysis

import ROOT
import math
import DMesonJetProjectors
from DMesonJetBase import *
from array import array
from copy import deepcopy

globalList = []

class DMesonJetAnalysisEngine:
    def __init__(self, trigger, dmeson, binSet, nMassBins, minMass, maxMass, jets, spectra, projector):
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fBinSet = binSet
        self.fNMassBins = nMassBins
        self.fMinMass = minMass
        self.fMaxMass = maxMass
        self.fSpectra = dict()
        for s in spectra:
            name = "{0}_{1}".format(self.fDMeson, s["name"])
            self.fSpectra[s["name"]] = Spectrum(s, name)
        
    def CreateMassFitter(self, name):
        if "D0" in self.fDMeson:
            minMass = self.fMinMass-(self.fMaxMass-self.fMinMass)/2;
            maxMass = self.fMaxMass+(self.fMaxMass-self.fMinMass)/2;
            minFitRange = self.fMinMass;
            maxFitRange = self.fMaxMass;
            startingSigma = 0.01;
            startingSigmaBkg = 0.01;
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
        
    def Start(self):
        self.fProjector.GetInvMassHisograms(self.fTrigger, self.fDMeson, self.fJetDefinitions, 
                                            self.fBinSet, self.fNMassBins, self.fMinMass, self.fMaxMass)
        
        if not "MCTruth" in self.fDMeson:
            self.FitInvMassPlots()
            self.PlotInvMassPlots()
        
        if not "BackgroundOnly" in self.fDMeson:
            self.GenerateSpectra()
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
        c = ROOT.TCanvas("{0}_canvas".format(s.fName), s.fName)
        c.SetLogy()
        c.cd()
        h = s.fHistogram.DrawCopy()
        
        h.SetMarkerColor(ROOT.kBlue+2)
        h.SetMarkerStyle(ROOT.kFullCircle)
        h.SetMarkerSize(0.9)
        h.SetLineColor(ROOT.kBlue+2)
        
        globalList.append(c)
        globalList.append(h)
        
    def PlotSpectrum2D(self, s):
        c = ROOT.TCanvas("{0}_canvas".format(s.fName), s.fName)
        c.SetLogz()
        c.cd()
        h = s.fHistogram.DrawCopy("colz")
        
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
            
    def GenerateSpectrum1D(self, s):
        s.fHistogram = ROOT.TH1D(s.fName, s.fName, len(s.fAxis[0].fBins)-1, array('d',s.fAxis[0].fBins))
        s.fHistogram.GetXaxis().SetTitle(s.fAxis[0].GetTitle())
        s.fHistogram.GetYaxis().SetTitle("counts")
        s.fHistogram.Sumw2()
        for binSetName in s.fBins:
            for bin in self.fBinSet.fBins[binSetName]:
                if bin.fMassFitter is None:
                    print("The bin printed below does not have a mass fitter!")
                    bin.Print()
                    continue
                xbin = s.fHistogram.GetXaxis().FindBin(bin.GetBinCenter(s.fAxis[0].fName))
                if "SignalOnly" in self.fDMeson:
                    s.fHistogram.SetBinContent(xbin, bin.fInvMassHisto.GetEntries())
                    s.fHistogram.SetBinError(xbin, math.sqrt(bin.fInvMassHisto.GetEntries())) 
                else:
                    s.fHistogram.SetBinContent(xbin, bin.fMassFitter.GetSignal())
                    s.fHistogram.SetBinError(xbin, bin.fMassFitter.GetSignalError())
                
    def GenerateSpectrum2D(self, s):
        s.fHistogram = ROOT.TH2D(s.fName, s.fName, len(s.fAxis[0].fBins)-1, array('d',s.fAxis[0].fBins), len(s.fAxis[1].fBins)-1, array('d',s.fAxis[1].fBins))
        s.fHistogram.GetXaxis().SetTitle(s.fAxis[0].GetTitle())
        s.fHistogram.GetYaxis().SetTitle(s.fAxis[1].GetTitle())
        s.fHistogram.GetZaxis().SetTitle("counts")
        s.fHistogram.Sumw2()
        for binSetName in s.fBins:
            for bin in self.fBinSet.fBins[binSetName]:
                if bin.fMassFitter is None:
                    print("The bin printed below does not have a mass fitter!")
                    bin.Print()
                    continue
                xbin = s.fHistogram.GetXaxis().FindBin(bin.GetBinCenter(s.fAxis[0].fName))
                ybin = s.fHistogram.GetYaxis().FindBin(bin.GetBinCenter(s.fAxis[1].fName))
                if "SignalOnly" in self.fDMeson:
                    s.fHistogram.SetBinContent(xbin, ybin, bin.fInvMassHisto.GetEntries())
                    s.fHistogram.SetBinError(xbin, ybin, math.sqrt(bin.fInvMassHisto.GetEntries()))
                else:
                    s.fHistogram.SetBinContent(xbin, ybin, bin.fMassFitter.GetSignal())
                    s.fHistogram.SetBinError(xbin, ybin, bin.fMassFitter.GetSignalError())
                
    def GenerateSpectrum3D(self, s):
        print("GenerateSpectrum3D not implemented!")
                
    def DrawFitResutls(self,bin):
        if bin.fMassFitter is None:
            return
        
        bin.fMassFitter.Draw("same");

        paveSig = ROOT.TPaveText(0.165, 0.710, 0.490, 0.92, "NB NDC")
        globalList.append(paveSig)
        paveSig.SetBorderSize(0)
        paveSig.SetFillStyle(0)
        paveSig.SetTextFont(43)
        paveSig.SetTextSize(14)
        paveSig.SetTextAlign(13)
        paveSig.AddText("{0}, {1}".format(bin.fMassFitter.GetSignalString().Data(), 
                                          bin.fMassFitter.GetBackgroundString().Data()))
        paveSig.AddText("{0}, {1}".format(bin.fMassFitter.GetSignalOverBackgroundString().Data(), 
                                          bin.fMassFitter.GetSignalOverSqrtSignalBackgroundString().Data()))
        fitStatus = int(bin.fMassFitter.GetFitStatus())
        if fitStatus == 0:
            paveSig.AddText(bin.fMassFitter.GetChisquareString().Data())
        else:
            paveSig.AddText("Fit failed")
            
        paveSig.Draw()

        paveFit = ROOT.TPaveText(0.50, 0.66, 0.99, 0.92, "NB NDC")
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
        for name,bins in self.fBinSet.fBins.iteritems():
            self.FitInvMassPlotsBinSet(name,bins)
            
    def FitInvMassPlotsBinSet(self,name,bins):
        pdgMass = ROOT.TDatabasePDG.Instance().GetParticle(421).Mass()
        
        for i,bin in enumerate(bins):
            if not bin.fInvMassHisto:
                continue
            fitter = self.CreateMassFitter("{0}_{1}_{2}_fitter".format(self.fDMeson, name, bin.GetName()))
            bin.SetMassFitter(fitter)
            integral = bin.fInvMassHisto.Integral(1, bin.fInvMassHisto.GetXaxis().GetNbins())
            fitter.GetFitFunction().FixParameter(0, integral) # total integral is fixed
            fitter.GetFitFunction().SetParameter(2, integral / 100) # signal integral (start with very small signal)
            fitter.GetFitFunction().SetParLimits(2, 0, integral) # signal integral has to be contained in the total integral
            fitter.GetFitFunction().SetParameter(3, pdgMass) # start fitting using PDG mass

            fitter.Fit(bin.fInvMassHisto, "0 L E S");
            
    def PlotInvMassPlots(self):
        for name,bins in self.fBinSet.fBins.iteritems():
            self.PlotInvMassPlotsBinSet(name,bins)
    
    def GenerateInvMassCanvas(self, name, n):
        rows = int(math.floor(math.sqrt(n)))
        cols = int(math.ceil(float(n) / rows))
        cname = "{0}_{1}".format(self.fDMeson, name)
        c = ROOT.TCanvas(cname, cname, cols*400, rows*400)
        c.Divide(cols, rows)
        globalList.append(c)
        return c
    
    def PlotInvMassPlotsBinSet(self, name, bins):
        c = self.GenerateInvMassCanvas(name, len(bins))
        for i,bin in enumerate(bins):
            if not bin.fInvMassHisto:
                continue
            pad = c.cd(i+1)
            pad.SetLeftMargin(0.12)
            pad.SetRightMargin(0.05)
            pad.SetTopMargin(0.08)
            pad.SetBottomMargin(0.13)
            h = bin.fInvMassHisto.DrawCopy()
            globalList.append(h)
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
            h.SetMaximum(h.GetMaximum()*1.8)
            htitle = ROOT.TPaveText(0.12, 0.99, 0.95, 0.93, "NB NDC")
            htitle.SetBorderSize(0)
            htitle.SetFillStyle(0)
            htitle.SetTextFont(43)
            htitle.SetTextSize(18)
            htitle.AddText(bin.GetTitle())
            htitle.Draw()
            globalList.append(htitle)
            self.DrawFitResutls(bin)

class DMesonJetAnalysis:
    def __init__(self, name):
        self.fName = name
        self.fAnalysisEngine = []

    def SetProjector(self, projector):
        self.fProjector = projector
        
    def StartAnalysis(self, config):
        binSet = DMesonJetProjectors.BinSet()
        for binLists in config["bins"]:
            binSet.AddBins(name = binLists["name"], jetPtLimits = binLists["jet_pt"], DPtLimits = binLists["d_pt"], ZLimits = binLists["d_z"], 
                           _showJetPt = binLists["show_jet_pt"], _showDPt = binLists["show_d_pt"], _showDZ = binLists["show_d_z"])
        
        for trigger in config["trigger"]:
            for d_meson in config["d_meson"]:
                binset_copy = deepcopy(binSet)
                eng = DMesonJetAnalysisEngine(trigger, d_meson, 
                                              binset_copy, config["n_mass_bins"], config["min_mass"], config["max_mass"],
                                              config["jets"], config["spectra"], self.fProjector)
                self.fAnalysisEngine.append(eng)
                eng.Start()

        spectraToCompare = []

        for s in config["spectra"]:
            dim = 0
            if s["jet_pt"]: dim += 1
            if s["d_pt"]: dim += 1
            if s["d_z"]: dim += 1
            if dim == 1:
                spectraToCompare.append(s["name"])

        self.CompareSpectra(spectraToCompare)

    def CompareSpectra(self, spectraNames):
        if not len(self.fAnalysisEngine) > 1:
            return []
        colors = [ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2]
        for spectrumName in spectraNames:
            cname = "{0}_{1}_SpectraComparison".format(self.fName, spectrumName)
            c = ROOT.TCanvas(cname, cname)
            c.cd()
            globalList.append(c)
            h0 = self.fAnalysisEngine[0].fSpectra[spectrumName].fHistogram.DrawCopy()
            h0.SetMarkerColor(colors[0])
            h0.SetLineColor(colors[0])
            h0.SetMarkerStyle(ROOT.kFullCircle)
            h0.SetMarkerSize(1.2)
            globalList.append(h0)

            cname = "{0}_{1}_SpectraComparison_Ratio".format(self.fName, spectrumName)
            cRatio = ROOT.TCanvas(cname, cname)
            cRatio.cd()
            globalList.append(cRatio)
            leg = ROOT.TLegend(0.15, 0.85, 0.45, 0.7)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.AddEntry(h0, self.fAnalysisEngine[0].fDMeson, "pe")
            globalList.append(leg)
            print(zip(colors[1:len(self.fAnalysisEngine)],self.fAnalysisEngine[1:]))
            for i, (color,eng) in enumerate(zip(colors[1:len(self.fAnalysisEngine)],self.fAnalysisEngine[1:])):
                print("{0} - Working on {1}".format(i, eng.fDMeson))
                c.cd()
                h = eng.fSpectra[spectrumName].fHistogram.DrawCopy("same")
                h.SetMarkerColor(color)
                h.SetLineColor(color)
                h.SetMarkerStyle(ROOT.kOpenCircle)
                h.SetMarkerSize(1.2)
                globalList.append(h)

                leg.AddEntry(h, eng.fDMeson, "pe")

                cRatio.cd()
                hRatio = h.Clone("{0}_Ratio".format(h.GetName()))
                hRatio.SetTitle("{0} Ratio".format(h.GetTitle()))
                hRatio.Divide(h0)
                globalList.append(hRatio)
                if i == 0:
                    hRatio.Draw()
                else:
                    hRatio.Draw("same")
            c.cd()
            leg.Draw()