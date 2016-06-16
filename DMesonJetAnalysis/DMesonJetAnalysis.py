#!/usr/bin/env python
#python program to perform a D meson jet analysis

import ROOT
import math
import DMesonJetProjectors

globalList = []

class DMesonJetAnalysisEngine:
    def __init__(self, trigger, dmeson, binSet, nMassBins, minMass, maxMass, jets, projector):
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fBinSet = binSet
        self.fNMassBins = nMassBins
        self.fMinMass = minMass
        self.fMaxMass = maxMass
        
    def CreateMassFitter(self, name):
        if self.fDMeson == "D0":
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
        globalList.append(fitter)
        
        return fitter
        
    def Start(self):
        self.fProjector.GetInvMassHisograms(self.fTrigger, self.fDMeson, self.fJetDefinitions, 
                                            self.fBinSet, self.fNMassBins, self.fMinMass, self.fMaxMass)
        
        self.FitInvMassPlots()
        self.PlotInvMassPlots()
        
        
    def DrawFitResutls(self,bin):
        if not bin.fMassFitter:
            return
        
        bin.fMassFitter.Draw("same");

        paveSig = ROOT.TPaveText(0.165, 0.710, 0.490, 0.925, "NB NDC")
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
            fitter = self.CreateMassFitter("{0}_fitter".format(bin.GetName()))
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
        c = ROOT.TCanvas(name, name, cols*400, rows*400)
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
         
        eng = DMesonJetAnalysisEngine(config["trigger"], config["d_meson"], 
                                     binSet, config["n_mass_bins"], config["min_mass"], config["max_mass"],
                                     config["jets"], self.fProjector)
        self.fAnalysisEngine.append(eng)
        eng.Start()
        