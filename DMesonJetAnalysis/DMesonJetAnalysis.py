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
    def __init__(self, trigger, dmeson, binSet, nMassBins, minMass, maxMass, jets, spectra, projector):
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fBinSet = binSet
        self.fNMassBins = nMassBins
        self.fMinMass = minMass
        self.fMaxMass = maxMass
        self.fSpectra = collections.OrderedDict()
        self.fCanvases = []
        for s in spectra:
            if "active" in s and not s["active"]:
                continue
            name = "{0}_{1}".format(self.fDMeson, s["name"])
            self.fSpectra[s["name"]] = Spectrum(s, name)
            
    def CompareSpectra(self):
        spectraToCompare = []
        for name,s in self.fSpectra.iteritems():
            if len(s.fAxis) != 1:
                continue

            h = s.fNormHistogram.Clone("{0}_copy".format(s.fNormHistogram.GetName()))
            if s.fTitle:
                h.SetTitle(s.fTitle)
            globalList.append(h)
            spectraToCompare.append(h)

        results = DMesonJetUtils.CompareSpectra(spectraToCompare[0], spectraToCompare[1:], "{0}_SpectraComparison".format(self.fDMeson), "", "hist")
        for obj in results:
            if obj and isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)
            globalList.append(obj)        

    def SaveRootFile(self, file):
        file.cd()
        for rlist in self.fBinSet.GenerateInvMassRootLists():
            rlist.Write("{0}_{1}".format(self.fDMeson, rlist.GetName()), ROOT.TObject.kSingleKey)
        for s in self.fSpectra.itervalues():
            if s.fHistogram:
                s.fHistogram.Write()
            if s.fNormHistogram:
                s.fNormHistogram.Write()

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))
        
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

        self.fEvents = self.fProjector.fTotalEvents
        
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
        pave.AddText("ALICE Work In Progress")
        pave.AddText("pp #sqrt{#it{s}} = 7 TeV")
        pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R}=0.4")
        pave.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
        pave.AddText("|#eta_{jet}| < 0.5, #it{p}_{T,D} > 0.5 GeV/#it{c}")
        pave.Draw()

        globalList.append(c)
        globalList.append(h)
        globalList.append(pave)

    def PlotSpectrum2D(self, s):
        c = ROOT.TCanvas("{0}_canvas".format(s.fName), s.fName)
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

    def GenerateSpectrum1D(self, s):
        s.fHistogram = ROOT.TH1D(s.fName, s.fName, len(s.fAxis[0].fBins)-1, array.array('d',s.fAxis[0].fBins))
        s.fHistogram.GetXaxis().SetTitle(s.fAxis[0].GetTitle())
        s.fHistogram.GetYaxis().SetTitle("counts")
        s.fHistogram.Sumw2()
        for binSetName in s.fBins:
            for bin in self.fBinSet.fBins[binSetName][0]:
                if bin.fMassFitter is None:
                    print("The bin printed below does not have a mass fitter!")
                    bin.Print()
                    continue
                xbin = s.fHistogram.GetXaxis().FindBin(bin.GetBinCenter(s.fAxis[0].fName))
                if "SignalOnly" in self.fDMeson:
                    s.fHistogram.SetBinContent(xbin, bin.fInvMassHisto.Integral())
                    s.fHistogram.SetBinError(xbin, math.sqrt(bin.fInvMassHisto.Integral()))
                else:
                    s.fHistogram.SetBinContent(xbin, bin.fMassFitter.GetSignal())
                    s.fHistogram.SetBinError(xbin, bin.fMassFitter.GetSignalError())
                
    def GenerateSpectrum2D(self, s):
        s.fHistogram = ROOT.TH2D(s.fName, s.fName, len(s.fAxis[0].fBins)-1, array.array('d',s.fAxis[0].fBins), len(s.fAxis[1].fBins)-1, array.array('d',s.fAxis[1].fBins))
        s.fHistogram.GetXaxis().SetTitle(s.fAxis[0].GetTitle())
        s.fHistogram.GetYaxis().SetTitle(s.fAxis[1].GetTitle())
        s.fHistogram.GetZaxis().SetTitle("counts")
        s.fHistogram.Sumw2()
        for binSetName in s.fBins:
            for bin in self.fBinSet.fBins[binSetName][0]:
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
        for name,(bins,_) in self.fBinSet.fBins.iteritems():
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

            fitter.Fit(bin.fInvMassHisto, "0 WL S");
            
    def PlotInvMassPlots(self):
        for name,(bins,_) in self.fBinSet.fBins.iteritems():
            self.PlotInvMassPlotsBinSet(name,bins)
    
    def GenerateInvMassCanvas(self, name, n):
        rows = int(math.floor(math.sqrt(n)))
        cols = int(math.ceil(float(n) / rows))
        cname = "{0}_{1}".format(self.fDMeson, name)
        c = ROOT.TCanvas(cname, cname, cols*400, rows*400)
        c.Divide(cols, rows)
        self.fCanvases.append(c)
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
        self.fCanvases = []

    def SetProjector(self, projector):
        self.fProjector = projector
        
    def StartAnalysis(self, config):
        binSet = DMesonJetProjectors.BinSet()
        for binLists in config["binLists"]:
            if "active" in binLists and not binLists["active"]:
                continue
            limitSetList = []
            for name, binList in binLists["bins"].iteritems():
                limitSetList.append((name, binList))
            if "cuts" in binLists:
                cuts = binLists["cuts"]
            binSet.AddBins(binLists["name"], limitSetList, cuts)
        
        for trigger in config["trigger"]:
            for d_meson in config["d_meson"]:
                binset_copy = copy.deepcopy(binSet)
                eng = DMesonJetAnalysisEngine(trigger, d_meson, 
                                              binset_copy, config["n_mass_bins"], config["min_mass"], config["max_mass"],
                                              config["jets"], config["spectra"], self.fProjector)
                self.fAnalysisEngine.append(eng)
                eng.Start()
                eng.CompareSpectra()

        self.CompareSpectra(config["spectra"])

    def CompareSpectra(self, spectra):
        if len(self.fAnalysisEngine) <= 1:
            return
        for s in spectra:
            if "active" in s and not s["active"]:
                continue
            if len(s["axis"]) != 1:
                continue
            baseline = self.fAnalysisEngine[0].fSpectra[s["name"]].fHistogram.Clone("{0}_copy".format(self.fAnalysisEngine[0].fSpectra[s["name"]].fHistogram.GetName()))
            baseline.SetTitle(self.fAnalysisEngine[0].fDMeson)
            globalList.append(baseline)
            spectraToCompare = []
            for eng in self.fAnalysisEngine[1:]:
                h = eng.fSpectra[s["name"]].fHistogram.Clone("{0}_copy".format(eng.fSpectra[s["name"]].fHistogram.GetName()))
                h.SetTitle(eng.fDMeson)
                globalList.append(h)
                spectraToCompare.append(h)

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
