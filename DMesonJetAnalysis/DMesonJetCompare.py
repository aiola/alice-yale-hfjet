#!/usr/bin/env python
# python utilities for D Meson jet analysis

import ROOT
import math

class DMesonJetCompare:
    def __init__(self, name):
        self.fName = name
        self.fBaselineHistogram = None
        self.fHistograms = None
        self.fRatios = []
        self.fOptSpectrum = ""
        self.fOptRatio = ""
        self.fYaxisRatio = "Ratio"
        self.fDoSpectraPlot = "logy"
        self.fDoRatioPlot = "lineary"
        self.fCanvasSpectra = None
        self.fCanvasRatio = None
        self.fLegendSpectra = None
        self.fLegendRatio = None
        self.fBaselineRatio = None
        self.fColors = [ROOT.kBlack, ROOT.kBlue + 2, ROOT.kRed + 2, ROOT.kGreen + 2, ROOT.kOrange + 2, ROOT.kAzure + 2, ROOT.kMagenta + 2, ROOT.kCyan + 2, ROOT.kPink + 1, ROOT.kTeal + 2]
        self.fMarkers = [ROOT.kOpenCircle, ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond, ROOT.kFullStar, ROOT.kStar, ROOT.kOpenCircle]
        self.fLines = [1, 2, 9, 5, 7, 10, 4, 3, 6, 8, 9]
        self.fFills = [3001, 3002, 3003, 3004, 3005, 3006, 3007]
        self.fMainHistogram = None
        self.fMainRatioHistogram = None
        self.fMaxSpectrum = -1e30
        self.fMinSpectrum = +1e30
        self.fMaxRatio = -1e30
        self.fMinRatio = +1e30
        self.fRatioRelativeUncertainty = None
        self.fResults = None

    def SetRatioRelativeUncertaintyFromHistogram(self, hist):
        self.fRatioRelativeUncertainty = hist.Clone("{0}_unc".format(hist.GetName()))
        self.fRatioRelativeUncertainty.SetTitle("Uncertainty")
        for ibin in range(0, self.fRatioRelativeUncertainty.GetNbinsX() + 2):
            if self.fRatioRelativeUncertainty.GetBinContent(ibin) != 0:
                self.fRatioRelativeUncertainty.SetBinError(ibin, self.fRatioRelativeUncertainty.GetBinError(ibin) / self.fRatioRelativeUncertainty.GetBinContent(ibin))
                self.fRatioRelativeUncertainty.SetBinContent(ibin, 1)
            else:
                self.fRatioRelativeUncertainty.SetBinError(ibin, 0)
                self.fRatioRelativeUncertainty.SetBinContent(ibin, 0)

    def PrepareSpectraCanvas(self):
        if not self.fCanvasSpectra:
            self.fCanvasSpectra = ROOT.TCanvas(self.fName, self.fName)

        self.fCanvasSpectra.cd()
        if self.fDoSpectraPlot == "logy":
            self.fCanvasSpectra.SetLogy()

        if self.fLegendSpectra:
            self.fLegendSpectra.SetY1(self.fLegendSpectra.GetY1() - 0.04 * (len(self.fHistograms) + 1))
        else:
            self.fLegendSpectra = ROOT.TLegend(0.25, 0.87 - 0.04 * (len(self.fHistograms) + 1), 0.85, 0.87)
            self.fLegendSpectra.SetName("{0}_legend".format(self.fCanvasSpectra.GetName()))
            self.fLegendSpectra.SetFillStyle(0)
            self.fLegendSpectra.SetBorderSize(0)
            self.fLegendSpectra.SetTextFont(43)
            self.fLegendSpectra.SetTextSize(20)

        if "hist" in self.fOptSpectrum:
            self.fBaselineHistogram.SetLineColor(self.fColors[0])
            self.fBaselineHistogram.SetLineWidth(2)
            self.fBaselineHistogram.SetLineStyle(self.fLines[0])
            leg.AddEntry(self.fBaselineHistogram, self.fBaselineHistogram.GetTitle(), "l")
        else:
            self.fBaselineHistogram.SetMarkerColor(self.fColors[0])
            self.fBaselineHistogram.SetLineColor(self.fColors[0])
            self.fBaselineHistogram.SetMarkerStyle(self.fMarkers[0])
            self.fBaselineHistogram.SetMarkerSize(1.2)
            self.fLegendSpectra.AddEntry(self.fBaselineHistogram, self.fBaselineHistogram.GetTitle(), "pe")

        self.fBaselineHistogram.Draw(self.fOptSpectrum)
        if "frac" in self.fBaselineHistogram.GetYaxis().GetTitle():
            self.fCanvasSpectra.SetLeftMargin(0.12)
            self.fBaselineHistogram.GetYaxis().SetTitleOffset(1.4)
        if not "same" in self.fOptSpectrum:
            self.fOptSpectrum += "same"

        if not self.fMainHistogram:
            self.fMainHistogram = self.fBaselineHistogram

        self.fMinSpectrum = min(self.fMinSpectrum, self.fBaselineHistogram.GetBinContent(self.fMainHistogram.GetMinimumBin()))
        self.fMaxSpectrum = max(self.fMaxSpectrum, self.fBaselineHistogram.GetBinContent(self.fMainHistogram.GetMaximumBin()))

    def PrepareRatioCanvas(self):
        cname = "{0}_Ratio".format(self.fName)
        if not self.fCanvasRatio:
            self.fCanvasRatio = ROOT.TCanvas(cname, cname)
        self.fCanvasRatio.cd()
        if self.fDoRatioPlot == "logy":
            self.fDoRatioPlot.SetLogy()

        n = len(self.fHistograms)
        if self.fRatioRelativeUncertainty: n += 1

        if self.fLegendRatio:
            self.fLegendRatio.SetY1(self.fLegendRatio.GetY1() - 0.05 * n)
        else:
            self.fLegendRatio = ROOT.TLegend(0.25, 0.87 - 0.05 * n, 0.85, 0.87)
            self.fLegendRatio.SetName("{0}_legend".format(self.fCanvasRatio.GetName()))
            self.fLegendRatio.SetFillStyle(0)
            self.fLegendRatio.SetBorderSize(0)
            self.fLegendRatio.SetTextFont(43)
            self.fLegendRatio.SetTextSize(20)

        if self.fRatioRelativeUncertainty:
            opt = "e2"
            if "same" in self.fOptRatio:
                opt += "same"
            h = self.fRatioRelativeUncertainty.DrawCopy(opt)
            h.SetFillColor(self.fColors[0])
            h.SetFillStyle(self.fFills[0])
            h.SetLineColor(self.fColors[0])
            h.GetYaxis().SetTitle(self.fYaxisRatio)
            self.fLegendRatio.AddEntry(h, h.GetTitle(), "f")
            self.fResults.append(h)
            if not "same" in self.fOptRatio:
                self.fOptRatio += "same"
            if not self.fMainRatioHistogram:
                self.fMainRatioHistogram = h

    def PlotHistogram(self, color, marker, line, h):
        self.fCanvasSpectra.cd()
        self.fMinSpectrum = min(self.fMinSpectrum, h.GetBinContent(h.GetMinimumBin()))
        self.fMaxSpectrum = max(self.fMaxSpectrum, h.GetBinContent(h.GetMaximumBin()))

        h.Draw(self.fOptSpectrum)
        if "hist" in self.fOptSpectrum:
            h.SetLineColor(color)
            h.SetLineWidth(3)
            h.SetLineStyle(line)
            self.fLegendSpectra.AddEntry(h, h.GetTitle(), "l")
        else:
            h.SetMarkerColor(color)
            h.SetLineColor(color)
            h.SetMarkerStyle(marker)
            h.SetMarkerSize(1.2)
            self.fLegendSpectra.AddEntry(h, h.GetTitle(), "pe")

    def PlotRatio(self, color, marker, line, h):
        self.fCanvasRatio.cd()
        hRatio = h.Clone("{0}_Ratio".format(h.GetName()))
        hRatio.GetYaxis().SetTitle(self.fYaxisRatio)
        if not self.fBaselineRatio:
            self.fBaselineRatio = hRatio
        if "hist" in self.fOptRatio:
            hRatio.SetLineColor(color)
            hRatio.SetLineWidth(3)
            hRatio.SetLineStyle(line)
            self.fLegendRatio.AddEntry(hRatio, h.GetTitle(), "l")
        else:
            hRatio.SetMarkerColor(color)
            hRatio.SetLineColor(color)
            hRatio.SetMarkerStyle(marker)
            hRatio.SetMarkerSize(1.2)
            self.fLegendRatio.AddEntry(hRatio, h.GetTitle(), "pe")
        self.fRatios.append(hRatio)
        hRatio.SetTitle("{0} Ratio".format(h.GetTitle()))
        hRatio.Divide(self.fBaselineHistogram)
        hRatio.Draw(self.fOptRatio)
        if not self.fMainRatioHistogram:
            self.fMainRatioHistogram = hRatio
        self.fMinRatio = min(self.fMinRatio, hRatio.GetBinContent(hRatio.GetMinimumBin()))
        self.fMaxRatio = max(self.fMaxRatio, hRatio.GetBinContent(hRatio.GetMaximumBin()))

        if not "same" in self.fOptRatio:
            self.fOptRatio += "same"

    def CompareSpectra(self, baseline, histos):
        self.fResults = []
        print("CompareSpectra: {0}".format(self.fName))
        self.fBaselineHistogram = baseline
        self.fHistograms = histos
        print("Baseline: {0}".format(self.fBaselineHistogram.GetName()))
        for s in self.fHistograms:
            print(s.GetName())

        if self.fDoSpectraPlot:
            self.PrepareSpectraCanvas()

        if self.fDoRatioPlot:
            self.PrepareRatioCanvas()

        for color, marker, line, h in zip(self.fColors[1:], self.fMarkers[1:], self.fLines[1:], self.fHistograms):
            if self.fDoSpectraPlot:
                self.PlotHistogram(color, marker, line, h)
            if self.fDoRatioPlot:
                self.PlotRatio(color, marker, line, h)
        self.AdjustYLimits()
        self.GenerateResults()
        return self.fResults

    def GenerateResults(self):
        self.fResults.extend(self.fRatios)
        if self.fCanvasSpectra: self.fResults.append(self.fCanvasSpectra)
        if self.fCanvasRatio: self.fResults.append(self.fCanvasRatio)
        if self.fLegendSpectra: self.fResults.append(self.fLegendSpectra)
        if self.fLegendRatio: self.fResults.append(self.fLegendRatio)
        if self.fBaselineRatio: self.fResults.append(self.fBaselineRatio)
        if self.fRatioRelativeUncertainty: self.fResults.append(self.fRatioRelativeUncertainty)

    def AdjustYLimits(self):
        if self.fDoRatioPlot:
            if self.fDoRatioPlot == "logy":
                max = self.fMaxRatio * 5
                min = self.fMinRatio / 2
            else:
                max = self.fMaxRatio + (self.fMaxRatio - self.fMinRatio) * 0.9
                if self.fMinRatio < 0.2:
                    min = 0
                else:
                    min = self.fMinRatio - (self.fMaxRatio - self.fMinRatio) * 0.9
            self.fMainRatioHistogram.SetMinimum(min)
            self.fMainRatioHistogram.SetMaximum(max)
            self.fCanvasRatio.cd()
            self.fLegendRatio.Draw()

        if self.fDoSpectraPlot:
            if self.fDoSpectraPlot == "logy":
                max = self.fMaxSpectrum * 5
                min = self.fMinSpectrum / 2
            else:
                max = self.fMaxSpectrum + (self.fMaxSpectrum - self.fMinSpectrum) * 0.9
                if self.fMinSpectrum < 0.2:
                    min = 0
                else:
                    min = self.fMinSpectrum - (self.fMaxSpectrum - self.fMinSpectrum) * 0.9
            self.fMainHistogram.SetMinimum(min)
            self.fMainHistogram.SetMaximum(max)
            self.fCanvasSpectra.cd()
            self.fLegendSpectra.Draw()
