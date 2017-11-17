#!/usr/bin/env python
# python utilities for D Meson jet analysis

import os
import ROOT
import math
import array

def binom(n, k):
    return math.factorial(n) / math.factorial(n - k) / math.factorial(k)

def GenerateMultiCanvas(cname, n):
    rows = int(math.floor(math.sqrt(n)))
    cols = int(math.ceil(float(n) / rows))
    c = ROOT.TCanvas(cname, cname, cols * 400, rows * 400)
    c.Divide(cols, rows)
    return c

class StatisticSet:
    def __init__(self, n=6, hname="nohist", htitle="nohist"):
        self.fRawMoments = [0.] * (n + 1)
        self.fSumOfWeights = 0.
        self.fSumOfWeights2 = 0.
        self.fEffectiveEntries = 0.
        self.fMedian = None
        self.fMedianError = None
        if hname == "nohist":
            self.fHistogram = None
        else:
            self.fHistogram = ROOT.TH1D(hname, htitle, 50, -1, 1)
            self.fHistogram.Sumw2()

    def Fill(self, y, w=1):
        for i in range(0, len(self.fRawMoments)):
            self.fRawMoments[i] += w * y ** i
        self.fSumOfWeights += w
        self.fSumOfWeights2 += w ** 2
        self.fEffectiveEntries = self.fSumOfWeights ** 2 / self.fSumOfWeights2
        if self.fHistogram:
            self.fHistogram.Fill(y, w)
        # print("effective entries are {0}".format(self.fEffectiveEntries))

    def GetRawMoment(self, n):
        if n < 0 or n > len(self.fRawMoments):
            print("Moment {0} not available!".format(n))

        if self.fSumOfWeights == 0:
            print("No entries to calculate the moment!")
            return 0

        return self.fRawMoments[n] / self.fSumOfWeights

    def GetCentralMoment(self, n):
        if n < 0 or n > len(self.fRawMoments):
            print("Moment {0} not available!".format(n))
        m = 0
        # see http://mathworld.wolfram.com/CentralMoment.html
        for k in range(0, n + 1):
            m += binom(n, k) * (-1) ** (n - k) * self.GetRawMoment(k) * self.GetRawMoment(1) ** (n - k)
        return m

    def GetMean(self):
        return self.GetRawMoment(1)

    def GetVariance(self):
        return self.GetCentralMoment(2)

    def GetStdDev(self):
        if self.GetVariance() > 0: return math.sqrt(self.GetVariance())
        else: return 0;

    def GetStdErrOfMean(self):
        if self.fEffectiveEntries == 0:
            print("No entries to calculate the standard error of the mean!")
            return 0
        # see https://en.wikipedia.org/wiki/Standard_error#Standard_error_of_the_mean
        return self.GetStdDev() / math.sqrt(self.fEffectiveEntries)

    def GetVarianceOfVariance(self):
        if self.fEffectiveEntries == 0:
            print("No entries to calculate the standard error of the mean!")
            return 0
        # see http://mathworld.wolfram.com/SampleVarianceDistribution.html
        return (self.fEffectiveEntries - 1) ** 2 / self.fEffectiveEntries ** 3 * self.GetCentralMoment(4) - (self.fEffectiveEntries - 1) * (self.fEffectiveEntries - 3) / self.fEffectiveEntries ** 3 * self.GetCentralMoment(2) ** 2

    def GetStdErrOfVariance(self):
        # see https://en.wikipedia.org/wiki/Standard_error
        # print("variance of variance is {0}".format(self.GetVarianceOfVariance()))
        return math.sqrt(self.GetVarianceOfVariance())

    def GetStdErrOfStdDev(self):
        if self.GetStdDev() == 0:
            return 0
        # error propagation: StdDev = sqrt(Var)
        return self.GetStdErrOfVariance() / self.GetStdDev() / 2

    def GetMedian(self):
        if not self.fMedian:
            self.RecalculateMedian()

        return self.fMedian

    def GetMedianError(self):
        if not self.fMedianError:
            self.RecalculateMedian()

        return self.fMedianError

    def RecalculateMedian(self):
        self.fMedian = 0
        self.fMedianError = 0
        if not self.fHistogram:
            return
        if self.fHistogram.Integral("width") <= 0:
            return
        self.fHistogram.Scale(1. / self.fHistogram.Integral("width"))
        q = array.array('d', [0.5])
        med = array.array('d', [0.])
        self.fHistogram.GetQuantiles(1, med, q)
        self.fMedian = med[0]
        # see https://web.williams.edu/Mathematics/sjmiller/public_html/BrownClasses/162/Handouts/MedianThm04.pdf
        den = (2 * self.fHistogram.Interpolate(self.fMedian) * math.sqrt(self.fEffectiveEntries))
        if den != 0:
            self.fMedianError = 1. / den
        else:
            self.fMedianError = 0

    def PrintSummary(self):
        print("Mean is {0} +/- {1}".format(self.GetMean(), self.GetStdErrOfMean()))
        print("Standard deviation is {0} +/- {1}".format(self.GetStdDev(), self.GetStdErrOfStdDev()))
        print("Variance of variance is {0}".format(self.GetVarianceOfVariance()))
        print("Central moments: mu0 = {0}, mu1 = {1}, mu2 = {2}, mu3 = {3}, mu4 = {4}, mu5 = {5}, mu6 = {6}".format(self.GetCentralMoment(0), self.GetCentralMoment(1), self.GetCentralMoment(2), self.GetCentralMoment(3), self.GetCentralMoment(4), self.GetCentralMoment(5), self.GetCentralMoment(6)))
        print("Raw moments: mu0 = {0}, mu1 = {1}, mu2 = {2}, mu3 = {3}, mu4 = {4}, mu5 = {5}, mu6 = {6}".format(self.GetRawMoment(0), self.GetRawMoment(1), self.GetRawMoment(2), self.GetRawMoment(3), self.GetRawMoment(4), self.GetRawMoment(5), self.GetRawMoment(6)))
        print("The number of effective entries is {0}".format(self.fEffectiveEntries))
        self.RecalculateMedian()
        if self.fMedian: print("The median is {0} +/- {1}".format(self.GetMedian(), self.GetMedianError()))

class StatisticMultiSet:
    def __init__(self, name, axis, variableName, n=6):
        self.fName = name
        self.fStatisticSets = []
        for min, max in zip(axis.fBins[:-1], axis.fBins[1:]):
            hname = "{0}_{1}_{2}".format(self.fName, int(min), int(max))
            if axis.fName == "jet_pt" or axis.fName == "d_pt":
                htitle = "{0} < {1} < {2} {3};{4};counts".format(int(min), axis.GetVariableName(), int(max), axis.GetVariableUnits(), variableName)
            elif axis.fName == "d_z":
                htitle = "{0:.1f} < {1} < {2:.1f} {3};{4};counts".format(min, axis.GetVariableName(), max, axis.GetVariableUnits(), variableName)
            self.fStatisticSets.append(StatisticSet(n, hname, htitle))
        self.fAxis = axis
        self.fVariableName = variableName
        hname = "{0}_{1}".format(self.fName, "Profile")
        self.fProfileCrossCheck = ROOT.TProfile(hname, hname, len(axis.fBins) - 1, array.array('d', axis.fBins), "i")
        self.fProfileCrossCheck.Sumw2()

    def Fill(self, x, y, w=1):
        self.fProfileCrossCheck.Fill(x, y, w)
        bin = self.fAxis.FindBin(x)
        if bin < 0:
            return
        self.fStatisticSets[bin].Fill(y, w)

    def GenerateMeanHistogram(self, name="Mean"):
        hname = "{0}_{1}".format(self.fName, name)
        hist = ROOT.TH1D(hname, hname, len(self.fAxis.fBins) - 1, array.array('d', self.fAxis.fBins))
        hist.GetXaxis().SetTitle(self.fAxis.GetTitle())
        hist.GetYaxis().SetTitle("<{0}>".format(self.fVariableName))
        for bin in range(1, len(self.fAxis.fBins)):
            hist.SetBinContent(bin, self.fStatisticSets[bin - 1].GetMean())
            hist.SetBinError(bin, self.fStatisticSets[bin - 1].GetStdErrOfMean())
        return hist

    def GenerateMedianHistogram(self, name="Median"):
        hname = "{0}_{1}".format(self.fName, name)
        hist = ROOT.TH1D(hname, hname, len(self.fAxis.fBins) - 1, array.array('d', self.fAxis.fBins))
        hist.GetXaxis().SetTitle(self.fAxis.GetTitle())
        hist.GetYaxis().SetTitle("median({0})".format(self.fVariableName))
        for bin in range(1, len(self.fAxis.fBins)):
            hist.SetBinContent(bin, self.fStatisticSets[bin - 1].GetMedian())
            hist.SetBinError(bin, self.fStatisticSets[bin - 1].GetMedianError())
        return hist

    def GenerateStdDevHistogram(self, name="StdDev"):
        hname = "{0}_{1}".format(self.fName, name)
        hist = ROOT.TH1D(hname, hname, len(self.fAxis.fBins) - 1, array.array('d', self.fAxis.fBins))
        hist.GetXaxis().SetTitle(self.fAxis.GetTitle())
        hist.GetYaxis().SetTitle("#sigma({0})".format(self.fVariableName))
        for bin in range(1, len(self.fAxis.fBins)):
            hist.SetBinContent(bin, self.fStatisticSets[bin - 1].GetStdDev())
            hist.SetBinError(bin, self.fStatisticSets[bin - 1].GetStdErrOfStdDev())
        return hist

    def PrintSummary(self, name):
        print("{0}: printing summary".format(name))
        for bin in range(1, len(self.fAxis.fBins)):
            print("Bin: {0}-{1}".format(self.fAxis.fBins[bin - 1], self.fAxis.fBins[bin]))
            self.fStatisticSets[bin - 1].PrintSummary()
            print("cont={0}, sum={1}, binError={2}, neff={3}".format(self.fProfileCrossCheck.GetBinContent(bin), self.fProfileCrossCheck.GetBinEntries(bin), self.fProfileCrossCheck.GetBinError(bin), self.fProfileCrossCheck.GetBinEffectiveEntries(bin)))
            if self.fStatisticSets[bin - 1].fHistogram:
                print("from histogram: mean = {0} +/- {1}, std dev = {2} +/- {3}".format(self.fStatisticSets[bin - 1].fHistogram.GetMean(), self.fStatisticSets[bin - 1].fHistogram.GetMeanError(), self.fStatisticSets[bin - 1].fHistogram.GetRMS(), self.fStatisticSets[bin - 1].fHistogram.GetRMSError()))
