#!/usr/bin/env python
#python program to do extract B feed down correction factors

import yaml
import ROOT
import DMesonJetUtils

class DMesonJetFDCorrection:
    def __init__(self, inputPath, config=None):
        self.fInputPath = inputPath
        self.fFDHistograms = dict()
        if config:
            self.fTrainName = config["train"]
            self.fGenerator = config["gen"]
            self.fProcess = config["proc"]
            self.fTimeStamp = config["ts"]
            self.fSudDir = config["subdir"]
            self.fSpectrumName = config["spectrum"]
            self.GetFDcounts = self.GetFDcountsInternal
            self.LoadFD()
        else:
            self.GetFDcounts = self.GetFDcountsDummy

    def LoadFD(self):
        filename = "{0}/{1}_{2}_{3}_{4}/{5}/{1}_{2}_{3}_{4}.root".format(
            self.fInputPath, self.fTrainName, self.fGenerator, self.fProcess, self.fTimeStamp, self.fSubDir)
        file = ROOT.TFile(filename)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(filename))
            return
        else:
            print("File {0} open for FD correction".format(filename))
        rlist = file.Get(self.fSpectrumName)
        if not rlist:
            print("Could not find list {0}".format(self.fSpectrumName))
            return
        else:
            print("List {0} loaded".format(self.fSpectrumName))
        h = rlist.FindObject(self.fSpectrumName)
        if not h:
            print("Could not find histogram {0}".format(self.fSpectrumName))
            return
        else:
            print("Histogram {0} loaded".format(self.fSpectrumName))
        self.fFDHistogram = h.Clone("FD")

    def GetFDcountsDummy(self, jetpt):
        return 0., 0.

    def GetFDcountsInternal(self, jetpt):
        return 0., 0.

    def GetFDhistogram(self, axis):
        key = "_".join([axis.fName]+axis.fBins)
        if key in self.fFDHistograms:
            return self.fFDHistograms[key]
        else:
            return self.GenerateFDhistogram(axis)

    def GenerateFDhistogram(self, axis):
        key = "_".join([axis.fName]+axis.fBins)
        print("Generating FD spectrum {0}".format(key))
        h = DMesonJetUtils.BuildHistogram(axis, key, "counts")
        if len(axis) == 1 and axis[0].fName == "jet_pt":
            self.Fill1D(h)
        else:
            print("Error - FD correction: cannot handle histograms with dimension > 1 or different from jet pt!")
        self.fFDHistograms[key] = h
        return h

    def Fill1D(self, h, cuts):
        for xbin in range(1, h.GetNbinsX()+1):
            v,e = self.GetFDcounts(h.GetXaxis().GetBinCenter(xbin))
            h.SetBinContent(xbin, v)
            h.SetBinError(xbin, e)
