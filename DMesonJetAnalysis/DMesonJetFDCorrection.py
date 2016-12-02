#!/usr/bin/env python
#python program to do extract B feed down correction factors

import yaml
import ROOT
import DMesonJetUtils

class DMesonJetFDCorrection:
    def __init__(self, config=None):
        if config:
            self.fTrainName = config["train"]
            self.fGenerator = config["gen"]
            self.fProcess = config["proc"]
            self.fTimeStamp = config["ts"]
            self.fSudDir = config["subdir"]
            self.fSpectrumName = config["spectrum"]
            self.GetFDcounts = self.GetFDcountsInternal
        else:
            self.GetFDcounts = self.GetFDcountsDummy

    def GetFDcountsDummy(self, jetpt):
        return 0., 0.

    def GetFDcountsInternal(self, jetpt):
        return 0., 0.

    def GenerateFDhistogram(self, axis, cuts):
        h = DMesonJetUtils.BuildHistogram(axis, "FD", "counts")
        if len(axis) == 1:
            self.Fill1D(h, cuts)
        elif len(axis) == 2:
            self.Fill2D(h, cuts)
        else:
            print("Error - FD correction: cannot handle histograms with dimension > 2!")
        return h

    def Fill1D(self, h, cuts):
        for xbin in range(1, h.GetNbinsX()+1):
            pass