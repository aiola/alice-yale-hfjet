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
        pass