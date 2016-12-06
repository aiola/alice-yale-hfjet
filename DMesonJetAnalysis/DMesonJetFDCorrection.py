#!/usr/bin/env python
#python program to do extract B feed down correction factors

import yaml
import ROOT
import DMesonJetUtils

class DMesonJetFDCorrection:
    def __init__(self, config, inputPath=None, dmeson=None, jetType=None, jetRadius=None):
        self.fInputPath = inputPath
        self.fDMeson = dmeson
        self.fJetType = jetType
        self.fJetRadius = jetRadius
        if config:
            self.fTrainName = config["train"]
            self.fGenerator = config["gen"]
            self.fProcess = config["proc"]
            self.fTimeStamp = config["ts"]
            self.fSubDir = config["subdir"]
            self.fSpectrumName = config["spectrum"]
            self.LoadFD()
        else:
            self.fFDHistogram = None

    def LoadFD(self):
        filename = "{0}/{1}_{2}_{3}_{4}/{5}/BFeedDown_{2}_{3}_{4}.root".format(
            self.fInputPath, self.fTrainName, self.fGenerator, self.fProcess, self.fTimeStamp, self.fSubDir)
        file = ROOT.TFile(filename)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(filename))
            return
        else:
            print("File {0} open for FD correction".format(filename))
#         rlist = file.Get(self.fSpectrumName)
#         if not rlist:
#             print("Could not find list {0}".format(self.fSpectrumName))
#             return
#         else:
#             print("List {0} loaded".format(self.fSpectrumName))
#         h = rlist.FindObject(self.fSpectrumName)
        hname = "_".join([a for a in [self.fDMeson, "MCTruth", self.fJetType, self.fJetRadius, self.fSpectrumName] if a ])
        h = file.Get(hname)
        if not h:
            print("Could not find histogram {0}".format(hname))
            return
        else:
            print("Histogram {0} loaded".format(hname))
        self.fFDHistogram = h.Clone("FD")
