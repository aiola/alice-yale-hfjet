#!/usr/bin/env python
# python program to do extract B feed down correction factors

import yaml
import ROOT
import DMesonJetUtils

class DMesonJetFDCorrection:
    def __init__(self, config, spectrumName, inputPath=None, dmeson=None, jetType=None, jetRadius=None):
        self.fInputPath = inputPath
        self.fSpectrumName = spectrumName
        # For the moment these parameters are not used
        #####
        self.fDMeson = dmeson
        self.fJetType = jetType
        self.fJetRadius = jetRadius
        #####
        if config:
            self.fFileName = config["file_name"]
            self.fCentralPoints = config["central_points"]
            self.fFDSpectrumName = config["spectrum"]
            self.LoadFD()
        else:
            self.fFDHistogram = None

    def LoadFD(self):
        filename = "{0}/{1}".format(self.fInputPath, self.fFileName)
        file = ROOT.TFile(filename)
        if not file or file.IsZombie():
            print("Could not open file {0}".format(filename))
            return
        else:
            print("File {0} open for FD correction".format(filename))
        rlist = file.Get(self.fCentralPoints)
        if not rlist:
            print("Could not find list {0}".format(self.fCentralPoints))
            return
        else:
            print("List {0} loaded".format(self.fCentralPoints))
        rlist2 = rlist.FindObject(self.fSpectrumName)
        if not rlist2:
            print("Could not find list {0}".format(self.fSpectrumName))
            return
        else:
            print("List {0} loaded".format(self.fSpectrumName))
        h = rlist2.FindObject(self.fFDSpectrumName)
        if not h:
            print("Could not find histogram {0}".format(self.fFDSpectrumName))
            return
        else:
            print("Histogram {0} loaded".format(self.fFDSpectrumName))
        self.fFDHistogram = h.Clone("FD")

        # Loading systematic uncertainty
        systListName = "SystematicUncertainty"
        systList = file.Get(systListName)
        if not systList:
            print("Could not find list {0}".format(systListName))
            return
        else:
            print("List {0} loaded".format(systListName))

        detLevListName = "{0}/{1}".format(self.fSpectrumName, self.fFDSpectrumName)
        detLevList = systList.FindObject(detLevListName)
        if not detLevList:
            print("Could not find list {0}".format(detLevListName))
            systList.ls()
            return
        else:
            print("List {0} loaded".format(detLevListName))

        detLevUpSystName = "{0}_UpperSyst".format(self.fFDSpectrumName)
        h = detLevList.FindObject(detLevUpSystName)
        if not h:
            print("Could not find hist {0}".format(detLevUpSystName))
            return
        else:
            print("Hist {0} loaded".format(detLevUpSystName))
        self.fFDUpSystUncHistogram = h.Clone("FD_DetectorLevelUpSyst")

        detLevLowSystName = "{0}_LowerSyst".format(self.fFDSpectrumName)
        h = detLevList.FindObject(detLevLowSystName)
        if not h:
            print("Could not find hist {0}".format(detLevLowSystName))
            return
        else:
            print("Hist {0} loaded".format(detLevLowSystName))
        self.fFDLowSystUncHistogram = h.Clone("FD_DetectorLevelLowSyst")

        detLevSystName = "{0}_CentralAsymmSyst".format(self.fFDSpectrumName)
        h = detLevList.FindObject(detLevSystName)
        if not h:
            print("Could not find hist {0}".format(detLevSystName))
            return
        else:
            print("Hist {0} loaded".format(detLevSystName))
        self.fFDSystUncGraph = h.Clone("FD_DetectorLevelSyst")

        file.Close()
