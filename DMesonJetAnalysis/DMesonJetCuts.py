#!/usr/bin/env python
# python Class to apply cuts

import sys

class DMesonJetCuts:
    def __init__(self, cutList):
        self.fCuts = cutList
        self.InitializeCuts()

    def Dummy(self, obj, min, max):
        return True

    def CutJet(self, dmeson, jet, varname, min, max):
        v = getattr(jet, varname)
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutJetPt(self, dmeson, jet, dummy, min, max):
        v = jet.fPt
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutJetCorrPt(self, dmeson, jet, dummy, min, max):
        v = jet.fCorrPt
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutJetEta(self, dmeson, jet, dummy, min, max):
        v = jet.fEta
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutJetPhi(self, dmeson, jet, dummy, min, max):
        v = jet.fPhi
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutJetZ(self, dmeson, jet, dummy, min, max):
        v = jet.fZ
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutD(self, dmeson, jet, varname, min, max):
        v = getattr(dmeson, varname)
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutDPt(self, dmeson, jet, dummy, min, max):
        v = dmeson.fPt
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutDEta(self, dmeson, jet, dummy, min, max):
        v = dmeson.fEta
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutDPhi(self, dmeson, jet, dummy, min, max):
        v = dmeson.fPhi
        if v < min:
            return False
        if v >= max:
            return False
        return True

    def CutJetInclude(self, dmeson, jet, varname, val, dummy):
        v = getattr(jet, varname)
        if v == val:
            return True
        else:
            return False

    def CutJetExclude(self, dmeson, jet, varname, val, dummy):
        v = getattr(jet, varname)
        if v == val:
            return False
        else:
            return True

    def CutDInclude(self, dmeson, jet, varname, val, dummy):
        v = getattr(dmeson, varname)
        if v == val:
            return True
        else:
            return False

    def CutDExclude(self, dmeson, jet, varname, val, dummy):
        v = getattr(dmeson, varname)
        if v == val:
            return False
        else:
            return True

    def InitializeCuts(self):
        self.fFastCuts = []
        for cut in self.fCuts:
            if "min" in cut: min = cut["min"]
            else: min = -sys.float_info.max
            if "max" in cut: max = cut["max"]
            else: max = +sys.float_info.max
            if "include" in cut:
                min = cut["include"]
                max = None
            elif "exclude" in cut:
                min = cut["exclude"]
                max = None
            if cut["object"] == "d":
                if cut["variable"] == "fPt":
                    self.fFastCuts.append((self.CutDPt, None, min, max))
                elif cut["variable"] == "fEta":
                    self.fFastCuts.append((self.CutDEta, None, min, max))
                elif cut["variable"] == "fPhi":
                    self.fFastCuts.append((self.CutDPhi, None, min, max))
                elif "include" in cut:
                    self.fFastCuts.append((self.CutDInclude, cut["variable"], min, max))
                elif "exclude" in cut:
                    self.fFastCuts.append((self.CutDExclude, cut["variable"], min, max))
                else:
                    self.fFastCuts.append((self.CutD, cut["variable"], min, max))
            elif cut["object"] == "jet":
                if cut["variable"] == "fPt":
                    self.fFastCuts.append((self.CutJetPt, None, min, max))
                elif cut["variable"] == "fCorrPt":
                    self.fFastCuts.append((self.CutJetCorrPt, None, min, max))
                elif cut["variable"] == "fEta":
                    self.fFastCuts.append((self.CutJetEta, None, min, max))
                elif cut["variable"] == "fPhi":
                    self.fFastCuts.append((self.CutJetPhi, None, min, max))
                elif cut["variable"] == "fZ":
                    self.fFastCuts.append((self.CutJetZ, None, min, max))
                elif "include" in cut:
                    self.fFastCuts.append((self.CutJetInclude, cut["variable"], min, max))
                elif "exclude" in cut:
                    self.fFastCuts.append((self.CutJetExclude, cut["variable"], min, max))
                else:
                    self.fFastCuts.append((self.CutJet, cut["variable"], min, max))

    def ApplyCuts(self, dmeson, jet):
        for funct, var, value1, value2 in self.fFastCuts:
            if not funct(dmeson, jet, var, value1, value2): return False
        return True
