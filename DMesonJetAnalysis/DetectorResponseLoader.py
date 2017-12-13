#!/usr/bin/env python
# python class to handle reconstruction efficiency of D-jets

import ROOT

import DMesonJetUtils


class DetectorResponseLoader:

    def __init__(self, responseFile, dmeson, jetName, spectrumName, suffix, nbins, bins):
        self.fResponseLoaded = False
        self.fResponseFile = responseFile
        self.fDMeson = dmeson
        self.fSuffix = suffix
        self.fJetName = jetName
        self.fSpectrumName = spectrumName
        self.fBins = bins
        self.fNBins = nbins
        self.fDebug = True

    @classmethod
    def fromConfigAdvanced(cls, conf, trigger, DMesonDef, jetName, nbins, bins):
        eff_file_name = conf["file_name"]
        responseFile = ROOT.TFile.Open(eff_file_name)
        if not responseFile or responseFile.IsZombie():
            print("Could not open file '{0}'! Efficiency could not be loaded!".format(eff_file_name))
            return None
        dmeson = "_".join([obj for obj in [trigger, DMesonDef] if obj])
        spectrumName = conf["list_name"]
        suffix = conf["object_name"].replace("Efficiency", "")
        return cls(responseFile, dmeson, jetName, spectrumName, suffix, nbins, bins)

    @classmethod
    def fromConfig(cls, conf, trigger, DMesonDef, jetName):
        eff_file_name = conf["file_name"]
        responseFile = ROOT.TFile.Open(eff_file_name)
        if not responseFile or responseFile.IsZombie():
            print("Could not open file '{0}'! Efficiency could not be loaded!".format(eff_file_name))
            return None
        dmeson = "_".join([obj for obj in [trigger, DMesonDef] if obj])
        spectrumName = conf["list_name"]
        suffix = conf["object_name"].replace("Efficiency", "")
        return cls(responseFile, dmeson, jetName, spectrumName, suffix, 0, None)

    def GetEfficiencyObject(self):
        if not self.fResponseLoaded: self.LoadResponse()
        if not self.fResponseLoaded: return None
        return self.fEfficiencyObject

    def GetResponseMatrixObject(self):
        if not self.fResponseLoaded: self.LoadResponse()
        if not self.fResponseLoaded: return None
        return self.fResponseMatrixObject

    def LoadResponse(self):
        rlistName = "{}_{}_{}".format(self.fDMeson, self.fJetName, self.fSpectrumName)
        rlist = self.fResponseFile.Get(rlistName)
        if not rlist:
            print("Could not get list {} from file {}".format(rlistName, self.fResponseFile.GetName()))
            exit(1)
        truthName = "{}_{}_{}_Truth{}".format(self.fDMeson, self.fJetName, self.fSpectrumName, self.fSuffix)
        truth = rlist.FindObject(truthName)
        if truth:
            print("Truth loaded from object '{}/{}' in file '{}'".format(rlistName, truthName, self.fResponseFile.GetName()))
        else:
            print("Could not get histogram {}".format(truthName))
            exit(1)
        if self.fNBins > 0:
            truth_coarse = truth.Rebin(self.fNBins, "{}_coarse".format(truthName), self.fBins)
        else:
            truth_coarse = truth.Clone()
        respName = "{}_{}_{}_DetectorResponse{}".format(self.fDMeson, self.fJetName, self.fSpectrumName, self.fSuffix)
        resp = rlist.FindObject(respName)
        if resp:
            print("Response loaded from object '{}/{}' in file '{}'".format(rlistName, respName, self.fResponseFile.GetName()))
        else:
            print("Could not get histogram {}".format(respName))
            exit(1)
        if self.fNBins > 0:
            self.fResponseMatrixObject = DMesonJetUtils.Rebin2D_fromBins(resp, "{}_coarse".format(respName), self.fNBins, self.fBins, self.fNBins, self.fBins)
        else:
            self.fResponseMatrixObject = resp.Clone()
        recoTruthName = "{}_{}_{}_ReconstructedTruth{}".format(self.fDMeson, self.fJetName, self.fSpectrumName, self.fSuffix)
        recoTruth = rlist.FindObject(recoTruthName)
        if recoTruth:
            print("Reconstructed Truth loaded from object '{}/{}' in file '{}'".format(rlistName, recoTruthName, self.fResponseFile.GetName()))
        else:
            print("Could not get histogram {}".format(recoTruthName))
            exit(1)
        if self.fNBins > 0:
            self.fEfficiencyObject = recoTruth.Rebin(self.fNBins, "{}_coarse".format(recoTruthName), self.fBins)
        else:
            self.fEfficiencyObject = recoTruth.Clone(truth_coarse.GetName().replace("Truth", "Efficiency"))
        self.fEfficiencyObject.Divide(truth_coarse)
        self.fResponseLoaded = True

        if self.fDebug:
            effName = "{}_{}_{}_Efficiency{}".format(self.fDMeson, self.fJetName, self.fSpectrumName, self.fSuffix)
            eff = rlist.FindObject(effName)
            if truth:
                print("Original efficiency (not used) loaded from object '{}/{}' in file '{}'".format(rlistName, effName, self.fResponseFile.GetName()))
            else:
                print("Could not get histogram {}".format(effName))
                exit(1)
            for ibin in range(1, eff.GetNbinsX() + 1):
                print("Efficiency {0} for bin {1},{2}".format(eff.GetBinContent(ibin),
                                                              eff.GetXaxis().GetBinLowEdge(ibin),
                                                              eff.GetXaxis().GetBinUpEdge(ibin)))


class DMesonJetEfficiency:

    def __init__(self, obj, doNotInvert=False):
        self.fVariables = [{"obj": "d", "name": "fPt"}]
        if not obj:
            self.GetEfficiencyWeightTH1ForPt = self.GetEfficiencyWeightDummy
            self._GetEfficiencyWeight = self.GetEfficiencyWeightDummy
        else:
            if not isinstance(obj, ROOT.TObject):
                print("EfficiencyWeightCalculator.__init__: Not a root object!!")
                print(obj)
                exit(1)
            self.fRootObject = obj.Clone()
            self.fDoNotInvert = doNotInvert
            if isinstance(self.fRootObject, ROOT.TGraph):
                print("The ROOT object is of type TGraph")
                self._GetEfficiencyWeight = self.GetEfficiencyWeightTGraph
            elif isinstance(self.fRootObject, ROOT.TF1):
                print("The ROOT object is of type TF1")
                self._GetEfficiencyWeight = self.GetEfficiencyWeightTF1
            elif isinstance(self.fRootObject, ROOT.TH2):
                print("The ROOT object is of type TH2")
                self._GetEfficiencyWeight = self.GetEfficiencyWeightTH2
            elif isinstance(self.fRootObject, ROOT.TH1):
                print("The ROOT object is of type TH1")
                self._GetEfficiencyWeight = self.GetEfficiencyWeightTH1
                self.PrintEfficiencyValues()
            else:
                print("The ROOT object type is not recognized!!")
                self.fRootObject.Print()
                exit(1)
            self.GetEfficiencyWeightTH1ForPt = self._GetEfficiencyWeightTH1ForPt

    def PrintEfficiencyValues(self):
        hist = self.fRootObject
        print("Efficiency from object '{}'".format(hist.GetName()))
        for ibin in range(1, hist.GetNbinsX() + 1):
            print("Efficiency {0} for bin {1},{2}".format(hist.GetBinContent(ibin),
                                                          hist.GetXaxis().GetBinLowEdge(ibin),
                                                          hist.GetXaxis().GetBinUpEdge(ibin)))

    @classmethod
    def fromFileName(cls, filename, listname, objectname, doNotInvert=False):
        file = ROOT.TFile.Open(filename)
        if not file or file.IsZombie():
            print("Could not open file '{0}'! Efficiency could not be loaded!".format(filename))
            obj = None
        if listname:
            rlist = file.Get(listname)
            if not rlist:
                print("Could not find list '{0}' in file '{1}'".format(listname, file.GetName()))
                file.ls()
                obj = None
            obj = rlist.FindObject(objectname)
            if not obj:
                print("Could not find object '{0}' in list '{1}'".format(objectname, listname))
                rlist.Print()
                obj = None
        else:
            obj = file.Get(objectname)
            if not obj:
                print("Could not find object '{0}' in file '{1}'".format(objectname, file.GetName()))
                rlist.Print()
                obj = None
        if file: file.Close()
        if obj: print("Efficiency weights loaded from object '{0}/{1}' in file '{2}'".format(listname, objectname, filename))
        return cls(obj, doNotInvert)

    @classmethod
    def fromFunctionName(cls, functionName, minX, maxX, doNotInvert=False):
        obj = ROOT.TF1(functionName, functionName, minX, maxX)
        print("Function '{}' created. Value at extrema [{}, {}] are [{}, {}]".format(functionName, minX, maxX, obj.Eval(minX), obj.Eval(maxX)))
        return cls(obj, doNotInvert)

    @classmethod
    def fromConfigAdvanced(cls, conf, trigger, DMesonDef, jetName, nbins, bins):
        if conf:
            doNotInvert = False
            if "do_not_invert" in conf: doNotInvert = conf["do_not_invert"]
            if "file_name" in conf:
                loader = DetectorResponseLoader.fromConfigAdvanced(conf, trigger, DMesonDef, jetName, nbins, bins)
                if loader:
                    effWeight = cls(loader.GetEfficiencyObject(), doNotInvert)
                else:
                    effWeight = cls(None)
            elif "function_name" in conf:
                effWeight = cls.fromFunctionName(conf["function_name"], conf["min"], conf["max"], doNotInvert)
            else:
                print("Configuration for efficiency not recognized")
                print(conf)
                exit(1)

            if effWeight and "variables" in conf:
                effWeight.fVariables = conf["variables"]
        else:
            effWeight = cls(None)

        return effWeight

    @classmethod
    def fromConfig(cls, conf, trigger, DMesonDef, jetName):
        if conf:
            doNotInvert = False
            if "do_not_invert" in conf: doNotInvert = conf["do_not_invert"]
            if "file_name" in conf:
                eff_file_name = conf["file_name"]
                eff_list_name = "_".join([obj for obj in [trigger, DMesonDef, jetName, conf["list_name"]] if obj])
                eff_obj_name = "_".join([obj for obj in [trigger, DMesonDef, jetName, conf["list_name"], conf["object_name"]] if obj])
                effWeight = cls.fromFileName(eff_file_name, eff_list_name, eff_obj_name, doNotInvert)
            elif "function_name" in conf:
                effWeight = cls.fromFunctionName(conf["function_name"], conf["min"], conf["max"], doNotInvert)
            else:
                print("Configuration for efficiency not recognized")
                print(conf)
                exit(1)

            if effWeight and "variables" in conf:
                effWeight.fVariables = conf["variables"]
        else:
            effWeight = cls(None)

        return effWeight

    def GetEfficiencyWeightDummy(self, dummy):
        return 1

    def GetEfficiencyWeight(self, dmeson, jet):
        values = []
        for var in self.fVariables:
            if var["obj"] == "d": values.append(getattr(dmeson, var["name"]))
            if var["obj"] == "jet": values.append(getattr(jet, var["name"]))
        return self._GetEfficiencyWeight(values)

    def GetEfficiencyWeightTGraph(self, var):
        eff = self.fRootObject.Eval(var[0])
        if eff <= 0:
            return 0
        elif self.fDoNotInvert:
            return eff
        else:
            return 1. / eff

    def GetEfficiencyWeightTF1(self, var):
        eff = self.fRootObject.Eval(var[0])
        if eff <= 0:
            return 0
        elif self.fDoNotInvert:
            return eff
        else:
            return 1. / eff

    def GetEfficiencyWeightTH2(self, var):
        eff = self.fRootObject.GetBinContent(self.fRootObject.FindBin(var[0], var[1]))
        if eff <= 0:
            return 0
        elif self.fDoNotInvert:
            return eff
        else:
            return 1. / eff

    def GetEfficiencyWeightTH1(self, var):
        eff = self.fRootObject.GetBinContent(self.fRootObject.FindBin(var[0]))
        if eff <= 0:
            return 0
        elif self.fDoNotInvert:
            return eff
        else:
            return 1. / eff

    def _GetEfficiencyWeightTH1ForPt(self, pt):
        eff = self.fRootObject.GetBinContent(self.fRootObject.FindBin(pt))
        if eff <= 0:
            return 0
        elif self.fDoNotInvert:
            return eff
        else:
            return 1. / eff
