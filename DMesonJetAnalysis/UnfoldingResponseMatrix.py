#!/usr/local/bin/python
#python program to perform a D meson jet unfolding

import ROOT

class ResponseMatrix:
    def __init__(self, name, rooUnfoldResponse, response, truth, smallBinResponse, smallBinTruth, normResp):
        self.fName = name
        self.fRooUnfoldResponse = rooUnfoldResponse

        self.fResponse = response.Clone("{0}_ResponseMatrix".format(self.fName))
        self.fResponse.SetTitle("{0} Response Matrix".format(self.fName))

        self.fNormResponse = normResp.Clone("{0}_NormalizedResponseMatrix".format(self.fName))
        self.fNormResponse.SetTitle("{0} Normalized Response Matrix".format(self.fName))

        self.fTruth = truth.Clone("{0}_Truth".format(self.fName))
        self.fTruth.SetTitle("{0} Truth".format(self.fName))

        self.fSmallBinResponse = smallBinResponse.Clone("{0}_SmallBinResponseMatrix".format(self.fName))
        self.fSmallBinResponse.SetTitle("{0} Small-Bin Response Matrix".format(self.fName))

        self.fSmallBinTruth = smallBinTruth.Clone("{0}_SmallBinTruth".format(self.fName))
        self.fSmallBinTruth.SetTitle("{0} Small-Bin Truth".format(self.fName))
        
        self.fReconstructedGenLev = self.fResponse.ProjectionY("{0}_ReconstructedGenLev".format(self.fName), 1, self.fResponse.GetXaxis().GetNbins())
        self.fReconstructedGenLev.SetTitle("{0} Reconstructed Generator Level".format(self.fName))

        self.fReconstructedDetLev = self.fResponse.ProjectionX("{0}_ReconstructedDetLev".format(self.fName), 1, self.fResponse.GetYaxis().GetNbins())
        self.fReconstructedDetLev.SetTitle("{0} Reconstructed Detector Level".format(self.fName))

        self.GenerateEfficiencies()

        self.GenerateResponseUncertainty()

    def GenerateResponseUncertainty(self):
        self.fResponseUncertainty = ROOT.TH2D("{0}_ResponseMatrixUncertainty".format(self.fName), "{0} Response Matrix Uncertainty".format(self.fName), self.fResponse.GetNbinsX(), self.fResponse.GetXaxis().GetXbins().GetArray(), self.fResponse.GetNbinsY(), self.fResponse.GetYaxis().GetXbins().GetArray())
        self.fResponseUncertainty.GetXaxis().SetTitle(self.fResponse.GetXaxis().GetTitle())
        self.fResponseUncertainty.GetYaxis().SetTitle(self.fResponse.GetYaxis().GetTitle())
        self.fResponseUncertainty.GetZaxis().SetTitle("relative statistical uncertainty")
        for x in range(0, self.fResponseUncertainty.GetNbinsX()+2):
            for y in range(0, self.fResponseUncertainty.GetNbinsY()+2):
                if self.fResponse.GetBinContent(x,y) == 0:
                    continue
                self.fResponseUncertainty.SetBinContent(x,y,self.fResponse.GetBinError(x,y) / self.fResponse.GetBinContent(x,y))

    def GenerateEfficiencies(self):
        # total efficiency
        self.fTotalEfficiency =  self.fReconstructedGenLev.Clone("{0}_TotalEfficiency".format(self.fName))
        self.fTotalEfficiency.Divide(self.fTruth)
        self.fTotalEfficiency.SetTitle("Total Efficiency")
        self.fTotalEfficiency.GetYaxis().SetTitle("Total Efficiency")

        # reconstruction efficiency
        self.fRecoEfficiency =  self.fResponse.ProjectionY("{0}_RecoEfficiency".format(self.fName), 0, -1)
        self.fRecoEfficiency.Divide(self.fTruth)
        self.fRecoEfficiency.SetTitle("Reconstruction Efficiency")
        self.fRecoEfficiency.GetYaxis().SetTitle("Reconstruction Efficiency")
        
        # kinematic efficiency
        self.fKineEfficiency =  self.fReconstructedGenLev.Clone("{0}_KineEfficiency".format(self.fName))
        self.fKineEfficiency.Divide(self.fResponse.ProjectionY("temp", 0, -1))
        self.fKineEfficiency.SetTitle("Kinematic Efficiency")
        self.fKineEfficiency.GetYaxis().SetTitle("Kinematic Efficiency")

    def GenerateRootList(self):
        rlist = ROOT.TList()
        rlist.SetName("Response_{0}".format(self.fName))
        rlist.Add(self.fResponse)
        rlist.Add(self.fNormResponse)
        rlist.Add(self.fResponseUncertainty)
        rlist.Add(self.fTruth)
        rlist.Add(self.fSmallBinResponse)
        rlist.Add(self.fSmallBinTruth)
        rlist.Add(self.fTotalEfficiency)
        rlist.Add(self.fRecoEfficiency)
        rlist.Add(self.fKineEfficiency)
        rlist.Add(self.fReconstructedGenLev)
        rlist.Add(self.fReconstructedDetLev)
        return rlist
