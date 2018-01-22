#!/usr/bin/env python
# python script to test the class MassFitter

import ROOT
import IPython
import subprocess
import math

globalList = []

# To mimic ROOT5 behavior
if ROOT.gROOT.GetVersionInt() >= 60000: ROOT.ROOT.Math.IntegratorOneDimOptions.SetDefaultIntegrator("Gauss")


def TestMassFitter():
    ROOT.TH1.AddDirectory(False)
    ROOT.gRandom = ROOT.TRandom3(0)
    subprocess.call("make")
    ROOT.gSystem.Load("MassFitter.so")

    invMassHist = GetMassHistogram()

    # some defualt params
    pdgMass = ROOT.TDatabasePDG.Instance().GetParticle(421).Mass()
    massWidth = 0.012

    minMass = invMassHist.GetXaxis().GetBinLowEdge(1)
    maxMass = invMassHist.GetXaxis().GetBinUpEdge(invMassHist.GetXaxis().GetNbins())
    totIntegral = invMassHist.Integral(1, invMassHist.GetXaxis().GetNbins(), "width")
    integral3sigma = invMassHist.Integral(invMassHist.GetXaxis().FindBin(pdgMass - 3 * massWidth), invMassHist.GetXaxis().FindBin(pdgMass + 3 * massWidth), "width")

    expoParBkg1 = math.log(invMassHist.GetBinContent(invMassHist.GetXaxis().GetNbins()) / invMassHist.GetBinContent(1)) / (maxMass - minMass)
    expoParBkg0 = totIntegral

    sig = integral3sigma - (math.exp(expoParBkg1 * (pdgMass - 3 * massWidth)) - math.exp(expoParBkg1 * (pdgMass + 3 * massWidth))) / (-expoParBkg1) * expoParBkg0
    GaussConst = sig

    refl_function = ROOT.TF1("refl_function", "[p0]/(TMath::Sqrt(2.*TMath::Pi())*[p2])*TMath::Exp(-(x-[p1])*(x-[p1])/(2.*[p2]*[p2]))+[p3]/(TMath::Sqrt(2.*TMath::Pi())*[p5])*TMath::Exp(-(x-[p4])*(x-[p4])/(2.*[p5]*[p5]))", minMass, maxMass)
    refl_function.SetParameter(0, 0.066022)
    refl_function.SetParameter(1, 0.043864)
    refl_function.SetParameter(2, 0.519630)
    refl_function.SetParameter(3, 0.013195)
    refl_function.SetParameter(4, 1.856881)
    refl_function.SetParameter(5, 0.088874)
    refl_hist = ROOT.TH1D("refl_hist", "refl_hist", invMassHist.GetXaxis().GetNbins(), minMass, maxMass)
    refl_hist.Add(refl_function)

    fitter = ROOT.MassFitter("MyMassFitter", ROOT.MassFitter.kDzeroKpi, minMass, maxMass)
    fitter.SetFitRange(minMass, maxMass)
    fitter.SetHistogram(invMassHist)
    fitter.SetReflectionTemplate(refl_hist, 0.3)
    fitter.GetFitFunction().SetParameter(0, expoParBkg0)
    fitter.GetFitFunction().SetParameter(1, expoParBkg1)
    fitter.GetFitFunction().SetParameter(2, GaussConst)
    fitter.GetFitFunction().SetParameter(3, pdgMass)  # start fitting using PDG mass
    fitter.GetFitFunction().SetParLimits(3, pdgMass * 0.9, pdgMass * 1.1)  # limiting mass parameter +/- 10% of PDG value
    fitter.GetFitFunction().SetParameter(4, massWidth)
    fitter.GetFitFunction().SetParLimits(4, 0, 1)  # limiting width to being positive

    # fitter.DisableRefl()

    c = ROOT.TCanvas()
    invMassHist.Draw()

    print("Expo index before fit: {0:.3f}".format(expoParBkg1))
    print("Gauss const before fit: {0:.3f}".format(GaussConst))

    fitInt = fitter.GetFitFunction().Integral(minMass, maxMass)
    print("before fitting: sig = {0:.3f}, hist tot integral = {1:.3f}, fit integral = {2:.3f}".format(sig, totIntegral, fitInt))

    fitter.Fit("0 L S V")
    fitInt = fitter.GetFitFunction().Integral(minMass, maxMass)
    print("after fitting: hist integral = {0:.3f}, fit integral = {1:.3f}".format(totIntegral, fitInt))

    fitter.Draw("same")

    bkg = fitter.GetBackground(3)
    sig = fitter.GetSignal()

    print("bkg = {0:.3f}, sig = {1:.3f}".format(bkg, sig))

    globalList.append(invMassHist)
    globalList.append(c)
    globalList.append(fitter)


def GetMassHistogram():
    file = ROOT.TFile("MassFitterTest.root")
    if not file or file.IsZombie():
        print("Could not open file MassFitterTest.root")
        exit(1)
    h = file.Get("hMass")
    if not h:
        print("Could not find histogram hMass")
        file.ls()
        exit(1)
    h_copy = h.Clone()
    file.Close()
    return h_copy


if __name__ == '__main__':
    TestMassFitter()

    IPython.embed()
