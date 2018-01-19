#!/usr/bin/env python
# python script to test the class MassFitter

import ROOT
import IPython
import subprocess
import math

globalList = []


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
    expoParBkg0 = totIntegral / (math.exp(expoParBkg1 * minMass) - math.exp(expoParBkg1 * maxMass)) * (-expoParBkg1)

    sig = integral3sigma - (math.exp(expoParBkg1 * (pdgMass - 3 * massWidth)) - math.exp(expoParBkg1 * (pdgMass + 3 * massWidth))) / (-expoParBkg1) * expoParBkg0
    GaussConst = sig / math.sqrt(2 * math.pi) / massWidth

    fitter = ROOT.MassFitter("MyMassFitter", ROOT.MassFitter.kDzeroKpi, minMass, maxMass)
    fitter.SetFitRange(minMass, maxMass)
    fitter.SetHistogram(invMassHist)
    fitter.GetFitFunction().SetParameter(0, expoParBkg0)
    fitter.GetFitFunction().SetParameter(1, expoParBkg1)
    fitter.GetFitFunction().SetParameter(2, GaussConst)
    fitter.GetFitFunction().SetParameter(3, pdgMass)  # start fitting using PDG mass
    fitter.GetFitFunction().SetParLimits(3, pdgMass * 0.9, pdgMass * 1.1)  # limiting mass parameter +/- 10% of PDG value
    fitter.GetFitFunction().SetParameter(4, massWidth)
    fitter.GetFitFunction().SetParLimits(4, 0, 1)  # limiting width to being positive

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
