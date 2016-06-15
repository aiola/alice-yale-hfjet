#!/usr/bin/env python
#python script to test the class MassFitter

import ROOT
import IPython
import subprocess

globalList = []

def TestMassFitter():
  ROOT.gRandom = ROOT.TRandom3(0)
  
  subprocess.call("make")
  ROOT.gSystem.Load("MassFitter.so")

  mpi = ROOT.TDatabasePDG.Instance().GetParticle(211).Mass()

  fitter = ROOT.MassFitter("MyMassFitter", ROOT.MassFitter.kGaus, ROOT.MassFitter.kExpo, 1.7, 2)
  
  fitter.GetFitFunction().SetParameter(0,   100.00)
  fitter.GetFitFunction().SetParameter(1,    -2.00)
  fitter.GetFitFunction().SetParameter(2,    60.00)
  fitter.GetFitFunction().SetParameter(3,     1.87)
  fitter.GetFitFunction().SetParameter(4,     0.01)
  
  c = ROOT.TCanvas()
  hist = ROOT.TH1F("hist", "hist", 100, 1.7, 2)
  hist.FillRandom(fitter.GetFitFunction().GetName(),100)
  hist.Sumw2()
  hist.Draw()

  fitter.GetFitFunction().SetParameter(0, 100)
  fitter.GetFitFunction().SetParameter(1, -5)

  fitter.Fit(hist, "0 E S")

  fitter.Draw("same")

  histInt = hist.Integral(hist.GetXaxis().FindBin(1.7), hist.GetXaxis().FindBin(2))
  fitInt = fitter.GetFitFunction().Integral(1.7,2)/0.3*100;

  print("hist integral = {0:.3f}, fit integral = {1:.3f}".format(histInt, fitInt))

  bkg = fitter.GetBackground()
  sig = fitter.GetSignal()

  print("bkg = {0:.3f}, sig = {1:.3f}".format(bkg, sig))
  
  globalList.append(hist)
  globalList.append(c)
  globalList.append(fitter)

if __name__ == '__main__':

    TestMassFitter()
    
    IPython.embed()
    