#!/usr/bin/env python
#python script to do extract B feed down correction factors

import argparse
import yaml
import IPython
import ROOT
import DMesonJetUtils

globalList = []

def main():
    ROOT.TH1.AddDirectory(False)
    path = "/Users/sa639/Documents/Work/ALICE/alice-yale-hfjet/fastSimulation"
    generators = ["pythia", "powheg"]
    
    for gen in generators:
        charm_filename = "{0}/FastSimAnalysis_{1}_charm_local.root".format(path, gen)
        beauty_filename = "{0}/FastSimAnalysis_{1}_beauty_local.root".format(path, gen)
        charm_jet, charm_d = GetSpectrum(charm_filename)
        charm_jet.SetTitle("c #rightarrow D^{0}-jet")
        charm_d.SetTitle("c #rightarrow D^{0}")
        beauty_jet, beauty_d = GetSpectrum(beauty_filename)
        beauty_jet.SetTitle("b #rightarrow D^{0}-jet")
        beauty_d.SetTitle("b #rightarrow D^{0}")
        r = DMesonJetUtils.CompareSpectra(charm_jet, [beauty_jet], "BFeedDown_Jet_{0}".format(gen))
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                obj.SaveAs("{0}/FastSimAnalysis_local_{1}.pdf".format(path, obj.GetName()))
        r = DMesonJetUtils.CompareSpectra(charm_d, [beauty_d], "BFeedDown_D_{0}".format(gen))
        for obj in r:
            globalList.append(obj)
            if isinstance(obj, ROOT.TCanvas):
                obj.SaveAs("{0}/FastSimAnalysis_local_{1}.pdf".format(path, obj.GetName()))

def GetSpectrum(filename):
    file = ROOT.TFile(filename)
    jet_list = file.Get("D0_MCTruth_D_Tagged_Jet_Pt_Spectrum")
    jet_h = jet_list.FindObject("D0_MCTruth_D_Tagged_Jet_Pt_Spectrum")
    jet_h.Scale(1, "width")
    jet_h.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [mb (GeV/#it{c})^{-1}]")
    d_list = file.Get("D0_MCTruth_D_Pt_Spectrum")
    d_h = d_list.FindObject("D0_MCTruth_D_Pt_Spectrum")
    d_h.Scale(1, "width")
    d_h.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [mb (GeV/#it{c})^{-1}]")
    return jet_h, d_h

if __name__ == '__main__':
    main()
    
    IPython.embed()