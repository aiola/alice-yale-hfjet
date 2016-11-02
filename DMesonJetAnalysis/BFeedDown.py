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
    ROOT.gStyle.SetOptTitle(0)
    path = "/Users/sa639/Documents/Work/ALICE/alice-yale-hfjet/fastSimulation"
    generators = ["powheg"]
    spectra = ["D0_MCTruth_D_Tagged_Jet_Pt_Spectrum", "D0_MCTruth_D_Pt_Spectrum", "D0_MCTruth_D_Tagged_Jet_Pt_Spectrum_PtD_2", "D0_MCTruth_D_Pt_Spectrum_JetPt_5"]
    
    for gen in generators:
        charm_filename = "{0}/FastSimAnalysis_{1}_charm_1478079449_local.root".format(path, gen)
        beauty_filename = "{0}/FastSimAnalysis_{1}_beauty_1478079771_local.root".format(path, gen)
        for spectrum in spectra:
            charm = GetSpectrum(charm_filename, spectrum)
            charm.SetTitle("c #rightarrow D^{0}")
            beauty = GetSpectrum(beauty_filename, spectrum)
            beauty.SetTitle("b #rightarrow D^{0}")
            r = DMesonJetUtils.CompareSpectra(charm, [beauty], "BFeedDown_{0}_{1}".format(gen, spectrum))
            for obj in r:
                globalList.append(obj)
                if isinstance(obj, ROOT.TCanvas):
                    obj.SaveAs("{0}/FastSimAnalysis_local_{1}.pdf".format(path, obj.GetName()))

def GetSpectrum(filename, spectrum):
    file = ROOT.TFile(filename)
    rlist = file.Get(spectrum)
    h = rlist.FindObject(spectrum)
    h.Scale(1, "width")
    h.GetYaxis().SetTitle("#frac{d#sigma}{d#it{p}_{T}} [mb (GeV/#it{c})^{-1}]")
    return h

if __name__ == '__main__':
    main()
    
    IPython.embed()
