#!/usr/bin/env python

import ROOT
import math
import IPython
import argparse
import DMesonJetCompare
import DMesonJetUtils

globalList = []

def main(tslist):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)

    fileNameTemp = "/Volumes/DATA/ALICE/JetResults/FastSim_pythia_charm_{ts}/FastSimAnalysis_Reduced_pythia_charm_{ts}.root"
    spectraNames = ["D0_MCTruth/Charged_R040/D0_MCTruth_Charged_R040_JetPtSpectrum_DPt_3/D0_MCTruth_Charged_R040_JetPtSpectrum_DPt_3_Normalized",
                    "D0_MCTruth/Charged_R040/D0_MCTruth_Charged_R040_DPtSpectrum_JetPt_5_80/D0_MCTruth_Charged_R040_DPtSpectrum_JetPt_5_80_Normalized",
                    "D0_MCTruth/D0_MCTruth_DPtSpectrum/D0_MCTruth_DPtSpectrum_Normalized"]

    spectra = dict()

    tslistStr = "_".join(tslist)

    for ts in tslist:
        fileName = fileNameTemp.format(ts=ts)
        print("Opening file {0}".format(fileName))
        file = ROOT.TFile(fileName)

        for sname in spectraNames:
            h = DMesonJetUtils.GetObject(file, sname)
            if not h:
                print("Could not find object {} in file {}".format(sname, fileName))
                exit(1)
            hname = h.GetName()
            h.SetName(ts)
            h.SetTitle(ts)
            if not hname in spectra: spectra[hname] = []
            spectra[hname].append(h)

    for hname, hlist in spectra.iteritems():
        print("Histogram: {}".format(hname))
        comp = DMesonJetCompare.DMesonJetCompare(hname)
        comp.fDoRatioPlot = True
        comp.fDoSpectraPlot = "logy"
        comp.fGridyRatio = True
        # comp.fOptSpectrumBaseline = "hist"
        # comp.fOptSpectrum = "hist"
        globalList.extend(hlist)
        r = comp.CompareSpectra(hlist[0], hlist[1:])
        globalList.extend(r)
        comp.fCanvasSpectra.SaveAs("/Volumes/DATA/ALICE/JetResults/{}_{}.pdf".format(hname, tslistStr))
        comp.fCanvasRatio.SaveAs("/Volumes/DATA/ALICE/JetResults/{}_{}_Ratio.pdf".format(hname, tslistStr))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Compare statistical uncertainties.')
    parser.add_argument('ts', metavar='ts',
                        help='Timestamps', nargs='*')
    args = parser.parse_args()

    main(args.ts)

    IPython.embed()
