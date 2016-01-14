#!/usr/bin/env python
#python script to use the class DJetCorrAnalysis

import ROOT
import argparse
from ROOT import gROOT
from ROOT import TGaxis

import IPython

import commonFunctions

def main(train, run=True, refit=False, plot=True, truthOnly=False, isMC=False, isBkgSub=False, loadLibs=True, inputPath="$JETRESULTS"):
    
    ROOT.TH1.AddDirectory(False)
    
    TGaxis.SetMaxDigits(3) 

    tracksName = "tracks"
    tracksD0Name = "DcandidatesAndTracksD0"
    tracksDStarName = "DcandidatesAndTracksDStar"

    if not isMC:
        tracksD0Name += "rec"
        tracksDStarName += "rec"

    if loadLibs:
        commonFunctions.LoadDJetCorrClasses()
  
    from ROOT import DJetCorrAnalysis
  
    projDjet = DJetCorrAnalysis(train);
  
    projDjet.SetOverwrite(True)
    projDjet.SetInputPath(inputPath)

    qaListName = "AliAnalysisTaskSAQA_{0}_TPC_histos".format(tracksName)
  
    projDjet.SetQAListName(qaListName)

    projDjet.SetPlotFormat("pdf")
    projDjet.SetSavePlots(True)

    #projDjet.AddAnalysisParams("D0", "Full", "R040", tracksD0Name, isMC, isBkgSub)
    #projDjet.AddAnalysisParams("DStar", "Charged", "R040", tracksDStarName, isMC, isBkgSub)

    param = projDjet.AddAnalysisParams("D0", "Full", "R060", tracksD0Name, isMC, isBkgSub)
    #param.BackgroundOnly(True);
    #param.SignalOnly(True);

    param = projDjet.AddAnalysisParams("DStar", "Charged", "R060", tracksDStarName, isMC, isBkgSub)
    #param.BackgroundOnly(True);
    #param.SignalOnly(True);

    if truthOnly:
         projDjet.ProjectTruthSpectrum()
  
    else:
        if run:
            projDjet.GenerateQAHistograms()
            projDjet.GenerateDJetCorrHistograms()
            projDjet.ProjectTruthSpectrum()

        if refit:
            projDjet.PlotTrackHistograms()
            projDjet.PlotDJetCorrHistograms(True)
        elif plot:
            projDjet.PlotTrackHistograms()
            projDjet.PlotDJetCorrHistograms(False)

    if truthOnly or run or refit or plot: 
        projDjet.SaveOutputFile()

    return projDjet

if __name__ == '__main__':
    # runDJetCorrAnalysis.py executed as script
    
    parser = argparse.ArgumentParser(description='D meson jet correlation response matrix.')
    parser.add_argument('train', metavar='train',
                        help='Train to be analyzed')
    parser.add_argument('-r', '--run', action='store_const',
                        default=False, const=True,
                        help='Run the analysis')
    parser.add_argument('--refit', action='store_const',
                        default=False, const=True,
                        help='Refit using previously projected his')
    parser.add_argument('-p', '--plot', action='store_const',
                        default=False, const=True,
                        help='Plot the results')
    parser.add_argument('--truth-only', action='store_const',
                        default=False, const=True,
                        help='Only work on truth histograms')
    parser.add_argument('--MC', action='store_const',
                        default=False, const=True,
                        help='Monte Carlo train')
    parser.add_argument('--bkg-sub', action='store_const',
                        default=False, const=True,
                        help='Monte Carlo background subtracted')
    parser.add_argument('--no-Libs', action='store_const',
                        default=False, const=True,
                        help='Load the DJetCorr libraries')
    parser.add_argument('--inputPath', metavar='inputPath',
                        default="$JETRESULTS",
                        help='Input path')
    args = parser.parse_args()
    
    main(args.train, args.run, args.refit, args.plot, args.truth_only, args.MC, args.bkg_sub, not args.no_Libs, args.inputPath)
    
    IPython.embed()
