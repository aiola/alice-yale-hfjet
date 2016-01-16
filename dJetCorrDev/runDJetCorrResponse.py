#!/usr/bin/env python
#python script to use the class DJetCorrResponse

import ROOT
import argparse
from ROOT import gROOT
from ROOT import TGaxis

import IPython

import commonFunctions

def main(train, run=True, plot=True, loadLibs=True, inputPath="$JETRESULTS"):
    
    ROOT.TH1.AddDirectory(False)
    
    TGaxis.SetMaxDigits(3)
  
    tracksD0Name = "DcandidatesAndTracksD0"
    tracksDStarName = "DcandidatesAndTracksDStar"
    
    if loadLibs:
        commonFunctions.LoadDJetCorrClasses()

    from ROOT import DJetCorrResponse
  
    projDjet = DJetCorrResponse(train);

    #projDjet->SetEfficiencyMode("TH1")
    projDjet.SetOverwrite(True)
    projDjet.SetInputPath(inputPath)
        
    projDjet.SetPlotFormat("pdf")
    projDjet.SetSavePlots(False)
        
    #projDjet.AddAnalysisParams("D0", "Full", "R040", tracksD0Name)
    projDjet.AddAnalysisParams("D0", "Full", "R060", tracksD0Name)
        
    #projDjet.AddAnalysisParams("DStar", "Charged", "R040", tracksDStarName)
    projDjet.AddAnalysisParams("DStar", "Charged", "R060", tracksDStarName)
          
    if run:
        projDjet.ProjectResponseMatrices()
        projDjet.SaveOutputFile()
        
    if plot:
        projDjet.PlotResponseMatrices()
        
    return projDjet;

if __name__ == '__main__':
    # runDJetCorrResponse.py executed as script
    
    parser = argparse.ArgumentParser(description='D meson jet correlation response matrix.')
    parser.add_argument('train', metavar='train',
                        help='Train to be analyzed')
    parser.add_argument('-r', '--run', action='store_const',
                        default=False, const=True,
                        help='Run the analysis')
    parser.add_argument('-p', '--plot', action='store_const',
                        default=False, const=True,
                        help='Plot the results')
    parser.add_argument('--no-Libs', action='store_const',
                        default=False, const=True,
                        help='Load the DJetCorr libraries')
    parser.add_argument('--inputPath', metavar='inputPath',
                        default="$JETRESULTS",
                        help='Input path')
    args = parser.parse_args()
    
    main(args.train, args.run, args.plot, not args.no_Libs, args.inputPath)
    
    IPython.embed()
