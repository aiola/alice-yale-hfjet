#!/usr/bin/env python
#python script to use the class DJetCorrAnalysis

import ROOT
import argparse
from ROOT import gROOT
from ROOT import TGaxis

import IPython

import commonFunctions

def main(train, run=True, refit=False, plot=True,
         D0Ana=True, DstarAna=True, ChargedAna=True, FullAna=True, radius="R040",
         tracksName = "tracks", caloName = "caloClusters", cellName = "emcalCells", 
         loadLibs=True, inputPath="$JETRESULTS"):
    
    ROOT.TH1.AddDirectory(False)
    
    TGaxis.SetMaxDigits(3) 

    if loadLibs:
        commonFunctions.LoadDJetCorrClasses()
  
    from ROOT import DJetCorrAnalysis
  
    projDjet = DJetCorrAnalysis(train);
  
    projDjet.SetOverwrite(True)
    projDjet.SetInputPath(inputPath)
    
    collName = "_"
    
    if tracksName:
        collName += tracksName + "_"
    if caloName:
        collName += caloName + "_"
    if cellName:
        collName += cellName + "_"

    qaListName = "AliAnalysisTaskEmcalJetQA{0}_AnyINT_histos".format(collName)
  
    #projDjet.SetQAListName(qaListName)

    projDjet.SetPlotFormat("pdf")
    projDjet.SetSavePlots(True)
    projDjet.SetInputFileName("AnalysisResults_DJets.root")


    if D0Ana:
        if FullAna:
            param = projDjet.AddAnalysisParams("D0", "Full", radius, "AnyINT")
        if ChargedAna:
            param = projDjet.AddAnalysisParams("D0", "Charged", radius, "AnyINT")
     
    if DstarAna:
        if FullAna:
            param = projDjet.AddAnalysisParams("DStar", "Full", radius, "AnyINT")
        if ChargedAna:
            param = projDjet.AddAnalysisParams("DStar", "Charged", radius, "AnyINT")   

    if run:
        projDjet.GenerateQAHistograms()
        projDjet.GenerateDJetCorrHistograms()

    if refit:
        projDjet.PlotTrackHistograms()
        projDjet.PlotDJetCorrHistograms(True)
    elif plot:
        projDjet.PlotTrackHistograms()
        projDjet.PlotDJetCorrHistograms(False)

    if run or refit or plot: 
        projDjet.SaveOutputFile()

    return projDjet

if __name__ == '__main__':
    # runDJetCorrAnalysis.py executed as script
    
    parser = argparse.ArgumentParser(description='D meson jet correlation response matrix.')
    parser.add_argument('train', metavar='train',
                        help='Train to be analyzed')
    parser.add_argument('--run', action='store_const',
                        default=False, const=True,
                        help='Run the analysis')
    parser.add_argument('--refit', action='store_const',
                        default=False, const=True,
                        help='Refit using previously projected his')
    parser.add_argument('--plot', action='store_const',
                        default=False, const=True,
                        help='Plot the results')
    parser.add_argument('--charged', action='store_const',
                        default=False, const=True,
                        help='Charged jet analysis')
    parser.add_argument('--full', action='store_const',
                        default=False, const=True,
                        help='Full jet analysis')
    parser.add_argument('--D0', action='store_const',
                        default=False, const=True,
                        help='D0 analysis')
    parser.add_argument('--DStar', action='store_const',
                        default=False, const=True,
                        help='D* analysis')
    parser.add_argument('--radius', metavar='radius',
                        default="R040",
                        help='Jet radius (e.g. R040)')
    parser.add_argument('--no-Libs', action='store_const',
                        default=False, const=True,
                        help='Load the DJetCorr libraries')
    parser.add_argument('--inputPath', metavar='inputPath',
                        default="$JETRESULTS",
                        help='Input path')
    parser.add_argument('--tracks', metavar='tracks',
                        default="tracks",
                        help='Track collection name')
    parser.add_argument('--calo', metavar='calo',
                        default="caloClusters",
                        help='Calorimeter cluster collection name')
    parser.add_argument('--cells', metavar='cells',
                        default="emcalCells",
                        help='Calorimeter cell collection name')
    args = parser.parse_args()
    
    main(args.train, args.run, args.refit, args.plot,
         args.D0, args.DStar, args.charged, args.full, args.radius,
         args.tracks, args.calo, args.cells, 
         not args.no_Libs, args.inputPath)
    
    IPython.embed()
