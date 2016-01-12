#!/usr/bin/env ipython -i
#python script to use the class DJetCorrUnfold

import ROOT
import argparse
from ROOT import gROOT
from ROOT import TGaxis

import IPython

import commonFunctions
import runDJetCorrAnalysis
import runDJetCorrResponse

def runDJetCorrUnfold(trainData, trainResp, loadLibs = True, inputPath = "$JETRESULTS"):	
	if loadLibs:
		commonFunctions.LoadDJetCorrClasses()

	from ROOT import DJetCorrAnalysis
	from ROOT import DJetCorrResponse
	from ROOT import DJetCorrUnfold
  
  	ana = runDJetCorrAnalysis.runDJetCorrAnalysis(trainData, False, False, False, False, True, False, False, inputPath)
   	resp = runDJetCorrResponse.runDJetCorrResponse(trainResp, False, False, False, inputPath)
  
  	unfold = DJetCorrUnfold(ana, resp)
   	unfold.SetDataParamIndex(0)
   	unfold.SetRespParamIndex(0)
   	unfold.SetSavePlots(True)
    #unfold.SetUseEfficiency(False)
    #unfold.SetUseKinEfficiency(False)
	unfold.Start()

	return unfold

if __name__ == '__main__':
    # runDJetCorrUnfold.py executed as script
    
    parser = argparse.ArgumentParser(description='D meson jet correlation response matrix.')
    parser.add_argument('trainData', metavar='trainData',
                        help='Train with data to be analyzed')
    parser.add_argument('trainResp', metavar='trainResp',
						help='MC train used for the response matrix')
    parser.add_argument('--no-Libs', action='store_const',
                        default=False, const=True,
                        help='Load the DJetCorr libraries')
    parser.add_argument('--inputPath', metavar='inputPath',
                        default="$JETRESULTS",
                        help='Input path')
    args = parser.parse_args()
    
    runDJetCorrUnfold(args.trainData, args.trainResp, not args.no_Libs, args.inputPath)
