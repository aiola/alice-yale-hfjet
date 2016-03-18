#!/usr/bin/env python
#python script to test EMCal trigger QA analysis

import argparse
import ROOT
import helperFunctions
from enum import Enum

def runEMCalJetAnalysisNew(fileList, nFiles, nEvents, runPeriod, strmode="AOD", taskName="JetAna",
						chargedJets=False, fullJets=False, mcJets=False, 
						local=True, gridMode="full", debugLevel=0):
	ROOT.gSystem.Load("libCGAL")
	
	#AliEmcalPhysicsSelection::kEmcalOk, AliEmcalPhysicsSelection::kEmcalH,
	#AliVEvent::kINT7, AliVEvent::kMB, AliVEvent::kCentral, AliVEvent::kSemiCentral,
	#AliVEvent::kEMCEGA, AliVEvent::kEMCEJE
	#AliVEvent::kAnyINT
	physSel = ROOT.AliVEvent.kMB

	if strmode == "ESD":
		dataType = helperFunctions.AnaMode.ESD
	elif strmode == "AOD":
		dataType = helperFunctions.AnaMode.AOD
	else:
		print("Incorrect data type option, check third argument of run macro.")
		printf("datatype = AOD or ESD")
		exit(1)
  
	if local:
		anaType = 0
		OCDBpath = "local:///Volumes/DATA/ALICE/OCDB/2011"
	else:
		anaType = 2
		OCDBpath = "raw://"
		
  	ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/runEMCalJetAnalysisNew.C")
  	pMgr = ROOT.runEMCalJetAnalysisNew(strmode, fileList, nFiles, nEvents, runPeriod, physSel, taskName, chargedJets, fullJets, mcJets, OCDBpath, anaType, gridMode)

	if local:
		if debugLevel == 0:
			pMgr.SetUseProgressBar(True, 25)
    
		else:
			pMgr.SetUseProgressBar(False, 25)
			pMgr.SetDebugLevel(debugLevel)
  
      	#pMgr.AddClassDebug("AliEmcalJetTask", AliLog::kDebug+100)
      	chain = None
      	if dataType is helperFunctions.AnaMode.AOD:
      		ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C")
      		chain = ROOT.CreateAODChain(fileList, nFiles, 0, False)
      	else:
			ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C")
			chain = ROOT.CreateESDChain(fileList, nFiles, 0, False)

	    #start analysis
    	print("Starting Analysis...")
    	pMgr.StartAnalysis("local", chain, nEvents)

if __name__ == '__main__':
   
    parser = argparse.ArgumentParser(description='EMCal jet analysis.')
    parser.add_argument('--mode', metavar='mode',
                        default='AOD',
                        help='Analysis mode (ESD or AOD)')
    parser.add_argument('fileList', metavar='fileList',
                        help='File list to be analyzed')
    parser.add_argument('-f', '--n-files', 
                        default=100,
                        type=int,
                        help='Number of files to be analyzed')
    parser.add_argument('-e', '--n-events',
                        type=int,
                        help='Number of events to be analyzed')
    parser.add_argument('-r', '--run-period', metavar='run_period',
                        help='Run period (e.g. LHC10b)')
    parser.add_argument('--task-name', metavar='task_name',
                        default="JetAna",
                        help='Task name')
    parser.add_argument('-d', '--debug-level', 
                        default=0,
                        type=int,
                        help='Debug level')
    parser.add_argument('--grid-mode',
						default="",
						help='Grid mode ("submit", "full", "terminate", etc.)')
    parser.add_argument('--charged-jets', action='store_const',
						default=False, const=True,
						help='Charged jet analysis')
    parser.add_argument('--full-jets', action='store_const',
						default=False, const=True,
						help='Full jet analysis')
    parser.add_argument('--mc-jets', action='store_const',
						default=False, const=True,
						help='MC jet analysis')
    args = parser.parse_args()
    
    if args.grid_mode:
    	local = False
    else:
    	local = True
    
    runEMCalJetAnalysisNew(args.fileList, args.n_files, args.n_events, args.run_period, args.mode, args.task_name,
						args.charged_jets, args.full_jets, args.mc_jets,
						local, args.grid_mode,
						args.debug_level)
