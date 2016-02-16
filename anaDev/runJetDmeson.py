#!/usr/bin/env python
#python script to test Jet D meson analysis

import argparse
import ROOT
import helperFunctions

def main(fileList, nFiles, nEvents, runPeriod, strmode="AOD", jetRadius=0.6, doHF=True, doChargedJets=True, doTrackQA=True, physSel=ROOT.AliVEvent.kMB, taskName="JetDmesonAna", debugLevel=0):

    ROOT.gSystem.Load("libCGAL")

    ROOT.AliParticleContainer.SetDefTrackCutsPeriod(runPeriod)
    
    mode = None
    if strmode == "AOD":
        mode = helperFunctions.AnaMode.AOD
    elif strmode == "ESD":
        mode = helperFunctions.AnaMode.ESD
    else:
        print "Error: mode has to be either ESD or AOD!"
        exit(1)

    print strmode, "analysis chosen."
    print "Setting local analysis for", nFiles, "files from list", fileList, "max events =", nEvents

    #AliVEvent::kINT7, AliVEvent::kMB, AliVEvent::kCentral, AliVEvent::kSemiCentral
    #AliVEvent::kEMCEGA, AliVEvent::kEMCEJ

    #Analysis manager
    mgr = ROOT.AliAnalysisManager(taskName)

    helperFunctions.LoadMacros()

    if mode is helperFunctions.AnaMode.AOD:
        helperFunctions.AddAODHandler()
    elif mode is AnaMode.ESD:
        helperFunctions.AddESDHandler()

    #Physics selection task
    if mode is helperFunctions.AnaMode.ESD:
        ROOT.AddTaskPhysicsSelection()

    #PID response
    if doHF:
        PIDtask = helperFunctions.AddTaskPIDResponse(False)

    if doTrackQA:
        pSpectraTask = ROOT.AddTaskEmcalJetQA("usedefault", "", "")
        pSpectraTask.GetParticleContainer(0).SetFilterHybridTracks(True)
        pSpectraTask.SetNeedEmcalGeom(False)
        pSpectraTask.SelectCollisionCandidates(physSel)
        pSpectraTask.SetHistoBins(150, 0, 150)
        
    #Charged jet analysis
    if doChargedJets:
        pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", 1, 0.4, 1, 0.15, 0., 0.1, 1, "Jet", 0., False, False, False)
        pChJetTask.SelectCollisionCandidates(physSel)
        pChJetTask.GetParticleContainer(0).SetFilterHybridTracks(True)

        pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", 1, 0.6, 1, 0.15, 0., 0.1, 1, "Jet", 0., False, False, False)
        pChJetTask.SelectCollisionCandidates(physSel)
        pChJetTask.GetParticleContainer(0).SetFilterHybridTracks(True)

        pSpectraTask = ROOT.AddTaskEmcalJetSpectraQA("usedefault", "")
        pSpectraTask.SelectCollisionCandidates(physSel)
        pSpectraTask.SetHistoBins(200, 0, 200)
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kTPCfid)
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.6, ROOT.AliJetContainer.kTPCfid)

    if doHF:
        pDMesonJetsTask = ROOT.AddTaskDmesonJets("usedefault", "");
        pDMesonJetsTask.GetParticleContainer(0).SetFilterHybridTracks(True)
        pDMesonJetsTask.SetShowJetConstituents(True)
        pDMesonJetsTask.SetShowPositionD(True)
        pDMesonJetsTask.SetShowDeltaR(True)
        pDMesonJetsTask.SetShowLeadingPt(True)
        pDMesonJetsTask.SetShowPositionJet(True)
        pDMesonJetsTask.SelectCollisionCandidates(physSel)

        # D*
        pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, 0.4)
        pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, 0.6)
        
        # D0
        pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, 0.4)
        pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, 0.6)
        
    tasks = mgr.GetTasks()
    for task in tasks:
        if isinstance(task, ROOT.AliAnalysisTaskEmcal):
            task.SetForceBeamType(ROOT.AliAnalysisTaskEmcal.kpp)
	
    res = mgr.InitAnalysis()

    if not res:
        print "Error initializing the analysis!"
        exit(1)
    
    mgr.PrintStatus()

    outFile = ROOT.TFile("train.root","RECREATE")
    outFile.cd()
    mgr.Write()
    outFile.Close()
    
    chain = None
    if mode is helperFunctions.AnaMode.AOD:
        ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
        chain = ROOT.CreateAODChain(fileList, nFiles, 0, False, "AliAOD.VertexingHF.root")
    elif mode is helperFunctions.AnaMode.ESD:
        ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
        chain = ROOT.CreateESDChain(fileList, nFiles, 0, False)
        
    if debugLevel == 0:
        mgr.SetUseProgressBar(1, 250)
        
    mgr.SetDebugLevel(debugLevel)

    #To have more debug info
    #pMgr->AddClassDebug("AliEmcalClusTrackMatcherTask", AliLog::kDebug+100);
    
    #start analysis
    print "Starting Analysis..."
    mgr.StartAnalysis("local", chain, nEvents)


if __name__ == '__main__':
    # runDJetCorrAnalysis.py executed as script
    
    parser = argparse.ArgumentParser(description='Jet D meson analysis.')
    parser.add_argument('--mode',
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
    parser.add_argument('-r', '--run-period',
                        help='Run period (e.g. LHC10b)')
    parser.add_argument('--jet-radius',
                        default=0.6,
                        type=float,
                        help='Jet radius')
    parser.add_argument('--no-hf', action='store_const',
                        default=False, const=True,
                        help='No HF analysis')
    parser.add_argument('--no-charged-jets', action='store_const',
                        default=False, const=True,
                        help='No charged jet analysis')
    parser.add_argument('--no-track-qa', action='store_const',
                        default=False, const=True,
                        help='No track QA')
    parser.add_argument('--phys-sel',
                        default=ROOT.AliVEvent.kMB, 
                        help='Physics selection')
    parser.add_argument('--task-name',
                        default="JetDmesonAna",
                        help='Task name')
    parser.add_argument('-d', '--debug-level', 
                        default=0,
                        type=int,
                        help='Debug level')
    args = parser.parse_args()
    
    main(args.fileList, args.n_files, args.n_events, args.run_period, args.mode, args.jet_radius, not args.no_hf, not args.no_charged_jets, not args.no_track_qa, not args.phys_sel, args.task_name, args.debug_level)
