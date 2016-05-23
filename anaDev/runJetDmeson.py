#!/usr/bin/env python
#python script to test Jet D meson analysis

import argparse
import ROOT
import helperFunctions

def main(fileList, nFiles, nEvents, runPeriod, strmode="AOD", doHF=True, doChargedJets=False, doFullJets=False, doTrackQA=True, 
         taskName="JetDmesonAna", debugLevel=0):

    physSel = 0
    ROOT.gSystem.Load("libCGAL")

    ROOT.AliTrackContainer.SetDefTrackCutsPeriod(runPeriod)
    
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
    elif mode is helperFunctions.AnaMode.ESD:
        helperFunctions.AddESDHandler()
        
    #Physics selection task
    if mode is helperFunctions.AnaMode.ESD:
        ROOT.AddTaskPhysicsSelection()

    #Setup task
    if doFullJets or mode is helperFunctions.AnaMode.ESD:
        OCDBpath = "local:///Volumes/DATA/ALICE/OCDB/2012";
        pSetupTask = ROOT.AliEmcalSetupTask("EmcalSetupTask");
        pSetupTask.SetNoOCDB(0)
        mgr.AddTask(pSetupTask)
        cinput = mgr.GetCommonInputContainer()
        mgr.ConnectInput(pSetupTask, 0,  cinput)
        pSetupTask.SetOcdbPath(OCDBpath)
        
    if doFullJets:
        helperFunctions.PrepareEMCAL(physSel)

    #PID response
    if doHF:
        helperFunctions.AddTaskPIDResponse(False)

    if doTrackQA:
        if doChargedJets and not doFullJets:
            pSpectraTask = ROOT.AddTaskEmcalJetQA("usedefault", "", "")
            pSpectraTask.SetNeedEmcalGeom(False)
        if doFullJets:
            pSpectraTask = ROOT.AddTaskEmcalJetQA("usedefault", "usedefault", "usedefault")
            pSpectraTask.SetNeedEmcalGeom(True)
        pSpectraTask.SelectCollisionCandidates(physSel)
        pSpectraTask.SetPtBin(1, 150)

    #Charged jet analysis
    if doChargedJets:
        pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", 1, 0.2, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pChJetTask.SelectCollisionCandidates(physSel)

        pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", 1, 0.4, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pChJetTask.SelectCollisionCandidates(physSel)

    #Full jet analysis
    if doFullJets:
        pJetTask = ROOT.AddTaskEmcalJet("usedefault", "usedefault", 1, 0.2, ROOT.AliJetContainer.kFullJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pJetTask.SelectCollisionCandidates(physSel)

        pJetTask = ROOT.AddTaskEmcalJet("usedefault", "usedefault", 1, 0.4, ROOT.AliJetContainer.kFullJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pJetTask.SelectCollisionCandidates(physSel)

    pSpectraTask = ROOT.AddTaskEmcalJetTree("usedefault", "usedefault")
    pSpectraTask.SelectCollisionCandidates(physSel)
    pSpectraTask.SetHistoType(ROOT.AliAnalysisTaskEmcalJetSpectraQA.kTTree)

    if doChargedJets:
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.2, ROOT.AliJetContainer.kTPCfid)
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kTPCfid)

    if doFullJets:
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.2, ROOT.AliJetContainer.kEMCALfid)
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kEMCALfid)

    if doHF:
        pDMesonJetsTask = ROOT.AddTaskDmesonJets("usedefault", "usedefault");
        pDMesonJetsTask.SetShowJetConstituents(True)
        pDMesonJetsTask.SetShowPositionD(True)
        pDMesonJetsTask.SetShowDeltaR(True)
        pDMesonJetsTask.SetShowPositionJet(True)
        pDMesonJetsTask.SelectCollisionCandidates(physSel)
        pDMesonJetsTask.SetTreeOutput(True)
        
        # D0
        if doChargedJets:
            pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.2)
            pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.4)

        if doFullJets:
            pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.2)
            pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.4)

        # D*
        if doChargedJets:
            pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.2)
            pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.4)

        if doFullJets:
            pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.2)
            pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.4)
        
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
    #mgr.AddClassDebug("AliAnalysisTaskDmesonJets", ROOT.AliLog.kDebug+100)
    
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
                        default=2000,
                        help='Number of events to be analyzed')
    parser.add_argument('-r', '--run-period',
                        help='Run period (e.g. LHC10b)')
    parser.add_argument('--no-hf', action='store_const',
                        default=False, const=True,
                        help='No HF analysis')
    parser.add_argument('--charged-jets', action='store_const',
                        default=False, const=True,
                        help='Charged jet analysis')
    parser.add_argument('--full-jets', action='store_const',
                        default=False, const=True,
                        help='Full jet analysis')
    parser.add_argument('--no-track-qa', action='store_const',
                        default=False, const=True,
                        help='No track QA')
    parser.add_argument('--task-name',
                        default="JetDmesonAna",
                        help='Task name')
    parser.add_argument('-d', '--debug-level', 
                        default=0,
                        type=int,
                        help='Debug level')
    args = parser.parse_args()
    
    main(args.fileList, args.n_files, args.n_events, args.run_period, args.mode, 
         not args.no_hf, args.charged_jets, args.full_jets, not args.no_track_qa, 
         args.task_name, args.debug_level)
