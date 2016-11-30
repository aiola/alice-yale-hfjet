#!/usr/bin/env python
#python script to test Jet D meson analysis

import argparse
import ROOT
import helperFunctions
import yaml

def main(configFileName, nFiles, nEvents, doRecLevel, doSignalOnly, doMCTruth, doResponse, noInclusiveJets,
         taskName="JetDmesonAna", debugLevel=0):
    
    f = open(configFileName, 'r')
    config = yaml.load(f)
    f.close()

    physSel = 0
    ROOT.gSystem.Load("libCGAL")

    ROOT.AliTrackContainer.SetDefTrackCutsPeriod(config["run_period"])
    
    mode = None
    if config["mode"] == "AOD":
        mode = helperFunctions.AnaMode.AOD
    elif config["mode"] == "ESD":
        mode = helperFunctions.AnaMode.ESD
    else:
        print "Error: mode has to be either ESD or AOD!"
        exit(1)

    print("{0} analysis chosen.".format(config["mode"]))
    print("Setting local analysis for {0} files from list {1} max events = {2}".format(nFiles, config["file_list"], nEvents))

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
    if mode is helperFunctions.AnaMode.ESD and not config["MC"]:
        ROOT.AddTaskPhysicsSelection()

    #Setup task
    if config["full_jets"] or mode is helperFunctions.AnaMode.ESD:
        OCDBpath = "raw://";
        pSetupTask = ROOT.AliEmcalSetupTask("EmcalSetupTask");
        pSetupTask.SetNoOCDB(0)
        mgr.AddTask(pSetupTask)
        cinput = mgr.GetCommonInputContainer()
        mgr.ConnectInput(pSetupTask, 0,  cinput)
        pSetupTask.SetOcdbPath(OCDBpath)
        
    if config["full_jets"]:
        helperFunctions.PrepareEMCAL("userQAconfiguration.yaml")

    if config["MC"]:
        print "Running on a MC production"
    else:
        print "Running on data"

    #PID response
    helperFunctions.AddTaskPIDResponse(config["MC"], True, True, config["reco_pass"])

    if config["full_jets"]:
        pSpectraTask = ROOT.AddTaskEmcalJetQA("usedefault", "usedefault", "usedefault")
        pSpectraTask.SetNeedEmcalGeom(True)
    else:
        pSpectraTask = ROOT.AddTaskEmcalJetQA("usedefault", "", "")
        pSpectraTask.SetNeedEmcalGeom(False)

    pSpectraTask.SelectCollisionCandidates(physSel)
    pSpectraTask.SetPtBin(1, 150)

    if not noInclusiveJets:
        #Charged jet analysis
        if config["charged_jets"]:
            pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", 1, 0.4, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
            pChJetTask.SelectCollisionCandidates(physSel)
    
            pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", 1, 0.6, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
            pChJetTask.SelectCollisionCandidates(physSel)
    
        #Full jet analysis
        if config["full_jets"]:
            pJetTask = ROOT.AddTaskEmcalJet("usedefault", "usedefault", 1, 0.2, ROOT.AliJetContainer.kFullJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
            pJetTask.SelectCollisionCandidates(physSel)
    
            pJetTask = ROOT.AddTaskEmcalJet("usedefault", "usedefault", 1, 0.4, ROOT.AliJetContainer.kFullJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
            pJetTask.SelectCollisionCandidates(physSel)
    
        if config["full_jets"]:
            pJetSpectraTask = ROOT.AddTaskEmcalJetTree("usedefault", "usedefault")
            pJetSpectraTask.SetNeedEmcalGeom(True)
        else:
            pJetSpectraTask = ROOT.AddTaskEmcalJetTree("usedefault", "")
            pJetSpectraTask.SetNeedEmcalGeom(False)
                    
        pJetSpectraTask.SelectCollisionCandidates(physSel)
    
        if config["charged_jets"]:
            pJetSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kTPCfid)
            pJetSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.6, ROOT.AliJetContainer.kTPCfid)
    
        if config["full_jets"]:
            pJetSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.2, ROOT.AliJetContainer.kEMCALfid)
            pJetSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kEMCALfid)

    nOutputTrees = 0

    if doRecLevel:
        nOutputTrees += 3
        
    if config["MC"] and doSignalOnly:
        nOutputTrees += 2
        
    if config["MC"] and doMCTruth:
        nOutputTrees += 2
        
    if nOutputTrees > 0:
        if config["MC"]:
            if doResponse:
                if config["full_jets"]:
                    pDMesonJetsTask = ROOT.AddTaskDmesonJetsDetectorResponse("usedefault", "usedefault", "usedefault", 2)
                    pDMesonJetsTask.SetNeedEmcalGeom(True)
                else:
                    pDMesonJetsTask = ROOT.AddTaskDmesonJetsDetectorResponse("usedefault", "", "usedefault", 2)
            else:
                if config["full_jets"]:
                    pDMesonJetsTask = ROOT.AddTaskDmesonJets("usedefault", "usedefault", "usedefault", nOutputTrees)
                    pDMesonJetsTask.SetNeedEmcalGeom(True)
                else:
                    pDMesonJetsTask = ROOT.AddTaskDmesonJets("usedefault", "", "usedefault", nOutputTrees)
                pDMesonJetsTask.SetOutputType(ROOT.AliAnalysisTaskDmesonJets.kTreeOutput)
        else:
            if config["full_jets"]:
                pDMesonJetsTask = ROOT.AddTaskDmesonJets("usedefault", "usedefault", "", nOutputTrees)
                pDMesonJetsTask.SetNeedEmcalGeom(True)
            else:
                pDMesonJetsTask = ROOT.AddTaskDmesonJets("usedefault", "", "", nOutputTrees)

        pDMesonJetsTask.SelectCollisionCandidates(physSel)
        pDMesonJetsTask.SetApplyKinematicCuts(False)

        if doRecLevel and not doResponse:
            # D0
            if config["charged_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.4)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.6)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpiLikeSign, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.4)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpiLikeSign, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.6)
        
            if config["full_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.2)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.4)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpiLikeSign, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.2)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpiLikeSign, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.4)
            # D*
            if config["charged_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.4)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kChargedJet, 0.6)
        
            if config["full_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.2)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kNoMC, ROOT.AliJetContainer.kFullJet, 0.4)
            
        if config["MC"] and doSignalOnly:
            # D0
            if config["charged_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kChargedJet, 0.4)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kChargedJet, 0.6)
        
            if config["full_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kFullJet, 0.2)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kFullJet, 0.4)
        
            # D*
            if config["charged_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kChargedJet, 0.4)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kChargedJet, 0.6)
        
            if config["full_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kFullJet, 0.2)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kSignalOnly, ROOT.AliJetContainer.kFullJet, 0.4)
                
        if config["MC"] and doMCTruth:
            # D0
            if config["charged_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kChargedJet, 0.4)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kChargedJet, 0.6)
        
            if config["full_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kFullJet, 0.2)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kD0toKpi, ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kFullJet, 0.4)
        
            # D*
            if config["charged_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kChargedJet, 0.4)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kChargedJet, 0.6)
        
            if config["full_jets"]:
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kFullJet, 0.2)
                pDMesonJetsTask.AddAnalysisEngine(ROOT.AliAnalysisTaskDmesonJets.kDstartoKpipi, ROOT.AliAnalysisTaskDmesonJets.kMCTruth, ROOT.AliJetContainer.kFullJet, 0.4)

    tasks = mgr.GetTasks()
    for task in tasks:
        if isinstance(task, ROOT.AliAnalysisTaskEmcal) or isinstance(task, ROOT.AliAnalysisTaskEmcalLight):
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
        chain = ROOT.CreateAODChain(config["file_list"], nFiles, 0, False, "AliAOD.VertexingHF.root")
    elif mode is helperFunctions.AnaMode.ESD:
        ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
        chain = ROOT.CreateESDChain(config["file_list"], nFiles, 0, False)
        
    if debugLevel == 0:
        mgr.SetUseProgressBar(1, 250)

    mgr.SetDebugLevel(debugLevel)

    #To have more debug info
    #mgr.AddClassDebug("AliAnalysisTaskDmesonJets", ROOT.AliLog.kDebug+100)

    #start analysis
    print "Starting Analysis..."
    mgr.StartAnalysis("local", chain, nEvents)

if __name__ == '__main__':
    # runJetDmeson.py executed as script

    parser = argparse.ArgumentParser(description='Jet D meson analysis.')
    parser.add_argument('config',
                        help='YAML configuration file')
    parser.add_argument('-f', '--n-files', 
                        default=100,
                        type=int,
                        help='Number of files to be analyzed')
    parser.add_argument('-e', '--n-events',
                        type=int,
                        default=2000,
                        help='Number of events to be analyzed')
    parser.add_argument('--no-signal-only',action='store_const',
                        default=False, const=True,
                        help='Signal only analysis')
    parser.add_argument('--no-rec-level',action='store_const',
                        default=False, const=True,
                        help='Reconstructed level analysis')
    parser.add_argument('--no-mc-truth', action='store_const',
                        default=False, const=True,
                        help='MC truth analysis')
    parser.add_argument('--no-incl-jets', action='store_const',
                        default=False, const=True,
                        help='No inclusive jets')
    parser.add_argument('--response', action='store_const',
                        default=False, const=True,
                        help='MC truth analysis')
    parser.add_argument('--task-name',
                        default="JetDmesonAna",
                        help='Task name')
    parser.add_argument('-d', '--debug-level', 
                        default=0,
                        type=int,
                        help='Debug level')
    args = parser.parse_args()

    main(args.config, args.n_files, args.n_events, not args.no_rec_level, not args.no_signal_only, not args.no_mc_truth, args.response, args.no_incl_jets,
         args.task_name, args.debug_level)
