#!/usr/bin/env python
#python script to test EMCal trigger QA analysis

import argparse
import ROOT
import helperFunctions

def main(fileList, nFiles, nEvents, runPeriod, strmode="AOD",
         doChargedJets=False, doFullJets=False, doNeutralJets=False, doTrackQA=False, doTriggerQA=False,
         physSel=0, taskName="EmcalTriggerQA", debugLevel=0):

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
    elif mode is AnaMode.ESD:
        helperFunctions.AddESDHandler()
        
    #Physics selection task
    if mode is helperFunctions.AnaMode.ESD:
        ROOT.AddTaskPhysicsSelection()

    #Setup task
    #OCDBpath = "local:///Volumes/DATA/ALICE/OCDB/2010";
    OCDBpath = "raw://";
    pSetupTask = ROOT.AliEmcalSetupTask("EmcalSetupTask");
    pSetupTask.SetNoOCDB(1)
    mgr.AddTask(pSetupTask)
    cinput = mgr.GetCommonInputContainer()
    mgr.ConnectInput(pSetupTask, 0,  cinput)
    pSetupTask.SetOcdbPath(OCDBpath)
    helperFunctions.PrepareEMCAL(physSel, True, True, False, False, False)

    if doTrackQA:
        pJetQATask = ROOT.AddTaskEmcalJetQA("usedefault", "usedefault", "usedefault")
    else:
        pJetQATask = ROOT.AddTaskEmcalJetQA("", "usedefault", "usedefault")
    pJetQATask.SetNeedEmcalGeom(True)
    pJetQATask.SelectCollisionCandidates(physSel)
    pJetQATask.SetHistoBins(150, 0, 150)
    
    #Trigger QA
    if doTriggerQA:
        pTriggerMakerTask = ROOT.AddTaskEmcalTriggerMakerNew()
        pTriggerMakerTask.SetUseL0Amplitudes(True)
        pTriggerQATask = ROOT.AddTaskEmcalTriggerQAPP()
        pTriggerQATask.GetTriggerQA().EnablePatchType(ROOT.AliEmcalTriggerQAPP.kOnlinePatch, False)
        pTriggerQATask.GetTriggerQA().EnablePatchType(ROOT.AliEmcalTriggerQAPP.kOfflinePatch, False)
        pTriggerQATask.SetADCperBin(5)
        pTriggerQATask.SetMinAmplitude(25)

    #Charged jet analysis
    if doChargedJets:
        pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pChJetTask.SelectCollisionCandidates(physSel)

        pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", ROOT.AliJetContainer.antikt_algorithm, 0.6, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pChJetTask.SelectCollisionCandidates(physSel)

    #Full jet analysis
    if doFullJets:
        pJetTask = ROOT.AddTaskEmcalJet("usedefault", "usedefault", ROOT.AliJetContainer.antikt_algorithm, 0.2, ROOT.AliJetContainer.kFullJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pJetTask.SelectCollisionCandidates(physSel)

        pJetTask = ROOT.AddTaskEmcalJet("usedefault", "usedefault", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kFullJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pJetTask.SelectCollisionCandidates(physSel)

    #Neutral jet analysis
    if doNeutralJets:
        pJetTask = ROOT.AddTaskEmcalJet("", "usedefault", ROOT.AliJetContainer.antikt_algorithm, 0.2, ROOT.AliJetContainer.kNeutralJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pJetTask.SelectCollisionCandidates(physSel)

        pJetTask = ROOT.AddTaskEmcalJet("", "usedefault", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kNeutralJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pJetTask.SelectCollisionCandidates(physSel)

    pSpectraTask = ROOT.AddTaskEmcalJetSpectraQA("usedefault", "usedefault")
    pSpectraTask.SelectCollisionCandidates(physSel)
    pSpectraTask.SetHistoBins(200, 0, 200)

    if doChargedJets:
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kTPCfid)
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.6, ROOT.AliJetContainer.kTPCfid)

    if doFullJets:
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.2, ROOT.AliJetContainer.kEMCALfid)
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kEMCALfid)

    if doNeutralJets:
        jetCont = pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kNeutralJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.2, ROOT.AliJetContainer.kEMCALfid)
        jetCont.SetLeadingHadronType(1)
        jetCont = pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kNeutralJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kEMCALfid)
        jetCont.SetLeadingHadronType(1)

    tasks = mgr.GetTasks()
    for task in tasks:
        if isinstance(task, ROOT.AliAnalysisTaskEmcal):
            task.SetForceBeamType(ROOT.AliAnalysisTaskEmcal.kpp)
            task.SetVzRange(-999, -999)
	
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
        chain = ROOT.CreateAODChain(fileList, nFiles, 0, False)
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
    parser.add_argument('--charged-jets', action='store_const',
                        default=False, const=True,
                        help='Charged jet analysis')
    parser.add_argument('--full-jets', action='store_const',
                        default=False, const=True,
                        help='Full jet analysis')
    parser.add_argument('--neutral-jets', action='store_const',
                        default=False, const=True,
                        help='Full jet analysis')
    parser.add_argument('--track-qa', action='store_const',
                        default=False, const=True,
                        help='Track QA')
    parser.add_argument('--trigger-qa', action='store_const',
                        default=False, const=True,
                        help='Track QA')
    parser.add_argument('--phys-sel',
                        default=0, 
                        help='Physics selection')
    parser.add_argument('--task-name',
                        default="JetDmesonAna",
                        help='Task name')
    parser.add_argument('-d', '--debug-level', 
                        default=0,
                        type=int,
                        help='Debug level')
    args = parser.parse_args()
    
    main(args.fileList, args.n_files, args.n_events, args.run_period, args.mode,
         args.charged_jets, args.full_jets, args.neutral_jets, args.track_qa, args.trigger_qa,
         args.phys_sel, args.task_name, args.debug_level)
