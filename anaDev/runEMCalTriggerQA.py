#!/usr/bin/env python
#python script to test EMCal trigger QA analysis

import argparse
import ROOT
import helperFunctions
import os
import yaml

def GenerateConfig(args):
    f = open(args.configYAML, 'r')
    config = yaml.load(f)
    f.close()

    if args.file_list:
        config["file_list"] = args.file_list

    if args.run_period:
        config["run_period"] = args.run_period

    if args.mode:
        config["mode"] = args.mode

    if args.bad_fastor:
        config["bad_fastor"] = args.bad_fastor

    if args.task_name:
        config["task_name"] = args.task_name

    config["n_files"] = args.n_files
    config["n_events"] = args.n_events
    config["charged_jets"] = args.charged_jets
    config["full_jets"] = args.full_jets
    config["neutral_jets"] = args.neutral_jets
    config["track_qa"] = args.track_qa
    config["cluster_qa"] = args.cluster_qa
    config["trigger_qa"] = args.trigger_qa
    config["debug_level"] = args.debug_level

    return config

def AddTriggerQATasks(config, trigger, physSel):

    if config["track_qa"] and config["cluster_qa"]:
        pJetQATask = ROOT.AddTaskEmcalJetQA("usedefault", "usedefault", "usedefault", trigger["label"])
        pJetQATask.SetNeedEmcalGeom(True)
    elif config["cluster_qa"]:
        pJetQATask = ROOT.AddTaskEmcalJetQA("", "usedefault", "usedefault", trigger["label"])
        pJetQATask.SetNeedEmcalGeom(True)
    elif config["track_qa"]:
        pJetQATask = ROOT.AddTaskEmcalJetQA("", "usedefault", "usedefault", trigger["label"])

    if config["track_qa"] or config["cluster_qa"]:
        if trigger.has_key("accept"):
            for acc in trigger["accept"]:
                pJetQATask.AddAcceptedTriggerClass(acc)
        if trigger.has_key("reject"):
            for rej in trigger["reject"]:
                pJetQATask.AddRejectedTriggerClass(rej)
        pJetQATask.SelectCollisionCandidates(physSel)
        pJetQATask.SetPtBin(1, 150)
    
    #Trigger QA
    if config["trigger_qa"]:
        if config["run_period"] == "LHC16q":
            pTriggerQATask = ROOT.AddTaskEmcalTriggerQA("EmcalTriggers", "", "", 0, False, trigger["label"])
            pTriggerQATask.EnableDCal(True)
        elif config["run_period"] == "LHC16d" or config["run_period"] == "LHC16e" or config["run_period"] == "LHC16f":
            pTriggerQATask = ROOT.AddTaskEmcalTriggerQA("EmcalTriggers", "", "", 0, False, trigger["label"])
            pTriggerQATask.EnableDCal(True)
        elif config["run_period"] == "LHC16c":
            pTriggerQATask = ROOT.AddTaskEmcalTriggerQA("EmcalTriggers", "", "", 0, False, trigger["label"])
            pTriggerQATask.EnableDCal(True)
        elif config["run_period"] == "LHC16b":
            pTriggerQATask = ROOT.AddTaskEmcalTriggerQA("EmcalTriggers", "", "", 0, False, trigger["label"])
            pTriggerQATask.EnableDCal(True)
        elif config["run_period"] == "LHC15o":
            pTriggerQATask = ROOT.AddTaskEmcalTriggerQA("EmcalTriggers", "", "", 5, True, trigger["label"])
            pTriggerQATask.SetUseNewCentralityEstimation(True)
            pTriggerQATask.EnableDCal(True)
        elif config["run_period"] == "LHC15j":
            pTriggerQATask = ROOT.AddTaskEmcalTriggerQA("EmcalTriggers", "", "", 0, False, trigger["label"])
            pTriggerQATask.EnableDCal(True)
        elif config["run_period"].startswith("LHC13"):
            pTriggerQATask = ROOT.AddTaskEmcalTriggerQA("EmcalTriggers", "", "", 0, False, trigger["label"])
            pTriggerQATask.EnableDCal(False)
        elif config["run_period"].startswith("LHC12"):
            pTriggerQATask = ROOT.AddTaskEmcalTriggerQA("EmcalTriggers", "", "", 0, False, trigger["label"])
            pTriggerQATask.EnableDCal(False)
        elif config["run_period"] == "LHC11h":
            pTriggerQATask = ROOT.AddTaskEmcalTriggerQA("EmcalTriggers", "", "", 4, True, trigger["label"])
            pTriggerQATask.EnableDCal(False)
        elif config["run_period"].startswith("LHC11") and config["run_period"] != "LHC11h":
            pTriggerQATask = ROOT.AddTaskEmcalTriggerQA("EmcalTriggers", "", "", 0, False, trigger["label"])
            pTriggerQATask.EnableDCal(False)
            
        if trigger.has_key("accept"):
            for acc in trigger["accept"]:
                pTriggerQATask.AddAcceptedTriggerClass(acc)
        if trigger.has_key("reject"):
            for rej in trigger["reject"]:
                pTriggerQATask.AddRejectedTriggerClass(rej)
            
        if config.has_key("bad_fastor"):
            pTriggerQATask.GetTriggerQA().ReadFastORBadChannelFromFile(config["bad_fastor"])
        #pTriggerQATask.EnableHistogramsByTimeStamp(120)
        pTriggerQATask.SetADCperBin(8)
        pTriggerQATask.SetMinAmplitude(0)
        pTriggerQATask.SelectCollisionCandidates(physSel)

    if config["charged_jets"] or config["full_jets"] or config["neutral_jets"]:
        pSpectraTask = ROOT.AddTaskEmcalJetSpectraQA("usedefault", "usedefault", 0.15, 0.30, trigger["label"])
        pSpectraTask.SelectCollisionCandidates(physSel)
        pSpectraTask.SetPtBin(1, 250)
        if trigger.has_key("accept"):
            for acc in trigger["accept"]:
                pSpectraTask.AddAcceptedTriggerClass(acc)
        if trigger.has_key("reject"):
            for rej in trigger["reject"]:
                pSpectraTask.AddRejectedTriggerClass(rej)

    if config["charged_jets"]:
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kTPCfid)
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kChargedJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.6, ROOT.AliJetContainer.kTPCfid)

    if config["full_jets"]:
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.2, ROOT.AliJetContainer.kEMCALfid)
        pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kFullJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kEMCALfid)

    if config["neutral_jets"]:
        jetCont = pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kNeutralJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.2, ROOT.AliJetContainer.kEMCALfid)
        jetCont.SetLeadingHadronType(1)
        jetCont = pSpectraTask.AddJetContainer(ROOT.AliJetContainer.kNeutralJet, ROOT.AliJetContainer.antikt_algorithm, ROOT.AliJetContainer.pt_scheme, 0.4, ROOT.AliJetContainer.kEMCALfid)
        jetCont.SetLeadingHadronType(1)

def main(config):
    ROOT.gSystem.Load("libCGAL")

    ROOT.gInterpreter.ProcessLine(os.path.expandvars('.L $ALICE_ROOT/include/AliEMCALTriggerConstants.h'))

    physSel = 0

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
    print("Setting local analysis for {0} files from list {1}, max events = {2}".format(config["n_files"], config["file_list"], config["n_events"]))

    #AliVEvent::kINT7, AliVEvent::kMB, AliVEvent::kCentral, AliVEvent::kSemiCentral
    #AliVEvent::kEMCEGA, AliVEvent::kEMCEJ

    #Analysis manager
    mgr = ROOT.AliAnalysisManager(config["task_name"])

    helperFunctions.LoadMacros()

    if mode is helperFunctions.AnaMode.AOD:
        helperFunctions.AddAODHandler()
    elif mode is helperFunctions.AnaMode.ESD:
        helperFunctions.AddESDHandler()

    #CDB connect
    #pCDBConnect = ROOT.AddTaskCDBconnect()
    pCDBConnect = ROOT.AliTaskCDBconnect("CDBconnect", "cvmfs://", 0)
    mgr.AddTask(pCDBConnect)
    cinput1 = mgr.GetCommonInputContainer()    
    mgr.ConnectInput(pCDBConnect,  0, cinput1)
    pCDBConnect.SetFallBackToRaw(True)

    #Physics selection task
    if mode is helperFunctions.AnaMode.ESD and physSel:
        ROOT.AddTaskPhysicsSelection()
    
    if config["run_period"] == "LHC15o":
        ROOT.AddTaskMultSelection(False)

    # EMCal prep
    if config["cluster_qa"] or config["trigger_qa"]:
        helperFunctions.PrepareEMCAL(physSel, True, True, config["cluster_qa"] or config["full_jets"] or config["neutral_jets"], config["full_jets"] or config["neutral_jets"], config["full_jets"] or config["neutral_jets"])

    #Trigger QA
    if config["trigger_qa"]:
        pTriggerMakerTask = ROOT.AddTaskEmcalTriggerMakerNew("EmcalTriggers")
        pTriggerMakerTask.SelectCollisionCandidates(physSel)
        pTriggerMakerTask.GetTriggerMaker().SetFastORandCellThresholds(0, 0, 0)
        
        if config.has_key("bad_fastor"):
            pTriggerMakerTask.GetTriggerMaker().ReadFastORBadChannelFromFile(config["bad_fastor"])

        if config["run_period"] == "LHC16q":
            pTriggerMakerTask.GetTriggerMaker().ConfigureForPP2015()
        elif config["run_period"] == "LHC16d" or config["run_period"] == "LHC16e" or config["run_period"] == "LHC16f":
            pTriggerMakerTask.GetTriggerMaker().ConfigureForPP2015()
        elif config["run_period"] == "LHC16c":
            pTriggerMakerTask.GetTriggerMaker().ConfigureForPP2015()
        elif config["run_period"] == "LHC16b":
            pTriggerMakerTask.GetTriggerMaker().ConfigureForPP2015()
        elif config["run_period"] == "LHC15o":
            pTriggerMakerTask.GetTriggerMaker().ConfigureForPbPb2015()
        elif config["run_period"] == "LHC15j":
            pTriggerMakerTask.GetTriggerMaker().ConfigureForPP2015()
            pTriggerMakerTask.SetUseL0Amplitudes(True)
        elif config["run_period"].startswith("LHC13"):
            pTriggerMakerTask.GetTriggerMaker().ConfigureForPPb2013()
        elif config["run_period"].startswith("LHC12"):
            pTriggerMakerTask.GetTriggerMaker().ConfigureForPP2012()
        elif config["run_period"] == "LHC11h":
            pTriggerMakerTask.GetTriggerMaker().ConfigureForPbPb2011()
        elif config["run_period"].startswith("LHC11") and config["run_period"] != "LHC11h":
            pTriggerMakerTask.GetTriggerMaker().ConfigureForPP2011()

    #Charged jet analysis
    if config["charged_jets"]:
        pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pChJetTask.SelectCollisionCandidates(physSel)

        pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", ROOT.AliJetContainer.antikt_algorithm, 0.6, ROOT.AliJetContainer.kChargedJet, 0.15, 0., 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pChJetTask.SelectCollisionCandidates(physSel)

    #Full jet analysis
    if config["full_jets"]:
        pJetTask = ROOT.AddTaskEmcalJet("usedefault", "usedefault", ROOT.AliJetContainer.antikt_algorithm, 0.2, ROOT.AliJetContainer.kFullJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pJetTask.SelectCollisionCandidates(physSel)

        pJetTask = ROOT.AddTaskEmcalJet("usedefault", "usedefault", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kFullJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pJetTask.SelectCollisionCandidates(physSel)

    #Neutral jet analysis
    if config["neutral_jets"]:
        pJetTask = ROOT.AddTaskEmcalJet("", "usedefault", ROOT.AliJetContainer.antikt_algorithm, 0.2, ROOT.AliJetContainer.kNeutralJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pJetTask.SelectCollisionCandidates(physSel)

        pJetTask = ROOT.AddTaskEmcalJet("", "usedefault", ROOT.AliJetContainer.antikt_algorithm, 0.4, ROOT.AliJetContainer.kNeutralJet, 0.15, 0.30, 0.1, ROOT.AliJetContainer.pt_scheme, "Jet", 0., False, False)
        pJetTask.SelectCollisionCandidates(physSel)

    for trigger in config["triggers"]:
        AddTriggerQATasks(config, trigger, physSel)

    if config["run_period"] == "LHC11h" or config["run_period"].startswith("LHC13") or config["run_period"] == "LHC15o":
        isPP = False
    else:
        isPP = True

    tasks = mgr.GetTasks()
    for task in tasks:
        if isinstance(task, ROOT.AliAnalysisTaskEmcal) or isinstance(task, ROOT.AliAnalysisTaskEmcalLight):
            if isPP:
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
        chain = ROOT.CreateAODChain(config["file_list"], config["n_files"], 0, False)
    elif mode is helperFunctions.AnaMode.ESD:
        ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
        chain = ROOT.CreateESDChain(config["file_list"], config["n_files"], 0, False)

    if config["debug_level"] == 0:
        mgr.SetUseProgressBar(1, 250)
    else:
        mgr.SetUseProgressBar(0, 0)

    mgr.SetDebugLevel(config["debug_level"])

    #To have more debug info
    #pMgr->AddClassDebug("AliEmcalClusTrackMatcherTask", AliLog::kDebug+100);

    #start analysis
    print "Starting Analysis..."
    mgr.StartAnalysis("local", chain, config["n_events"])

if __name__ == '__main__':
    # runDJetCorrAnalysis.py executed as script
    
    parser = argparse.ArgumentParser(description='EMCal trigger QA analysis.')
    parser.add_argument('configYAML', metavar='configYAML',
                        help='YAML configuation file')
    
    #these are defined in the YAML file but can be overridden
    parser.add_argument('--mode',
                        help='Analysis mode (ESD or AOD)')
    parser.add_argument('--file-list',
                        help='File list to be analyzed')
    parser.add_argument('-r', '--run-period',
                        help='Run period (e.g. LHC10b)')
    parser.add_argument('--task-name',
                        help='Task name')
    parser.add_argument('--bad-fastor',
                        help='Name of the file containing a bad FastOR list')
    
    #these are not defined in the YAML file
    parser.add_argument('-f', '--n-files', 
                        default=100,
                        type=int,
                        help='Number of files to be analyzed')
    parser.add_argument('-e', '--n-events',
                        default=1234567890,
                        type=int,
                        help='Number of events to be analyzed')
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
    parser.add_argument('--cluster-qa', action='store_const',
                        default=False, const=True,
                        help='Cluster QA')
    parser.add_argument('--trigger-qa', action='store_const',
                        default=False, const=True,
                        help='Track QA')
    parser.add_argument('-d', '--debug-level', 
                        default=0,
                        type=int,
                        help='Debug level')

    args = parser.parse_args()

    config = GenerateConfig(args)

    main(config)
