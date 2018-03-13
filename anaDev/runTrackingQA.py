#!/usr/bin/env python
# python script to test Jet D meson analysis

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import helperFunctions
import yaml


def ExtractTriggerSelection(triggerList):
    r = 0
    for t in triggerList:
        if hasattr(ROOT.AliVEvent, t):
            r = r | getattr(ROOT.AliVEvent, t)
            print("Trigger '{}' added".format(t))
        else:
            print("Error: could not parse trigger name '{}'".format(t))
            exit(1)
    return r


def main(configFileName, nFiles, nEvents, debugLevel=0):
    f = open(configFileName, 'r')
    config = yaml.load(f)
    f.close()

    physSel = ExtractTriggerSelection(config["trigger"])
    print("Trigger selection is {}".format(physSel))

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

    # Analysis manager
    mgr = ROOT.AliAnalysisManager("TrackingQA")

    helperFunctions.LoadMacros()

    if mode is helperFunctions.AnaMode.AOD:
        helperFunctions.AddAODHandler()
    elif mode is helperFunctions.AnaMode.ESD:
        helperFunctions.AddESDHandler()
        if config["MC"]:
            helperFunctions.AddMCHandler()

    if mode is helperFunctions.AnaMode.ESD:
        task = helperFunctions.AddTaskCDBConnect()
        # task.SetFallBackToRaw(True)

    # Physics selection task
    if not config["MC"]: ROOT.AddTaskPhysicsSelection()

    if mode is helperFunctions.AnaMode.ESD and config["MC"]:
        ROOT.AddTaskMCTrackSelector("mcparticles", False, False, -1, False)

    if config["cent_type"] == "new":
        ROOT.AddTaskMultSelection(False)
    elif config["cent_type"] == "old":
        print("Not configured for old centrality framework!")
        return

    if config["MC"]:
        print "Running on a MC production"
    else:
        print "Running on data"

    trackQAtask = ROOT.AliEmcalTrackingQATask.AddTaskTrackingQA(config["MC"])
    trackQAtask.SelectCollisionCandidates(physSel)

    tasks = mgr.GetTasks()
    for task in tasks:
        if isinstance(task, ROOT.AliAnalysisTaskEmcalLight):
            task.SetIsPythia(config["is_pythia"])
            if config["is_pythia"]: task.SetMaxMinimumBiasPtHard(5)
            if config["run_period"] == "LHC15o": task.SetSwitchOffLHC15oFaultyBranches(True)
            if config["beam_type"] == "pp":
                task.SetForceBeamType(ROOT.AliAnalysisTaskEmcalLight.kpp)
            elif config["beam_type"] == "PbPb":
                task.SetForceBeamType(ROOT.AliAnalysisTaskEmcalLight.kAA)
            elif config["beam_type"] == "pPb":
                task.SetForceBeamType(ROOT.AliAnalysisTaskEmcalLight.kpA)
            if config["cent_type"] == "new":
                task.SetCentralityEstimation(ROOT.AliAnalysisTaskEmcalLight.kNewCentrality)
            elif config["cent_type"] == "old":
                task.SetCentralityEstimation(ROOT.AliAnalysisTaskEmcalLight.kOldCentrality)
            else:
                task.SetCentralityEstimation(ROOT.AliAnalysisTaskEmcalLight.kNoCentrality)

        if isinstance(task, ROOT.AliAnalysisTaskEmcal):
            if config["beam_type"] == "pp":
                task.SetForceBeamType(ROOT.AliAnalysisTaskEmcal.kpp)
            elif config["beam_type"] == "PbPb":
                task.SetForceBeamType(ROOT.AliAnalysisTaskEmcal.kAA)
                task.SetUseNewCentralityEstimation(True)

    res = mgr.InitAnalysis()

    if not res:
        print "Error initializing the analysis!"
        exit(1)

    mgr.PrintStatus()

    outFile = ROOT.TFile("train.root", "RECREATE")
    outFile.cd()
    mgr.Write()
    outFile.Close()

    chain = None
    if mode is helperFunctions.AnaMode.AOD:
        ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C")
        chain = ROOT.CreateAODChain(config["file_list"], nFiles, 0, False)
    elif mode is helperFunctions.AnaMode.ESD:
        ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
        chain = ROOT.CreateESDChain(config["file_list"], nFiles, 0, False)

    if debugLevel == 0:
        mgr.SetUseProgressBar(1, 250)
    else:
        # To have more debug info
        mgr.AddClassDebug("AliEmcalTrackingQATask", ROOT.AliLog.kDebug + 100)
        mgr.AddClassDebug("AliEmcalContainer", ROOT.AliLog.kDebug + 100)
        mgr.AddClassDebug("AliMCParticleContainer", ROOT.AliLog.kDebug + 100)
        mgr.AddClassDebug("AliEmcalMCTrackSelector", ROOT.AliLog.kDebug + 100)

    mgr.SetDebugLevel(debugLevel)

    # start analysis
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
    parser.add_argument('-d', '--debug-level',
                        default=0,
                        type=int,
                        help='Debug level')
    args = parser.parse_args()

    main(args.config, args.n_files, args.n_events, args.debug_level)
