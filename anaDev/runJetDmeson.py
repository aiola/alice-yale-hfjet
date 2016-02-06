#!/usr/bin/env python
#python script to test Jet D meson analysis

import argparse
import ROOT
import helperFunctions

def main(fileList, nFiles, nEvents, runPeriod, strmode="AOD", jetRadius=0.6, doHF=True, doChargedJets=True, physSel=ROOT.AliVEvent.kMB, taskName="JetDmesonAna", debugLevel=0):

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
  
    #Charged jet analysis
    if doChargedJets:
        pChJetTask = ROOT.AddTaskEmcalJet("usedefault", "", 1, jetRadius, 1, 0.15, 0., 0.1, 1, "Jet", 0., False, False, False)
        pChJetTask.SelectCollisionCandidates(physSel)
        pChJetTask.GetParticleContainer(0).SetFilterHybridTracks(True)

        sChJetsName = pChJetTask.GetName()

        pSpectraChTask = ROOT.AddTaskEmcalJetSpectraQA("usedefault", "", sChJetsName, "",  jetRadius, 1., 0., "TPC")
        pSpectraChTask.SetNLeadingJets(1)
        pSpectraChTask.SelectCollisionCandidates(physSel)
        pSpectraChTask.SetHistoBins(200, 0, 200)

    if doHF:
        # D*
        
        pDStarMesonFilterTask = ROOT.AddTaskSEDmesonsFilterCJ(ROOT.AliAnalysisTaskSEDmesonsFilterCJ.kDstartoKpipi,
                                                              "DStartoKpipiCuts.root", False, True, "rec")
        pDStarMesonFilterTask.SelectCollisionCandidates(physSel)
        pDStarMesonFilterTask.SetCombineDmesons(True)
        trackContDStar = pDStarMesonFilterTask.AddParticleContainer("tracks")
        trackContDStar.SetFilterHybridTracks(False)

        sTracksDStarName = "DcandidatesAndTracksDStarrec"
    
        pChJetDStarTask = ROOT.AddTaskEmcalJet(sTracksDStarName, "", 1, jetRadius, 1, 0.15, 0., 0.1, 1, "Jet", 0., False, False, 0)
        pChJetDStarTask.SelectCollisionCandidates(physSel)

        sChJetsDStarName = pChJetDStarTask.GetName()

        pDStarMesonJetCorr = ROOT.AddTaskDmesonJetCorr(ROOT.AliAnalysisTaskDmesonJetCorrelations.kDstartoKpipi, "",
                                                       sTracksDStarName, "", sChJetsDStarName, "",
                                                       jetRadius, 1, 0., "TPC", 0, "Dcandidates", False,
                                                       "AliAnalysisTaskDmesonJetCorrelations", "rec")
        pDStarMesonJetCorr.SetMaxR(jetRadius)
        pDStarMesonJetCorr.SetMatchingType(ROOT.AliAnalysisTaskDmesonJetCorrelations.kCandidateConstituentMatching)
        pDStarMesonJetCorr.SetPlotOnlyAcceptedJets(True)
        pDStarMesonJetCorr.SetShowDeltaEta(True)
        pDStarMesonJetCorr.SetShowDeltaPhi(True)
        pDStarMesonJetCorr.SetShow2ProngInvMass(True)
        pDStarMesonJetCorr.SetShowInvMass(True)
        pDStarMesonJetCorr.SetShowJetConstituents(True)
        pDStarMesonJetCorr.SelectCollisionCandidates(physSel)

        # D0

        pD0mesonFilterTask = ROOT.AddTaskSEDmesonsFilterCJ(ROOT.AliAnalysisTaskSEDmesonsFilterCJ.kD0toKpi,
                                                           "D0toKpiCuts.root", False, True, "rec")
        pD0mesonFilterTask.SelectCollisionCandidates(physSel)
        pD0mesonFilterTask.SetCombineDmesons(True)
        trackContD0 = pD0mesonFilterTask.AddParticleContainer("tracks")
        trackContD0.SetFilterHybridTracks(False)
    
        sTracksD0Name = "DcandidatesAndTracksD0rec";
    
        pChJetD0Task = ROOT.AddTaskEmcalJet(sTracksD0Name, "", 1, jetRadius, 1, 0.15, 0., 0.1, 1, "Jet", 0., False, False, 0);
        pChJetD0Task.SelectCollisionCandidates(physSel)
        
        sChJetsD0Name = pChJetD0Task.GetName()

        pD0MesonJetCorr = ROOT.AddTaskDmesonJetCorr(ROOT.AliAnalysisTaskDmesonJetCorrelations.kD0toKpi, "",
                                                    sTracksD0Name, "", sChJetsD0Name, "",
                                                    jetRadius, 1, 0., "TPC", 0, "Dcandidates", False,
                                                    "AliAnalysisTaskDmesonJetCorrelations", "rec")
        pD0MesonJetCorr.SetMaxR(jetRadius)
        pD0MesonJetCorr.SetMatchingType(ROOT.AliAnalysisTaskDmesonJetCorrelations.kCandidateConstituentMatching);
        pD0MesonJetCorr.SetPlotOnlyAcceptedJets(True)
        pD0MesonJetCorr.SetShowDeltaEta(True)
        pD0MesonJetCorr.SetShowDeltaPhi(True)
        pD0MesonJetCorr.SetShowJetConstituents(True)
        pD0MesonJetCorr.SelectCollisionCandidates(physSel)
        
    tasks = mgr.GetTasks()
    for task in tasks:
        if isinstance(task, ROOT.AliAnalysisTaskEmcal):
            task.SetForceBeamType(ROOT.AliAnalysisTaskEmcal.kpp)
	
    res = mgr.InitAnalysis()

    if not res:
        print "Error initializing the analysis!"
        exit(1)
    
    mgr.PrintStatus()

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

    outFile = ROOT.TFile("train.root","RECREATE")
    outFile.cd()
    mgr.Write()
    outFile.Close()

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
    parser.add_argument('--phys-sel',
                        default=ROOT.AliVEvent.kMB, 
                        help='Physics selection')
    parser.add_argument('--task-name',
                        default="JetDmesonAna",
                        help='Task name')
    args = parser.parse_args()
    
    main(args.fileList, args.n_files, args.n_events, args.run_period, args.mode, args.jet_radius, not args.no_hf, not args.no_charged_jets, not args.phys_sel, args.task_name)
