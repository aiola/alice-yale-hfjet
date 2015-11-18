// runJetDmeson.C

class AliESDInputHandler;
class AliAODInputHandler;
class AliAODOutputHandler;
class AliAnalysisManager;

void LoadLibs();
void LoadMacros();

//______________________________________________________________________________
void runJetDmeson(const char   *cDataType     = "AOD",                                    // set the analysis type, AOD or ESD
		  const char   *cLocalFiles   = "fileLists/files_LHC10c_pass4_AOD.txt",   // set the local list file
		  UInt_t        iNumFiles     = 100,                                      // number of files analyzed locally
		  UInt_t        iNumEvents    = 50000,                                    // number of events to be analyzed
		  const char   *cRunPeriod    = "LHC11h",                                 // set the run period
		  const char   *cTaskName     = "DmesonJetAna"                            // sets name of analysis manager
		  )
{

  const Bool_t   bDoChargedJets          = kTRUE;
  const Bool_t   bDoHF                   = kTRUE;
  const Bool_t   bIsPbPb                 = kFALSE;
  const Double_t kJetRadius              = 0.6;
  const Double_t kTrackPtCut             = 0.15;
  const Double_t kJetAreaCut             = 0.;
  const Double_t kJetPtCut               = 1.;
  const Double_t kGhostArea              = 0.01;
  const Int_t    kHistoType              = 1;
  const Int_t    eFlavourJetMatchingType = AliAnalysisTaskDmesonJetCorrelations::kCandidateConstituentMatching;
  // kGeometricalMatching, kDaughterConstituentMatching, kCandidateConstituentMatching, kJetLoop

  TString sTracksName("tracks");
  TString sChJetsName;
  
  enum eDataType { kAod, kEsd };

  eDataType iDataType;
  if (!strcmp(cDataType, "ESD")) {
    iDataType = kEsd;
  }
  else if (!strcmp(cDataType, "AOD")) {
    iDataType = kAod;
  }
  else {
    Printf("Incorrect data type option, check first argument of run macro.");
    Printf("datatype = AOD or ESD");
    return;
  }

  Printf("%s analysis chosen.", cDataType);

  TString sLocalFiles(cLocalFiles);
  if (sLocalFiles == "") {
    Printf("You need to provide the list of local files!");
    return;
  }
  Printf("Setting local analysis for %d files from list %s, max events = %d", iNumFiles, sLocalFiles.Data(), iNumEvents);

  LoadLibs();
  LoadMacros();

  // AliEmcalPhysicsSelection::kEmcalOk, AliEmcalPhysicsSelection::kEmcalH,
  // AliVEvent::kINT7, AliVEvent::kMB, AliVEvent::kCentral, AliVEvent::kSemiCentral,
  // AliVEvent::kEMCEGA, AliVEvent::kEMCEJE
  UInt_t kPrePhysSel = AliVEvent::kMB;
  UInt_t kPhysSel = AliVEvent::kMB;

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);

  if (iDataType == kAod) {
    AliAODInputHandler* pAODHandler = AddAODHandler();
  }
  else {   //ESD or skimmed ESD
    AliESDInputHandler* pESDHandler = AddESDHandler();
  }

  if (0) {
    AliAODHandler* pAODOutHandler = AddAODOutputHandler();
  }

  // Physics selection task
  if (1) {
    AliEmcalPhysicsSelectionTask *pPhysSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE,
                                                                              kPrePhysSel,
                                                                              5, 5, 10, kTRUE);
  }

  // Centrality task
  if (iDataType == kEsd && bIsPbPb) {
    AliCentralitySelectionTask *pCentralityTask = AddTaskCentrality(kTRUE);
    pCentralityTask->SelectCollisionCandidates(kPhysSel);
  }
    
  // PID response
  if (bDoHF) {
    AliAnalysisTaskPIDResponse* pPIDtask = (AliAnalysisTaskPIDResponse*)AddTaskPIDResponse(kFALSE);
    pPIDtask->SelectCollisionCandidates(kPhysSel);
  }
  
  // Charged jet analysis
  if (bDoChargedJets) {
    AliEmcalJetTask *pChJetTask = AddTaskEmcalJet(sTracksName, "", 1, kJetRadius, 1, kTrackPtCut, 0., kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, kFALSE);
    pChJetTask->SelectCollisionCandidates(kPhysSel);
    pChJetTask->SetFilterHybridTracks(kTRUE);
    sChJetsName = pChJetTask->GetName();

    AliAnalysisTaskSAJF *pSpectraChTask = AddTaskSAJF(sTracksName, "", sChJetsName, "",  kJetRadius, kJetPtCut, kJetAreaCut, "TPC");
    pSpectraChTask->SetNLeadingJets(1);
    pSpectraChTask->SelectCollisionCandidates(kPhysSel);
    pSpectraChTask->SetHistoType(kHistoType);
  }


  if (bDoHF) {
    AliAnalysisTaskSEDmesonsFilterCJ* pDStarMesonFilterTask = AddTaskSEDmesonsFilterCJ(AliAnalysisTaskSEDmesonsFilterCJ::kDstartoKpipi,
										       "DStartoKpipiCuts.root",
										       kFALSE,  //   Bool_t theMCon
										       kTRUE,   //   Bool_t reco
										       "rec");
    pDStarMesonFilterTask->SelectCollisionCandidates(kPhysSel);
    pDStarMesonFilterTask->SetCombineDmesons(kTRUE);
    AliParticleContainer* trackContDStar = pDStarMesonFilterTask->AddParticleContainer(sTracksName);
    trackContDStar->SetClassName("AliAODTrack");
    trackContDStar->SetFilterHybridTracks(kTRUE);

    TString sTracksDStarName = "DcandidatesAndTracksDStarrec";
    
    AliEmcalJetTask *pChJetDStarTask = AddTaskEmcalJet(sTracksDStarName, "", 1, kJetRadius, 1, kTrackPtCut, 0., kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, 1);
    pChJetDStarTask->SelectCollisionCandidates(kPhysSel);
    TString sChJetsDStarName = pChJetDStarTask->GetName();

    if (1) {
      AliAnalysisTaskSAJF *pSpectraChDStarTask = AddTaskSAJF(sTracksDStarName, "", sChJetsDStarName, "",  kJetRadius, kJetPtCut, kJetAreaCut, "TPC");
      pSpectraChDStarTask->SetNLeadingJets(1);
      pSpectraChDStarTask->SelectCollisionCandidates(kPhysSel);
      pSpectraChDStarTask->SetHistoType(kHistoType);
    }

    AliAnalysisTaskDmesonJetCorrelations* pDStarMesonJetCorr = AddTaskDmesonJetCorr(AliAnalysisTaskDmesonJetCorrelations::kDstartoKpipi, "", 
										    sTracksDStarName, "", sChJetsDStarName, "",
										    kJetRadius, kJetPtCut, kJetAreaCut, "TPC", 0, "Dcandidates", kFALSE,
										    "AliAnalysisTaskDmesonJetCorrelations", "rec");
    pDStarMesonJetCorr->SetMaxR(kJetRadius);
    pDStarMesonJetCorr->SetMatchingType(eFlavourJetMatchingType);
    pDStarMesonJetCorr->SetPlotOnlyAcceptedJets(kTRUE);
    pDStarMesonJetCorr->SetShowDeltaEta(kTRUE);
    pDStarMesonJetCorr->SetShowDeltaPhi(kTRUE);
    pDStarMesonJetCorr->SetShow2ProngInvMass(kTRUE);
    pDStarMesonJetCorr->SetShowInvMass(kTRUE);
    pDStarMesonJetCorr->SetShowJetConstituents(kTRUE);
    pDStarMesonJetCorr->SelectCollisionCandidates(kPhysSel);

    // QA task
    AliAnalysisTaskSAQA *pQADStarTask = AddTaskSAQA(sTracksDStarName, "", "", "", "", 0.2, 1, 0, kTrackPtCut, 0., "TPC");
    pQADStarTask->SelectCollisionCandidates(kPhysSel);
    pQADStarTask->SetHistoBins(200, 0, 30);
  }

  if (1) {
    AliAnalysisTaskSEDmesonsFilterCJ* pD0mesonFilterTask = AddTaskSEDmesonsFilterCJ(AliAnalysisTaskSEDmesonsFilterCJ::kD0toKpi,
										    "DStartoKpipiCuts.root",
										    kFALSE,  //   Bool_t theMCon
										    kTRUE,   //   Bool_t reco
										    "rec");
    pD0mesonFilterTask->SetCombineDmesons(kTRUE);
    pD0mesonFilterTask->SelectCollisionCandidates(kPhysSel);
    AliParticleContainer* trackContD0 = pD0mesonFilterTask->AddParticleContainer(sTracksName);
    trackContD0->SetClassName("AliAODTrack");
    trackContD0->SetFilterHybridTracks(kTRUE);

    
    TString sTracksD0Name = "DcandidatesAndTracksD0rec";
    
    AliEmcalJetTask *pChJetD0Task = AddTaskEmcalJet(sTracksD0Name, "", 1, kJetRadius, 0, kTrackPtCut, 0., kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, 1);
    pChJetD0Task->SelectCollisionCandidates(kPhysSel);
    TString sChJetsD0Name = pChJetD0Task->GetName();

    AliAnalysisTaskSAJF *pSpectraChD0Task = AddTaskSAJF(sTracksD0Name, "", sChJetsD0Name, "",  kJetRadius, kJetPtCut, kJetAreaCut, "TPC");
    pSpectraChD0Task->SetNLeadingJets(1);
    pSpectraChD0Task->SelectCollisionCandidates(kPhysSel);
    pSpectraChD0Task->SetHistoType(kHistoType);
    
    AliAnalysisTaskDmesonJetCorrelations* pD0MesonJetCorr = AddTaskDmesonJetCorr(AliAnalysisTaskDmesonJetCorrelations::kD0toKpi, "", 
										 sTracksD0Name, "", sChJetsD0Name, "",
										 kJetRadius, kJetPtCut, kJetAreaCut, "TPC", 0, "Dcandidates", kFALSE,
										 "AliAnalysisTaskDmesonJetCorrelations", "rec");
    pD0MesonJetCorr->SetMaxR(kJetRadius);
    pD0MesonJetCorr->SetMatchingType(eFlavourJetMatchingType);
    pD0MesonJetCorr->SetPlotOnlyAcceptedJets(kTRUE);
    pD0MesonJetCorr->SetShowDeltaEta(kTRUE);
    pD0MesonJetCorr->SetShowDeltaPhi(kTRUE);
    pD0MesonJetCorr->SetShowJetConstituents(kTRUE);
    pD0MesonJetCorr->SelectCollisionCandidates(kPhysSel);

    AliAnalysisTaskSAQA *pQAD0Task = AddTaskSAQA(sTracksD0Name, "", "", "", "", 0.2, 1, 0, kTrackPtCut, 0., "TPC");
    pQAD0Task->SelectCollisionCandidates(kPhysSel);
    pQAD0Task->SetHistoBins(200, 0, 30);
  }

  if (1) {
    TObjArray *pTopTasks = pMgr->GetTasks();
    for (Int_t i = 0; i < pTopTasks->GetEntries(); ++i) {
      AliAnalysisTaskSE *pTask = dynamic_cast<AliAnalysisTaskSE*>(pTopTasks->At(i));
      if (!pTask) continue;
      if (pTask->InheritsFrom("AliAnalysisTaskEmcal")) {
        AliAnalysisTaskEmcal *pTaskEmcal = static_cast<AliAnalysisTaskEmcal*>(pTask);
        if (!bIsPbPb) {
          Printf("Setting beam type for task %s", pTaskEmcal->GetName());
          pTaskEmcal->SetForceBeamType(0);
        }
      }
    }
  }
	
  if (!pMgr->InitAnalysis()) return;
  pMgr->PrintStatus();

  TChain* pChain = 0;
  if (iDataType == kAod) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
    if (bDoHF) {
      pChain = CreateAODChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE, "AliAOD.VertexingHF.root");
    }
    else {
      pChain = CreateAODChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE, "");
    }
  }
  else { 
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
    pChain = CreateESDChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
  }
    
  // start analysis
  Printf("Starting Analysis...");
  pMgr->SetUseProgressBar(1, 250);
  //pMgr->SetDebugLevel(2);

  // To have more debug info
  //pMgr->AddClassDebug("AliEmcalClusTrackMatcherTask", AliLog::kDebug+100);
  //pMgr->AddClassDebug("AliAnalysisTaskSEDmesonsFilterCJ", AliLog::kDebug+100);
  //pMgr->AddClassDebug("AliAnalysisTaskDmesonJetCorrelations", AliLog::kDebug+100);
  //pMgr->AddClassDebug("AliEmcalJetTask", AliLog::kDebug+100);

  TFile *pOutFile = new TFile("train.root","RECREATE");
  pOutFile->cd();
  pMgr->Write();
  pOutFile->Close();
  delete pOutFile;

  pMgr->StartAnalysis("local", pChain, iNumEvents);
}

//______________________________________________________________________________
void LoadLibs()
{
  // load fastjet libraries 3.x
  gSystem->Load("libCGAL");
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical");
  gSystem->Load("$FASTJET/lib/libfastjetplugins");
  gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");
}

void LoadMacros()
{
  // Aliroot macros
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODOutputHandler.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskTrackingQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAJF.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskSEDmesonsFilterCJ.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskDmesonJetCorr.C");
}
