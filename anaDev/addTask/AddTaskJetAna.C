// AddTaskJetAna.C

void LoadMacros();

void AddTaskJetAna(const char *cDataType = "AOD", const char *cRunType = "local",
                   const UInt_t kPhysSel = AliEmcalPhysicsSelection::kEmcalOk,
                   Bool_t bDoHF = kTRUE, Bool_t bDoChargedJets = kTRUE, Bool_t bDoFullJets = kFALSE,
                   Bool_t bDoTender = kFALSE, Bool_t bDoReclusterize = kFALSE, Bool_t bDoHadCorr = kFALSE, Bool_t bDoTrackingQA = kFALSE)
{
  enum eDataType { kAod, kEsd };
  enum eRunType  { kLocal, kGrid, kPdsf, kCarver };

  eRunType iRunType;
  if (!strcmp(cRunType, "grid")) 
    iRunType = kGrid;
  else if (!strcmp(cRunType, "local")) 
    iRunType = kLocal;
  else if (!strcmp(cRunType, "pdsf")) 
    iRunType = kPdsf;
  else if (!strcmp(cRunType, "carver")) 
    iRunType = kCarver;
  else {
    Printf("Incorrect run option, check first argument of run macro.");
    Printf("runtype = local, grid, pdsf, carver.");
    return;
  }

  eDataType iDataType;
  if (!strcmp(cDataType, "ESD"))
    iDataType = kEsd;
  else if (!strcmp(cDataType, "AOD"))
    iDataType = kAod;
  else {
    Printf("Incorrect data type option, check third argument of run macro.");
    Printf("cDataType = AOD or ESD.");
    return;
  }

  const Double_t kPropDist            = 440.;
  const Double_t kMinJetPt            = 0.;
  const Double_t kJetRadius           = 0.2;
  const Double_t kClusPtCut           = 0.30;
  const Double_t kTrackPtCut          = 0.15;
  const Double_t kPartLevPtCut        = 0.;
  const Double_t kJetAreaCut          = 0.557;
  const Double_t kJetPtCut            = 1.;
  const Double_t kGhostArea           = 0.005;
  const Double_t kHadCorrF            = 2.;
  const Int_t    kHistoType           = 1;
  const UInt_t   kClusterizerType     = AliEMCALRecParam::kClusterizerv2;
  const Double_t kJetLeadingTrackBias = 0.;
  const Double_t kEmcalMinEta         = -0.7;
  const Double_t kEmcalMaxEta         = 0.7;
  const Double_t kEmcalMinPhi         = 80*TMath::DegToRad();
  const Double_t kEmcalMaxPhi         = 180*TMath::DegToRad();
  const Bool_t   bForcePP             = kTRUE;
  const Double_t kEMCtimeMin          = -50e-6;
  const Double_t kEMCtimeMax          = 100e-6;
  const Double_t kEMCtimeCut          =  75e-6;

  TString sTracksName("HybridTracks");
  TString sClusName("EmcCaloClusters");

  TString sCellName;
  TString sOrigClusName;
  TString sCorrClusName;

  TString sChJetsName;
  TString sFuJetsName;

  if (iDataType == kAod) {
    sCellName = "emcalCells";
    sOrigClusName = "caloClusters";
  }
  else {
    sCellName = "EMCALCells";
    sOrigClusName = "CaloClusters";
  }

  LoadMacros();

  AliAnalysisManager *pMgr = AliAnalysisManager::GetAnalysisManager();
  if (!pMgr) {
    ::Error("AddTaskJetAna", "No analysis manager to connect to.");
    return;
  }

  if (iDataType == kEsd) {
    // Hybrid tracks maker for ESD
    AliEmcalEsdTrackFilterTask *pHybTask = AddTaskEmcalEsdTrackFilter(sTracksName, "Hybrid_LHC10b");
    pHybTask->SetDoPropagation(bDoHadCorr);
    pHybTask->SetDist(kPropDist);
    pHybTask->SelectCollisionCandidates(kPhysSel);
  }
  else if (iDataType == kAod) {
    // Hybrid tracks maker for AOD
    AliEmcalAodTrackFilterTask *pHybTask = AddTaskEmcalAodTrackFilter(sTracksName, "tracks", "LHC10b");
    pHybTask->SelectCollisionCandidates(kPhysSel);
    pHybTask->SetAttemptProp(bDoHadCorr);
    pHybTask->SetIncludeNoITS(kTRUE);
  }

  if (bDoTrackingQA) {
    AliEmcalTrackingQATask* pTrackingQAtask = AddTaskTrackingQA(sTracksName, "", kFALSE);
    pTrackingQAtask->SelectCollisionCandidates(kPhysSel);
  }

  if (1) {
    if (bDoTender || bDoReclusterize) {
      // QA task
      AliAnalysisTaskSAQA *pQATaskBefore = AddTaskSAQA("", sOrigClusName, sCellName, "", "",
						      0, 0, 0, 0.15, 0.15, "TPC", "AliAnalysisTaskSAQA_BeforeTender");
      pQATaskBefore->SetHistoBins(150, 0, 150);
      pQATaskBefore->SelectCollisionCandidates(kPhysSel);
    }
  }
  
  // Tender Supplies
  if (bDoTender) {
    Bool_t  bDistBC         = kFALSE; //switch for recalculation cluster position from bad channel 
    Bool_t  bRecalibClus    = kFALSE;
    Bool_t  bRecalcClusPos  = kFALSE;
    Bool_t  bNonLinearCorr  = kFALSE;
    Bool_t  bRemExoticCell  = kFALSE;
    Bool_t  bRemExoticClus  = kFALSE;
    Bool_t  bFidRegion      = kFALSE;
    Bool_t  bCalibEnergy    = kTRUE;
    Bool_t  bCalibTime      = kTRUE;
    Bool_t  bRemBC          = kTRUE;
    UInt_t  iNonLinFunct    = 0;
    Bool_t  bReclusterize   = kFALSE;
    Float_t fSeedThresh     = 0.1;      // 100 MeV
    Float_t fCellThresh     = 0.05;     // 50 MeV 
    UInt_t  iClusterizer    = 0;
    Bool_t  bTrackMatch     = kFALSE;
    Bool_t  bUpdateCellOnly = kFALSE;
    Float_t fTimeMin        = -50e-6;   // minimum time of physical signal in a cell/digit
    Float_t fTimeMax        =  50e-6;   // maximum time of physical signal in a cell/digit
    Float_t fTimeCut        =  25e-6;
    const char *cPass       = 0;

    AliAnalysisTaskSE *pTenderTask = AddTaskEMCALTender(bDistBC, bRecalibClus, bRecalcClusPos, bNonLinearCorr, bRemExoticCell, bRemExoticClus,
                                                        bFidRegion, bCalibEnergy, bCalibTime, bRemBC, iNonLinFunct, bReclusterize, fSeedThresh,
                                                        fCellThresh, iClusterizer, bTrackMatch, bUpdateCellOnly, fTimeMin, fTimeMax, fTimeCut, cPass);
    pTenderTask->SelectCollisionCandidates(kPhysSel);
  }
  
  if (bDoReclusterize) {
    AliAnalysisTaskEMCALClusterizeFast *pClusterizerTask = AddTaskClusterizerFast("ClusterizerFast", "", "", kClusterizerType, 
                                                                                  0.05, 0.1, kEMCtimeMin, kEMCtimeMax, kEMCtimeCut, 
                                                                                  kTRUE, kFALSE, AliAnalysisTaskEMCALClusterizeFast::kFEEData);
    
    pClusterizerTask->SelectCollisionCandidates(kPhysSel);
  }

  if (bDoFullJets) {
    AliEmcalClusterMaker *pClusterMakerTask = AddTaskEmcalClusterMaker(AliEMCALRecoUtils::kBeamTestCorrected, kFALSE, 0, sClusName, 0.15, kFALSE);
    pClusterMakerTask->SelectCollisionCandidates(kPhysSel);
  }
    
  if (1) {
    if (bDoTender || bDoReclusterize) {
      // QA task
      AliAnalysisTaskSAQA *pQATaskAfter = AddTaskSAQA("", sClusName, sCellName, "", "",
						      0, 0, 0, 0.15, 0.15, "TPC", "AliAnalysisTaskSAQA_AfterTender");
      pQATaskAfter->SetHistoBins(150, 0, 150);
      pQATaskAfter->SelectCollisionCandidates(kPhysSel);
    }
  }
  
  if (bDoHadCorr) {    
    sCorrClusName = "CaloClustersCorr";
    
    TString sEmcalTracksName("EmcalTracks_");
    TString sEmcalClusName("EmcalClusters_");

    sEmcalTracksName += sTracksName;
    sEmcalClusName += sClusName;
    
    // EmcalParticles maker
    AliEmcalParticleMaker *pPartMakerTask = AddTaskEmcalParticleMaker(sTracksName, sClusName, sEmcalTracksName, sEmcalClusName, "AliEmcalParticleMaker");
    pPartMakerTask->SelectCollisionCandidates(kPhysSel);

    // Cluster-track matcher task
    AliEmcalClusTrackMatcherTask *pMatcherTask = AddTaskEmcalClusTrackMatcher(sEmcalTracksName, sEmcalClusName, 0.1, kFALSE, kTRUE);
    pMatcherTask->SelectCollisionCandidates(kPhysSel);
      
    // Hadronic correction task
    AliHadCorrTask *pHadCorrTask = AddTaskHadCorr(sEmcalTracksName, sEmcalClusName, sCorrClusName, 
                                                  kHadCorrF, 0.15, 0.030, 0.015, 0, kTRUE, kTRUE);
    pHadCorrTask->SetHistoBins(150,0,150);
    pHadCorrTask->SelectCollisionCandidates(kPhysSel);
  }
  else {
    sCorrClusName = sClusName;
  }
  
  if (1) {
    // QA task
    AliAnalysisTaskSAQA *pQATask = 0;

    if (bDoFullJets) {
      pQATask = AddTaskSAQA(sTracksName, sCorrClusName, sCellName, "", "", 0.2, 1, 0, kTrackPtCut, kClusPtCut, "TPC");
    }
    else {
      pQATask = AddTaskSAQA(sTracksName, "", "", "", "", 0.2, 1, 0, kTrackPtCut, kClusPtCut, "TPC");
    }
    
    pQATask->SelectCollisionCandidates(kPhysSel);
    pQATask->SetHistoBins(200, 0, 30);
  }

  // Charged jet analysis
  if (bDoChargedJets) {
    AliEmcalJetTask *pChJetTask = AddTaskEmcalJet(sTracksName, "", 1, kJetRadius, 1, kTrackPtCut, kClusPtCut, kGhostArea);
    pChJetTask->SelectCollisionCandidates(kPhysSel);
    sChJetsName = pChJetTask->GetName();

    AliAnalysisTaskSAJF *pSpectraChTask = AddTaskSAJF(sTracksName, "", sChJetsName, "",  kJetRadius, kJetPtCut, kJetAreaCut, "TPC");
    pSpectraChTask->SetNLeadingJets(1);
    pSpectraChTask->SelectCollisionCandidates(kPhysSel);
    pSpectraChTask->SetHistoType(kHistoType);
  }

  // Full jet analysis
  if (bDoFullJets) {
    AliEmcalJetTask *pFuJetTask = AddTaskEmcalJet(sTracksName, sCorrClusName, 1, kJetRadius, 0, kTrackPtCut, kClusPtCut, kGhostArea);
    pFuJetTask->SelectCollisionCandidates(kPhysSel);   
    sFuJetsName = pFuJetTask->GetName();

    AliAnalysisTaskSAJF *pSpectraFuTask = AddTaskSAJF(sTracksName, sCorrClusName, sFuJetsName, "", kJetRadius, kJetPtCut, kJetAreaCut, "EMCAL"); 
    pSpectraFuTask->SetNLeadingJets(1);
    pSpectraFuTask->SelectCollisionCandidates(kPhysSel);
    pSpectraFuTask->SetHistoType(kHistoType);
  }

  // HF-jet analysis
  if (bDoHF) {
    AliAnalysisTaskSEDmesonsFilterCJ* pDStarMesonFilterTask = AddTaskSEDmesonsFilterCJ(AliAnalysisTaskSEDmesonsFilterCJ::kDstartoKpipi,
                                                                                       "DStartoKpipiCuts.root",
                                                                                       kFALSE,  //   Bool_t theMCon
                                                                                       kTRUE,   //   Bool_t reco
                                                                                       "");
    pDStarMesonFilterTask->SelectCollisionCandidates(kPhysSel);

    AliAnalysisTaskSEDmesonsFilterCJ* pD0mesonFilterTask = AddTaskSEDmesonsFilterCJ(AliAnalysisTaskSEDmesonsFilterCJ::kD0toKpi,
                                                                                   "DStartoKpipiCuts.root",
                                                                                   kFALSE,  //   Bool_t theMCon
                                                                                   kTRUE,   //   Bool_t reco
                                                                                   "");
    pD0mesonFilterTask->SelectCollisionCandidates(kPhysSel);
  }
  
  if (1) {
    TObjArray *pTopTasks = pMgr->GetTasks();
    for (Int_t i = 0; i < pTopTasks->GetEntries(); ++i) {
      AliAnalysisTaskSE *pTask = dynamic_cast<AliAnalysisTaskSE*>(pTopTasks->At(i));
      if (!pTask) continue;
      if (pTask->InheritsFrom("AliAnalysisTaskEmcal")) {
        AliAnalysisTaskEmcal *pTaskEmcal = static_cast<AliAnalysisTaskEmcal*>(pTask);
        if (bForcePP) {
          Printf("Setting beam type for task %s", pTaskEmcal->GetName());
          pTaskEmcal->SetForceBeamType(0);
        }
      }
    }
  }
}

void LoadMacros()
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskTrackingQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAJF.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskSEDmesonsFilterCJ.C");
}
