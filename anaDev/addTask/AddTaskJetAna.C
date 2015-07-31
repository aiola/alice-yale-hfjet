// AddTaskJetAna.C

void LoadMacros();

void AddTaskJetAna(const char *cDataType = "AOD", const char *cRunType = "local",
                   const UInt_t kPhysSel = AliEmcalPhysicsSelection::kEmcalOk,
                   Bool_t bDoHF = kTRUE, Bool_t bDoChargedJets = kFALSE, Bool_t bDoFullJets = kFALSE,
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
  const Double_t kJetRadius           = 0.6;
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
  const Int_t eFlavourJetMatchingType = AliAnalysisTaskDmesonJetCorrelations::kGeometricalMatching;    // kGeometricalMatching, kConstituentMatching, kJetLoop

  TString sTracksName("AODFilterTracks");
  //TString sTracksName("tracks");
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
  else if (iDataType == kAod) { // for the moment disabled
    // Hybrid tracks maker for AOD
    AliEmcalAodTrackFilterTask *pHybTask = AddTaskEmcalAodTrackFilter(sTracksName, "tracks", "LHC10b");
    pHybTask->SelectCollisionCandidates(kPhysSel);
    pHybTask->SetAttemptProp(kFALSE);
  }

  if (bDoTrackingQA) {
    AliEmcalTrackingQATask* pTrackingQAtask = AddTaskTrackingQA(sTracksName, "", kFALSE);
    pTrackingQAtask->SelectCollisionCandidates(kPhysSel);
  }

  if (1) {
    if (bDoTender || bDoReclusterize) {
      // QA task
      AliAnalysisTaskSAQA *pQATaskBefore = AddTaskSAQA("", sOrigClusName, sCellName, "", "",
						      0, 0, 0, 0.15, 0., "TPC", "AliAnalysisTaskSAQA_BeforeTender");
      pQATaskBefore->GetClusterContainer(0)->SetClusECut(0.15);
      pQATaskBefore->SetHistoBins(200, 0, 30);
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

    AliEmcalClusterMaker *pClusterMakerTask = AddTaskEmcalClusterMaker(AliEMCALRecoUtils::kBeamTestCorrected, kFALSE, 0, sClusName, 0.15, kFALSE);
    pClusterMakerTask->SelectCollisionCandidates(kPhysSel);
  }
    
  if (1) {
    if (bDoTender || bDoReclusterize) {
      // QA task
      AliAnalysisTaskSAQA *pQATaskAfter = AddTaskSAQA("", sClusName, sCellName, "", "",
						      0, 0, 0, 0.15, 0., "TPC", "AliAnalysisTaskSAQA_AfterTender");
      pQATaskAfter->GetClusterContainer(0)->SetClusECut(kClusPtCut);
      pQATaskAfter->SetHistoBins(200, 0, 30);
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
  
  if (0) {
    // QA task
    AliAnalysisTaskSAQA *pQATask = AddTaskSAQA(sTracksName, sCorrClusName, sCellName, "", "", 0.2, 1, 0, kTrackPtCut, 0., "TPC");
    pQATask->GetParticleContainer(0)->SetClassName("AliAODTrack");
    //pQATask->GetParticleContainer(0)->SetFilterHybridTracks(kTRUE);
    //pQATask->SetAODfilterBits(256, 512);
    pQATask->GetClusterContainer(0)->SetClusECut(kClusPtCut);
    pQATask->SelectCollisionCandidates(kPhysSel);
    pQATask->SetHistoBins(200, 0, 30);
  }

  // Charged jet analysis
  if (bDoChargedJets) {
    AliEmcalJetTask *pChJetTask = AddTaskEmcalJet(sTracksName, "", 1, kJetRadius, 1, kTrackPtCut, kClusPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, kFALSE);
    pChJetTask->SelectCollisionCandidates(kPhysSel);
    pChJetTask->SetFilterHybridTracks(kTRUE);
    sChJetsName = pChJetTask->GetName();

    AliAnalysisTaskSAJF *pSpectraChTask = AddTaskSAJF(sTracksName, "", sChJetsName, "",  kJetRadius, kJetPtCut, kJetAreaCut, "TPC");
    pSpectraChTask->SetNLeadingJets(1);
    pSpectraChTask->SelectCollisionCandidates(kPhysSel);
    pSpectraChTask->SetHistoType(kHistoType);
  }

  // Full jet analysis
  if (0 && bDoFullJets) {
    AliEmcalJetTask *pFuJetTask = AddTaskEmcalJet(sTracksName, sCorrClusName, 1, kJetRadius, 0, kTrackPtCut, kClusPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, kFALSE);
    pFuJetTask->SelectCollisionCandidates(kPhysSel);   
    sFuJetsName = pFuJetTask->GetName();

    AliAnalysisTaskSAJF *pSpectraFuTask = AddTaskSAJF(sTracksName, sCorrClusName, sFuJetsName, "", kJetRadius, kJetPtCut, kJetAreaCut, "EMCAL"); 
    pSpectraFuTask->SetNLeadingJets(1);
    pSpectraFuTask->SelectCollisionCandidates(kPhysSel);
    pSpectraFuTask->SetHistoType(kHistoType);
  }

  // HF-jet analysis
  if (bDoHF) {
    if (1) {
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
    
      AliEmcalJetTask *pChJetDStarTask = AddTaskEmcalJet(sTracksDStarName, "", 1, kJetRadius, 1, kTrackPtCut, kClusPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, 1);
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
                                                                                      kJetRadius, kJetPtCut, kJetAreaCut, "TPC", 0,
                                                                                      "AliAnalysisTaskDmesonJetCorrelations", "rec");
      pDStarMesonJetCorr->SetMaxR(kJetRadius);
      pDStarMesonJetCorr->SetMatchingType(eFlavourJetMatchingType);
      pDStarMesonJetCorr->SetPlotOnlyAcceptedJets(kTRUE);
      pDStarMesonJetCorr->SetShowDeltaEta(kTRUE);
      pDStarMesonJetCorr->SetShowDeltaPhi(kTRUE);
      pDStarMesonJetCorr->SetShow2ProngInvMass(kTRUE);
      pDStarMesonJetCorr->SetShowInvMass(kTRUE);
      //pDStarMesonJetCorr->SetShowDaughterDistance(3);
      pDStarMesonJetCorr->SetShowJetConstituents(kTRUE);
      pDStarMesonJetCorr->SelectCollisionCandidates(kPhysSel);

      // QA task
      AliAnalysisTaskSAQA *pQADStarTask = AddTaskSAQA(sTracksDStarName, "", "", "", "", 0.2, 1, 0, kTrackPtCut, kClusPtCut, "TPC");
      //pQADStarTask->GetParticleContainer(0)->SetClassName("AliAODTrack");
      //pQADStarTask->GetParticleContainer(0)->SetFilterHybridTracks(kTRUE);
      //pQADStarTask->SetAODfilterBits(256, 512);
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
    
      AliEmcalJetTask *pChJetD0Task = AddTaskEmcalJet(sTracksD0Name, "", 1, kJetRadius, 0, kTrackPtCut, kClusPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, 1);
      pChJetD0Task->SelectCollisionCandidates(kPhysSel);
      TString sChJetsD0Name = pChJetD0Task->GetName();

      AliAnalysisTaskSAJF *pSpectraChD0Task = AddTaskSAJF(sTracksD0Name, "", sChJetsD0Name, "",  kJetRadius, kJetPtCut, kJetAreaCut, "TPC");
      pSpectraChD0Task->SetNLeadingJets(1);
      pSpectraChD0Task->SelectCollisionCandidates(kPhysSel);
      pSpectraChD0Task->SetHistoType(kHistoType);
    
      AliAnalysisTaskDmesonJetCorrelations* pD0MesonJetCorr = AddTaskDmesonJetCorr(AliAnalysisTaskDmesonJetCorrelations::kD0toKpi, "", 
                                                                                   sTracksD0Name, "", sChJetsD0Name, "",
                                                                                   kJetRadius, kJetPtCut, kJetAreaCut, "TPC", 0,
                                                                                   "AliAnalysisTaskDmesonJetCorrelations", "rec");
      pD0MesonJetCorr->SetMaxR(kJetRadius);
      pD0MesonJetCorr->SetMatchingType(eFlavourJetMatchingType);
      pD0MesonJetCorr->SetPlotOnlyAcceptedJets(kTRUE);
      pD0MesonJetCorr->SetShowDeltaEta(kTRUE);
      pD0MesonJetCorr->SetShowDeltaPhi(kTRUE);
      //pD0MesonJetCorr->SetShowDaughterDistance(2);
      pD0MesonJetCorr->SetShowJetConstituents(kTRUE);
      pD0MesonJetCorr->SelectCollisionCandidates(kPhysSel);

      AliAnalysisTaskSAQA *pQAD0Task = AddTaskSAQA(sTracksD0Name, "", "", "", "", 0.2, 1, 0, kTrackPtCut, kClusPtCut, "TPC");
      //pQAD0Task->GetParticleContainer(0)->SetClassName("AliAODTrack");
      //pQAD0Task->GetParticleContainer(0)->SetFilterHybridTracks(kTRUE);
      //pQAD0Task->SetAODfilterBits(256, 512);
      pQAD0Task->SelectCollisionCandidates(kPhysSel);
      pQAD0Task->SetHistoBins(200, 0, 30);
    }
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
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskDmesonJetCorr.C");
}
