TF1* GetScaleFunc(const Double_t kClusterizerType, const Double_t kHadCorrF);
TString GetScaleFuncName(const Double_t kClusterizerType, const Double_t kHadCorrF);
void LoadMacros();

void AddTaskJetAnaSA(Bool_t bIsEmcalTrain = kTRUE, const char *cDataType = "ESD", const char *cRunType = "local",
                     const UInt_t kPhysSel = AliEmcalPhysicsSelection::kEmcalOk, Bool_t bDoEmbedding = kFALSE,
                     Bool_t bDoChargedJets = kTRUE, Bool_t bDoFullJets = kTRUE, 
                     Bool_t bDoTender = kTRUE, Bool_t bDoReclusterize = kTRUE, Bool_t bDoHadCorr = kTRUE)
{
  TString sPeriod(gSystem->Getenv("ETRAIN_DATASET"));
  if (bIsEmcalTrain) {
    if (!sPeriod.BeginsWith("lhc11hs") && !sPeriod.BeginsWith("lhc12a17b")) {
      Printf("AddTaskJetAnaSA.C ignored because dataset is not lhc11hs or lhc12a17b.");
      return;
    }
  }

  enum eDataType { kAod, kEsd, kSesd };
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
  else if (!strcmp(cDataType, "sESD"))
    iDataType = kSesd;
  else if (!strcmp(cDataType, "AOD"))
    iDataType = kAod;
  else {
    Printf("Incorrect data type option, check third argument of run macro.");
    Printf("cDataType = AOD, ESD or sESD (skimmed ESD).");
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
  const Int_t    kNcent               = 1;
  const Double_t kRhoMinEta           = -0.5;
  const Double_t kRhoMaxEta           = 0.5;
  const Double_t kRhoMinPhi           = 30*TMath::DegToRad()+0.2;
  const Double_t kRhoMaxPhi           = 230*TMath::DegToRad()-0.2;
  const Double_t kEmcalMinEta         = -0.7;
  const Double_t kEmcalMaxEta         = 0.7;
  const Double_t kEmcalMinPhi         = 80*TMath::DegToRad();
  const Double_t kEmcalMaxPhi         = 180*TMath::DegToRad();
  const Bool_t   bForcePP             = kFALSE;
  const UInt_t   kMatching            = 1; //1=geometrical, 2=MClabel, 3=same collections
  const Double_t kMaxDistance         = 1.;
  const Int_t    kNcent               = 1;
  const Double_t kTrackEfficiency     = 0.95;
  const Double_t kEMCtimeMin          = -50e-6;
  const Double_t kEMCtimeMax          = 100e-6;
  const Double_t kEMCtimeCut          =  75e-6;

  const Int_t kNpTHardBins = 11;
  //const Int_t kNpTHardBins = 0;
  Double_t kPtHardBinsScaling[11] = {0, 
                                     5.135193e-05, 5.859497e-06, 4.444755e-07, 4.293118e-08, 5.154750e-09, 
                                     6.958612e-10, 1.149828e-10, 2.520137e-11, 6.222240e-12, 2.255832e-12};
  kPtHardBinsScaling[0] = 0;
  kPtHardBinsScaling[1] = 0;
  kPtHardBinsScaling[2] = 0;

  for (Int_t i = 1; i < 11; i++) kPtHardBinsScaling[i] = 1;

  const Double_t kMinPythiaJetPt = 0;

  //TString sPYTHIAPath("alien:///alice/sim/2012/LHC12a15e_fix/%d/%d/AOD149/%04d/AliAOD.root");
  //const Int_t kNAODFiles = 140;
  TString sPYTHIAPath("alien:///alice/sim/2012/LHC12a15e_fix/%d/%d/AOD149/%04d/AliAOD.root");
  const Int_t kNAODFiles = 140;
  TString sPYTHIAProd("lhc12a15e");

  TString sTracksName("PicoTracks");
  TString sClusName("EmcCaloClusters");
  TString sMCTracksName("MCParticlesSelected");

  TString sCellName;
  TString sOrigClusName;
  TString sCorrClusName;

  TString sOrigClusEmbOnlyName("CaloClustersEmbOnly");
  TString sClusEmbOnlyName("EmcCaloClustersEmbOnly");
  TString sCorrClusEmbOnlyName;

  TString sChKtJetsName;
  TString sChRhoName;
  TString sChRhoSmallName;
  TString sFuRhoName;

  TString sChJetsName;
  TString sChEmbOnlyJetsName;
  TString sChGenLevJetsName;
  TString sFuJetsName;
  TString sFuEmbOnlyJetsName;
  TString sFuGenLevJetsName;

  if (iDataType == kAod) {
    sCellName = "emcalCells";
    sOrigClusName = "caloClusters";
  }
  else {
    sCellName = "EMCALCells";
    sOrigClusName = "CaloClusters";
  }

  //TF1* pScaleFunc = GetScaleFunc(kClusterizerType, kHadCorrF);
  //TString sScaleFuncFilePath = "alien:///alice/cern.ch/user/s/saiola/LHC11h_ScaleFactorFunctions.root";
  TString sScaleFuncFilePath = "/Users/saiola/Documents/Work/ALICE/aliemcaldev/saiola/MyMacros/LHC11h_ScaleFactorFunctions.root";
  TString sScaleFuncName = GetScaleFuncName(kClusterizerType, kHadCorrF);

  if (sScaleFuncFilePath.Contains("alien://")) {
    if (!TGrid::Connect("alien://")) {
      Printf("Unable to connect to alien! Returning...");
      return;
    }
  }

  LoadMacros();

  AliAnalysisManager *pMgr = AliAnalysisManager::GetAnalysisManager();
  if (!pMgr) {
    ::Error("AddTaskJetEmbedding", "No analysis manager to connect to.");
    return;
  }  

  if (iDataType == kEsd || iDataType == kAod) {
    const char *cInputTracks = "HybridTracks";

    if (iDataType == kEsd) {
      // Hybrid tracks maker for ESD
      AliEmcalEsdTrackFilterTask *pHybTask = AddTaskEmcalEsdTrackFilter(cInputTracks, "Hybrid_LHC11h");
      pHybTask->SetDoPropagation(kTRUE);
      pHybTask->SetDist(kPropDist);
      pHybTask->SelectCollisionCandidates(kPhysSel);
    }
    else if (iDataType == kAod) {
      // Hybrid tracks maker for AOD
      AliEmcalAodTrackFilterTask *pHybTask = AddTaskEmcalAodTrackFilter(cInputTracks, "tracks", "lhc11h");
      pHybTask->SelectCollisionCandidates(kPhysSel);
      pHybTask->SetAttemptProp(kTRUE);
    }

    // PicoTracks maker
    AliEmcalPicoTrackMaker *pPicoTrackTask = AddTaskEmcalPicoTrackMaker(sTracksName, cInputTracks);
    pPicoTrackTask->SelectCollisionCandidates(kPhysSel);
  }

  if (0) {
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
      if (bIsEmcalTrain) RequestMemory(pQATaskBefore, 100 * 1024);
    }
  }

  if (bDoEmbedding) {
    AliJetEmbeddingFromPYTHIATask *pEmbeddingTask = AddTaskJetEmbeddingFromPYTHIA(sTracksName, "", sCellName, sMCTracksName, 
                                                                                  sPYTHIAPath, kNpTHardBins, kPtHardBinsScaling, 
                                                                                  "aodTree", "tracks", "", "emcalCells", "mcparticles", sPYTHIAProd, 
                                                                                  kFALSE, -1, -1, 0, kMinPythiaJetPt, kFALSE, kTRUE);
    pEmbeddingTask->SelectCollisionCandidates(kPhysSel);
    pEmbeddingTask->SetRandomAccess(kTRUE);
    pEmbeddingTask->SetTrackEfficiency(kTrackEfficiency);
    pEmbeddingTask->SetAODMC(kTRUE);
    pEmbeddingTask->SetIncludeNoITS(kFALSE);
    pEmbeddingTask->SetTotalFiles(kNAODFiles);
    pEmbeddingTask->SetAODfilterBits(256,512);
    pEmbeddingTask->SetMarkMC(99999);
    pEmbeddingTask->SetMinEntriesPerPtHardBin(-1);
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
    if (bIsEmcalTrain) RequestMemory(pTenderTask, 250 * 1024);
  }
  
  if (bDoReclusterize) {
    AliAnalysisTaskEMCALClusterizeFast *pClusterizerTask = AddTaskClusterizerFast("ClusterizerFast", "", "", kClusterizerType, 
                                                                                  0.05, 0.1, kEMCtimeMin, kEMCtimeMax, kEMCtimeCut, 
                                                                                  kTRUE, kFALSE, AliAnalysisTaskEMCALClusterizeFast::kFEEData);
    
    pClusterizerTask->SelectCollisionCandidates(kPhysSel);
    if (bIsEmcalTrain) RequestMemory(pClusterizerTask,100*1024);
  }

  AliEmcalClusterMaker *pClusterMakerTask = AddTaskEmcalClusterMaker(AliEMCALRecoUtils::kBeamTestCorrected, kFALSE, 0, sClusName, 0.15, kFALSE);
  pClusterMakerTask->SelectCollisionCandidates(kPhysSel);
  if (bIsEmcalTrain) RequestMemory(pClusterMakerTask, 100*1024);

  if (bDoEmbedding) {
      AliAnalysisTaskEMCALClusterizeFast *pClusterizerEmbOnlyTask = AddTaskClusterizerFast("ClusterizerFastEmbOnly", "", sOrigClusEmbOnlyName, kClusterizerType, 
                                                                                           0.05, 0.1, kEMCtimeMin, kEMCtimeMax, kEMCtimeCut, 
                                                                                           kTRUE, kFALSE, AliAnalysisTaskEMCALClusterizeFast::kFEEDataMCOnly);
    
      pClusterizerEmbOnlyTask->SelectCollisionCandidates(kPhysSel);
      if (bIsEmcalTrain) RequestMemory(pClusterizerEmbOnlyTask,100*1024);
    
    AliEmcalClusterMaker *pClusterMakerEmbOnlyTask = AddTaskEmcalClusterMaker(AliEMCALRecoUtils::kBeamTestCorrected, kFALSE, 0, sClusEmbOnlyName, 0.15, kFALSE);
    pClusterMakerEmbOnlyTask->GetClusterContainer(0)->SetArrayName(sOrigClusEmbOnlyName);
    pClusterMakerEmbOnlyTask->SetName(Form("EmcalClusterMaker_%s_%s", sOrigClusEmbOnlyName.Data(), sClusEmbOnlyName.Data()));
    pClusterMakerEmbOnlyTask->SelectCollisionCandidates(kPhysSel);
    if (bIsEmcalTrain) RequestMemory(pClusterMakerEmbOnlyTask, 100*1024);
  }
    
  if (1) {
    if (bDoTender || bDoReclusterize) {
      // QA task
      AliAnalysisTaskSAQA *pQATaskAfter = AddTaskSAQA("", sClusName, sCellName, "", "",
						      0, 0, 0, 0.15, 0.15, "TPC", "AliAnalysisTaskSAQA_AfterTender");
      pQATaskAfter->SetHistoBins(150, 0, 150);
      pQATaskAfter->SelectCollisionCandidates(kPhysSel);
      if (bIsEmcalTrain) RequestMemory(pQATaskAfter, 100 * 1024);
    }
  }

  if (1) {
    if ((bDoTender || bDoReclusterize) && bDoEmbedding) {
      // QA task
      AliAnalysisTaskSAQA *pQATaskAfterEmbOnly = AddTaskSAQA("", sClusEmbOnlyName, "", "", "",
                                                             0, 0, 0, 0.15, 0.15, "TPC", "AliAnalysisTaskSAQA_AfterTenderEmbOnly");
      pQATaskAfterEmbOnly->SetHistoBins(150, 0, 150);
      pQATaskAfterEmbOnly->SelectCollisionCandidates(kPhysSel);
      if (bIsEmcalTrain) RequestMemory(pQATaskAfterEmbOnly, 100 * 1024);
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
    pHadCorrTask->SetIsEmbedded(bDoEmbedding);
    pHadCorrTask->SetHistoBins(150,0,150);
    pHadCorrTask->SelectCollisionCandidates(kPhysSel);
    if (bIsEmcalTrain) RequestMemory(pHadCorrTask,400*1024);

    if (bDoEmbedding) {
      sCorrClusEmbOnlyName = "CaloClustersCorrEmbOnly";
    
      TString sEmcalTracksEmbOnlyName("EmcalTracks_EmbOnly");
      TString sEmcalClusEmbOnlyName("EmcalClusters_");

      sEmcalTracksEmbOnlyName += sTracksName;
      sEmcalClusEmbOnlyName += sClusEmbOnlyName;
    
      // EmcalParticles maker
      AliEmcalParticleMaker *pPartMakerEmbOnlyTask = AddTaskEmcalParticleMaker(sTracksName, sClusEmbOnlyName, 
                                                                               sEmcalTracksEmbOnlyName, sEmcalClusEmbOnlyName, "AliEmcalParticleMaker");
      pPartMakerEmbOnlyTask->SelectCollisionCandidates(kPhysSel);

      // Cluster-track matcher task
      AliEmcalClusTrackMatcherTask *pMatcherEmbOnlyTask = AddTaskEmcalClusTrackMatcher(sEmcalTracksEmbOnlyName, sEmcalClusEmbOnlyName, 0.1, kFALSE, kTRUE);
      pMatcherEmbOnlyTask->GetParticleContainer(0)->SetTrackBitMap(TObject::kBitMask);
      pMatcherEmbOnlyTask->GetParticleContainer(1)->SetTrackBitMap(TObject::kBitMask);
      pMatcherEmbOnlyTask->GetParticleContainer(0)->SetMCTrackBitMap(0);
      pMatcherEmbOnlyTask->GetParticleContainer(1)->SetMCTrackBitMap(0);
      pMatcherEmbOnlyTask->SelectCollisionCandidates(kPhysSel);

      // Hadronic correction task
      AliHadCorrTask *pHadCorrEmbOnlyTask = AddTaskHadCorr(sEmcalTracksEmbOnlyName, sEmcalClusEmbOnlyName, sCorrClusEmbOnlyName, 
                                                           kHadCorrF, 0.15, 0.030, 0.015, 0, kTRUE, kTRUE);
      pHadCorrEmbOnlyTask->GetParticleContainer(0)->SetTrackBitMap(TObject::kBitMask);
      pHadCorrEmbOnlyTask->GetParticleContainer(1)->SetTrackBitMap(TObject::kBitMask);
      pHadCorrEmbOnlyTask->SetHistoBins(150,0,150);
      pHadCorrEmbOnlyTask->SetIsEmbedded(kFALSE);
      pHadCorrEmbOnlyTask->SelectCollisionCandidates(kPhysSel);

      if (bIsEmcalTrain) RequestMemory(pHadCorrEmbOnlyTask,400*1024);
    }
  }
  else {
    sCorrClusName = sClusName;
    sCorrClusEmbOnlyName = sClusEmbOnlyName;
  }

  if (1) {
    AliEmcalJetTask *pChKtJetTask = AddTaskEmcalJet(sTracksName, "", 0, 0.2, 1, kTrackPtCut, kClusPtCut, kGhostArea);
    pChKtJetTask->SelectCollisionCandidates(kPhysSel);
    if (bIsEmcalTrain) RequestMemory(pChKtJetTask, 500 * 1024);

    sChKtJetsName = pChKtJetTask->GetName();
    sChRhoName = "Rho_ExLJ";
    sChRhoSmallName = "Rho_Small";
    sFuRhoName = sChRhoSmallName;
    sFuRhoName += "_Scaled";

    // Rho task charged jet
    AliAnalysisTaskRho *pChRhoTask = AddTaskRho(sChKtJetsName, sTracksName, "", sChRhoName, 
                                                0.2, "TPC", 0.01, 0, 0, 2, kTRUE, "Rho_ExLJ");
    pChRhoTask->SetHistoBins(50,0,250);
    pChRhoTask->SelectCollisionCandidates(kPhysSel);
    
    AliAnalysisTaskRho *pChRhoSmallTask = AddTaskRho(sChKtJetsName, sTracksName, sCorrClusName, sChRhoSmallName, 
                                                     0.2, "USER", 0.01, 0, 
                                                     sScaleFuncFilePath, sScaleFuncName,
                                                     0, kTRUE, "Rho_Small");
    pChRhoSmallTask->SetJetEtaLimits(kRhoMinEta,kRhoMaxEta);
    pChRhoSmallTask->SetJetPhiLimits(kRhoMinPhi,kRhoMaxPhi);
    //pChRhoSmallTask->SetScaleFunction(pScaleFunc);
    pChRhoSmallTask->SetHistoBins(50,0,250);
    pChRhoSmallTask->SelectCollisionCandidates(kPhysSel);
  }
  
  if (1) {
    // QA task
    AliAnalysisTaskSAQA *pQATask = AddTaskSAQA(sTracksName, sCorrClusName, sCellName, sChKtJetsName, sChRhoName,
                                               0.2, 1, 0, kTrackPtCut, kClusPtCut, "TPC");
    pQATask->SelectCollisionCandidates(kPhysSel);
    pQATask->SetHistoBins(200, 0, 30);
    pQATask->SetAdditionalCentEst("TRK");
    pQATask->SetDoV0QA(1);
    pQATask->SetDoLeadingObjectPosition(1);
    if (bIsEmcalTrain) RequestMemory(pQATask, 400 * 1024);
  }

  // Scale factor task
  if (1) {
    AliAnalysisTaskScale *pScaleTask = AddTaskScale(sTracksName, sCorrClusName, kTrackPtCut, kClusPtCut, "Scale", 
                                                    sScaleFuncFilePath, sScaleFuncName);
    //pScaleTask->SetScaleFunction(pScaleFunc);
    pScaleTask->SelectCollisionCandidates(kPhysSel);
  }

  // Charged jet analysis
  if (bDoChargedJets) {
    AliEmcalJetTask *pChJetTask = AddTaskEmcalJet(sTracksName, "", 1, kJetRadius, 1, kTrackPtCut, kClusPtCut, kGhostArea);
    pChJetTask->SelectCollisionCandidates(kPhysSel);
    sChJetsName = pChJetTask->GetName();
    if (bIsEmcalTrain) RequestMemory(pChJetTask, 200 * 1024);

    if (bDoEmbedding) {
      AliEmcalJetTask *pChEmbOnlyJetTask = AddTaskEmcalJet(sTracksName, "", 1, kJetRadius, 1, kTrackPtCut, kClusPtCut, kGhostArea, 1, "JetEmbOnly");
      pChEmbOnlyJetTask->SelectCollisionCandidates(kPhysSel);
      pChEmbOnlyJetTask->SelectConstituents(TObject::kBitMask, 0);
      sChEmbOnlyJetsName = pChEmbOnlyJetTask->GetName();
      if (bIsEmcalTrain) RequestMemory(pChEmbOnlyJetTask, 200 * 1024);

      AliEmcalJetTask *pChGenLevJetTask = AddTaskEmcalJet(sMCTracksName, "", 1, kJetRadius, 1, kPartLevPtCut, kPartLevPtCut, kGhostArea);
      pChGenLevJetTask->SelectCollisionCandidates(kPhysSel);
      sChGenLevJetsName = pChGenLevJetTask->GetName();
      if (bIsEmcalTrain) RequestMemory(pChGenLevJetTask, 200 * 1024);
    }

    // DeltaPt task
    if (1) {
      AliAnalysisTaskDeltaPt *pDeltaPtChTask = AddTaskDeltaPt(sTracksName, "", "", 
							      "", "", "", 
							      "", "", sChRhoName, 
							      kJetRadius, kJetAreaCut, 
							      kTrackPtCut, kClusPtCut, "TPC");
      pDeltaPtChTask->SetMCJetPtThreshold(1);
      pDeltaPtChTask->SetHistoBins(300,0,150);
      pDeltaPtChTask->SelectCollisionCandidates(kPhysSel);
      if (bIsEmcalTrain) RequestMemory(pDeltaPtChTask, 100 * 1024);
    }
    
    if (1) {
      AliAnalysisTaskSAJF *pSpectraChTask = AddTaskSAJF(sTracksName, "", sChJetsName, sChRhoName, 
                                                        kJetRadius, kJetPtCut, kJetAreaCut, "TPC");
      pSpectraChTask->SetNLeadingJets(1);
      pSpectraChTask->SelectCollisionCandidates(kPhysSel);
      pSpectraChTask->SetHistoType(kHistoType);
      if (bIsEmcalTrain) RequestMemory(pSpectraChTask, 100 * 1024);
    }

    if (bDoEmbedding) {
      if (1) {
        AliJetResponseMaker *pChRespTask = AddTaskJetRespPtHard(sTracksName, "", sChJetsName, sChRhoName, kJetRadius, 
                                                                sTracksName, "", sChEmbOnlyJetsName, "", kJetRadius, 
                                                                kJetPtCut, kJetAreaCut, 0., 0, 
                                                                kMatching, kMaxDistance, kMaxDistance, "TPC", -999, -999, kNcent);
        for (Int_t i = 0; i < kNcent; i++) {
          pChRespTask[i].SelectCollisionCandidates(kPhysSel);
          pChRespTask[i].GetJetContainer(0)->SetIsParticleLevel(kFALSE);
          pChRespTask[i].GetJetContainer(1)->SetIsParticleLevel(kFALSE);
          pChRespTask[i].GetJetContainer(1)->SetMaxTrackPt(100);
          pChRespTask[i].SetIsPythia(kFALSE);
          pChRespTask[i].SetIsEmbedded(kTRUE);
          pChRespTask[i].SetHistoType(kHistoType);
        }
      }

      if (1) {
        AliJetResponseMaker *pChGenLevRespTask = AddTaskJetRespPtHard(sTracksName, "", sChJetsName, sChRhoName, kJetRadius, 
                                                                      sMCTracksName, "", sChGenLevJetsName, "", kJetRadius, 
                                                                      kJetPtCut, kJetAreaCut, 0., 0, 
                                                                      kMatching, kMaxDistance, kMaxDistance, "TPC", -999, -999, kNcent);
        for (Int_t i = 0; i < kNcent; i++) {
          pChGenLevRespTask[i].SelectCollisionCandidates(kPhysSel);
          pChGenLevRespTask[i].SetIsPythia(kFALSE);
          pChGenLevRespTask[i].SetIsEmbedded(kTRUE);
          pChGenLevRespTask[i].SetHistoType(kHistoType);
        }
      }
    }
  }

  // Full jet analysis
  if (bDoFullJets) {
    AliEmcalJetTask *pFuJetTask = AddTaskEmcalJet(sTracksName, sCorrClusName, 1, kJetRadius, 0, kTrackPtCut, kClusPtCut, kGhostArea);
    pFuJetTask->SelectCollisionCandidates(kPhysSel);
    AliEmcalJetUtilityGenSubtractor* genSub = pFuJetTask->AddUtility(new AliEmcalJetUtilityGenSubtractor("GenSubtractor"));
    genSub->SetUseExternalBkg(kTRUE);
    genSub->SetRhoName(sFuRhoName);
    genSub->SetRhomName(sFuRhoName);
    AliEmcalJetUtilityConstSubtractor* constSub = pFuJetTask->AddUtility(new AliEmcalJetUtilityConstSubtractor("ConstSubtractor"));
    constSub->SetUseExternalBkg(kTRUE);
    constSub->SetRhoName(sFuRhoName);
    constSub->SetJetsSubName("JetsSub");
    constSub->SetParticlesSubName("PartSub");
    
    sFuJetsName = pFuJetTask->GetName();
    if (bIsEmcalTrain) RequestMemory(pFuJetTask, 200 * 1024);

    if (bDoEmbedding) {
      AliEmcalJetTask *pFuEmbOnlyJetTask = AddTaskEmcalJet(sTracksName, sCorrClusEmbOnlyName, 1, kJetRadius, 0, kTrackPtCut, kClusPtCut, kGhostArea, 1, "JetEmbOnly");
      pFuEmbOnlyJetTask->SelectConstituents(TObject::kBitMask, 0);
      pFuEmbOnlyJetTask->SelectCollisionCandidates(kPhysSel);
      sFuEmbOnlyJetsName = pFuEmbOnlyJetTask->GetName();
      if (bIsEmcalTrain) RequestMemory(pFuEmbOnlyJetTask, 200 * 1024);

      AliEmcalJetTask *pFuGenLevJetTask = AddTaskEmcalJet(sMCTracksName, "", 1, kJetRadius, 0, kPartLevPtCut, kPartLevPtCut, kGhostArea);
      pFuGenLevJetTask->SelectCollisionCandidates(kPhysSel);
      sFuGenLevJetsName = pFuGenLevJetTask->GetName();
      if (bIsEmcalTrain) RequestMemory(pFuGenLevJetTask, 200 * 1024);
    }

    // DeltaPt task
    if (1) {
      AliAnalysisTaskDeltaPt *pDeltaPtFuTask = AddTaskDeltaPt(sTracksName, sCorrClusName, "", 
                                                              "", "", "", 
                                                              "", "", sFuRhoName, 
                                                              kJetRadius, kJetAreaCut, 
                                                              kTrackPtCut, kClusPtCut, "EMCAL"); 

      pDeltaPtFuTask->SetMCJetPtThreshold(1);
      pDeltaPtFuTask->SelectCollisionCandidates(kPhysSel);
      if (bIsEmcalTrain) RequestMemory(pDeltaPtFuTask, 100 * 1024);
    }

    if (1) {
      AliAnalysisTaskSAJF *pSpectraFuTask = AddTaskSAJF(sTracksName, sCorrClusName, sFuJetsName, sFuRhoName, 
                                                        kJetRadius, kJetPtCut, kJetAreaCut, 
                                                        "EMCAL"); 
      pSpectraFuTask->SetNLeadingJets(1);
      pSpectraFuTask->SelectCollisionCandidates(kPhysSel);
      pSpectraFuTask->SetHistoType(kHistoType);
      if (bIsEmcalTrain) RequestMemory(pSpectraFuTask, 150 * 1024);
    }

    if (bDoEmbedding) {
      if (1) {
        AliJetResponseMaker *pFuRespTask = AddTaskJetRespPtHard(sTracksName, sCorrClusName, sFuJetsName, sFuRhoName, kJetRadius, 
                                                                sTracksName, sCorrClusEmbOnlyName, sFuEmbOnlyJetsName, "", kJetRadius, 
                                                                kJetPtCut, kJetAreaCut, 0., 0, 
                                                                kMatching, kMaxDistance, kMaxDistance, "EMCAL", -999, -999, kNcent);
        for (Int_t i = 0; i < kNcent; i++) {
          pFuRespTask[i].SelectCollisionCandidates(kPhysSel);
          pFuRespTask[i].GetJetContainer(0)->SetIsParticleLevel(kFALSE);
          pFuRespTask[i].GetJetContainer(1)->SetIsParticleLevel(kFALSE);
          pFuRespTask[i].GetJetContainer(1)->SetMaxTrackPt(100);
          pFuRespTask[i].SetIsPythia(kFALSE);
          pFuRespTask[i].SetIsEmbedded(kTRUE);
          pFuRespTask[i].SetHistoType(kHistoType);
        }
      }

      if (1) {
        AliJetResponseMaker *pFuGenLevRespTask = AddTaskJetRespPtHard(sTracksName, sCorrClusName, sFuJetsName, sFuRhoName, kJetRadius, 
                                                                      sMCTracksName, "", sFuGenLevJetsName, "", kJetRadius, 
                                                                      kJetPtCut, kJetAreaCut, 0., 0, 
                                                                      kMatching, kMaxDistance, kMaxDistance, "EMCAL", -999, -999, kNcent);
        for (Int_t i = 0; i < kNcent; i++) {
          pFuGenLevRespTask[i].SelectCollisionCandidates(kPhysSel);
          pFuGenLevRespTask[i].SetIsPythia(kFALSE);
          pFuGenLevRespTask[i].SetIsEmbedded(kTRUE);
          pFuGenLevRespTask[i].SetHistoType(kHistoType);
        }
      }
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
        if (bIsEmcalTrain) RequestMemory(task, 250 * 1024);	
      }
      else if (pTask->InheritsFrom("AliEmcalJetTask")) {
        if (bIsEmcalTrain) RequestMemory(pTask, 200 * 1024);	
      }
    }
  }
}

TF1* GetScaleFunc(const Double_t kClusterizerType, const Double_t kHadCorrF)
{
  TF1* pScaleFunc = new TF1("oScaleFunc","[0]*x*x+[1]*x+[2]",-1,100);
  
  if (kHadCorrF > 1.9) {
    if (kClusterizerType == AliEMCALRecParam::kClusterizerv2) {
      // 300 MeV, had corr = 2
      pScaleFunc->SetParameter(2,1.8078);
      pScaleFunc->SetParameter(1,-0.0115442);
      pScaleFunc->SetParameter(0,0.000109722);
    }
    else {
      // 3x3, 300 MeV, had corr = 2
      pScaleFunc->SetParameter(2, 2.02415);
      pScaleFunc->SetParameter(1, -8.13464e-03);
      pScaleFunc->SetParameter(0,2.28091e-05);
    }
  }
  else if (kHadCorrF > 1.6) {
    // 300 MeV, had corr = 1.7
    pScaleFunc->SetParameter(2, 1.79148);
    pScaleFunc->SetParameter(1, -1.13795e-02);
    pScaleFunc->SetParameter(0, 1.14528e-04);
  }
  else if (kHadCorrF > 1.2) {
    // 300 MeV, had corr = 1.3
    pScaleFunc->SetParameter(2,1.90655);
    pScaleFunc->SetParameter(1,-0.0159274);
    pScaleFunc->SetParameter(0,0.0001706891);
  }
  else {
    // 300 MeV, had corr = 0
    pScaleFunc->SetParameter(2,2.01265);
    pScaleFunc->SetParameter(1,-0.0128864);
    pScaleFunc->SetParameter(0,0.000100894);
  }
  return pScaleFunc;
}

TString GetScaleFuncName(const Double_t kClusterizerType, const Double_t kHadCorrF)
{  
  TString name("LHC11h");

  if (kHadCorrF > 1.9) {
    name += "_HadCorr20";
  }
  else if (kHadCorrF > 1.6) {
    name += "_HadCorr17";
  }
  else if (kHadCorrF > 1.2) {
    name += "_HadCorr13";
  }
  else {
    name += "_HadCorr00";
  }

  if (kClusterizerType == AliEMCALRecParam::kClusterizerv2) {
    name += "_ClustersV2";    
  }
  else {
    name += "_Clusters3x3";
  }

  return name;
}

void LoadMacros()
{
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbeddingFromPYTHIA.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEMCALTender.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskClusterizerFast.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); 
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskHadCorr.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskScale.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRho.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskDeltaPt.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskSAQA.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskSAJF.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetRespPtHard.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskTrackingQA.C");
}
