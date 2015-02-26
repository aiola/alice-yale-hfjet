// $Id$

void AddTaskJetResp(Bool_t pdsftrain = kTRUE, const char *datatype = "AOD", const char *runtype = "local",
		    Bool_t dotpconly = kFALSE, Bool_t doemcal = kTRUE, Bool_t doemb = kFALSE,
		    Bool_t dohadcorr = kTRUE)
{
  enum eDataType { kAod, kEsd };
  enum eRunType  { kLocal, kGrid };

  eRunType rType;
  if (!strcmp(runtype, "grid")) 
    rType = kGrid;
  else if (!strcmp(runtype, "local")) 
    rType = kLocal;
  else {
    cout << "Incorrect run option, check first argument of run macro" << endl;
    cout << "runtype = local or grid" << endl;
    return;
  }
  
  eDataType dType;
  if (!strcmp(datatype, "ESD"))
    dType = kEsd;
  else if (!strcmp(datatype, "AOD"))
    dType = kAod;
  else {
    cout << "Incorrect data type option, check third argument of run macro" << endl;
    cout << "datatype = AOD or ESD" << endl;
    return;
  }

  if (pdsftrain) {
    TString period(gSystem->Getenv("ETRAIN_PERIOD"));
    if (!period.Contains("lhc12a15a")) {
      cout << "AddTaskJetResp.C ignored because dataset is not lhc12a15a" << endl;
      return;
    }
  }

  const char *tracksName    = "PicoTracks";
  const char *clusName      = "EmcCaloClusters";
  const char *mcTracksName  = "mcparticles";
  Double_t jetRadius = 0.2;

  char cellsName[100];
  if (dType == kAod) {
    strcpy(cellsName, "emcalCells");
  }
  else {
    strcpy(cellsName, "EMCALCells");
  }
  
  char corrClusName[100];

  Float_t     trackPtCut          = 0.15;
  Float_t     clusPtCut           = 0.30;
  Float_t     partLevPtCut        = 0.;
  Float_t     jetBias             = 0;
  Int_t       biasType            = 1;   //  0 = charged, 1 = neutral, 2 = both
  Float_t     jetPtCut            = 1;
  Float_t     jetAreaCut          = 0.;
  Float_t     hadCorr             = 2;
  Int_t       histoType           = 1;  // 0 = TH2, 1 = THnSparse
  Float_t     propDist            = 440.;
  Bool_t      selectHijingPart    = kFALSE;
  Bool_t      isPythia            = kTRUE;
  Float_t     trackingEff         = 1.;
  Bool_t      forcePP             = kTRUE;

  TF1* trackingEff_function = new TF1("f","[0] * (1 - [1] / (0.5+x))",0.5,15);
  trackingEff_function->SetParameter(0,9.86064e-01);
  trackingEff_function->SetParameter(1,2.46923e-02);

  Double_t kRhoMinEta = -0.5;
  Double_t kRhoMaxEta = 0.5;
  Double_t kRhoMinPhi = 30*TMath::DegToRad()+0.2;
  Double_t kRhoMaxPhi = 230*TMath::DegToRad()-0.2;

  UInt_t   matching = 1; //1=geometrical, 2=MClabel
  Double_t matchingLevel = 0.99;

  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPicoTrackMaker.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbeddingFromAOD.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskClusterizerFast.C"); 
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); 
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskHadCorr.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskSAQA.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRho.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskDeltaPt.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskSAJF.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskTrackingQA.C");

  if (dType == kEsd && strcmp(mcTracksName,"")!=0) {
    // MC particle selector
    AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(mcTracksName, kFALSE, kFALSE, 1);
    mcPartTask->SetOnlyPhysPrim(kFALSE);
    mcPartTask->SetOnlyHIJING(kFALSE);
  }

  if (dType == kEsd || dType == kAod) {
    const char *inputTracks = "HybridTracks";

    if (dType == kEsd) {
      // Hybrid tracks maker for ESD
      AliEmcalEsdTrackFilterTask *hybTask = AddTaskEmcalEsdTrackFilter(inputTracks, "Hybrid_LHC12a15e");
      //hybTask->SetTrackEfficiency(trackingEff);
      hybTask->SetTrackEfficiency(trackingEff_function);
      hybTask->SetDoPropagation(kTRUE);
      hybTask->SetDist(propDist);
      hybTask->SetMC(kTRUE);
    }
    else if (dType == kAod) {
      // Hybrid tracks maker for AOD
      AliEmcalAodTrackFilterTask *hybTask = AddTaskEmcalAodTrackFilter(inputTracks, "tracks", "LHC12a15e");
      //hybTask->SetTrackEfficiency(trackingEff);
      hybTask->SetTrackEfficiency(trackingEff_function);
      hybTask->SetMC(kTRUE);
    }

    // PicoTracks maker
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(tracksName, inputTracks);
    if (strcmp(mcTracksName,"")!=0) pTrackTask->SetCopyMCFlag(kTRUE, mcTracksName);
  }

  if (1) {
    AliEmcalTrackingQATask* trackingQAtask = AddTaskTrackingQA(tracksName, mcTracksName, selectHijingPart);
  }

  if (doemb) {
    // Embedding
    AliJetEmbeddingFromAODTask *embFromAOD = AddTaskJetEmbeddingFromAOD(tracksName, "", cellsName, "", "files_LHC11h_2_AOD145.txt",
									"aodTree", "tracks", "", "emcalCells", "", "lhc11h", kFALSE, 
									0, 10, AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral,
                                                                        kFALSE, kTRUE, 1234567890);
    embFromAOD->SetEmbedCentrality(kTRUE);
  }

  if (doemcal) {
    // Tender Supplies
    if (1) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEMCALTender.C");

      //----------------------- Add tender -------------------------------------------------------
      Bool_t distBC         = kFALSE; //switch for recalculation cluster position from bad channel 
      Bool_t recalibClus    = kFALSE;
      Bool_t recalcClusPos  = kFALSE;
      Bool_t nonLinearCorr  = kFALSE;
      Bool_t remExoticCell  = kFALSE;
      Bool_t remExoticClus  = kFALSE;
      Bool_t fidRegion      = kFALSE;
      Bool_t calibEnergy    = kFALSE;
      Bool_t calibTime      = kFALSE;
      Bool_t remBC          = kTRUE;
      UInt_t nonLinFunct    = 0;
      Bool_t reclusterize   = kFALSE;
      Float_t seedthresh    = 0.1;      // 100 MeV
      Float_t cellthresh    = 0.05;     // 50 MeV 
      UInt_t clusterizer    = 0;
      Bool_t trackMatch     = kFALSE;
      Bool_t updateCellOnly = kFALSE;
      Float_t timeMin       = -50e-6;   // minimum time of physical signal in a cell/digit
      Float_t timeMax       =  50e-6;   // maximum time of physical signal in a cell/digit
      Float_t timeCut       =  25e-6;
      const char *pass = "pass2";
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEMCALTender.C");//tendertasks
      AliAnalysisTaskSE *tender = AddTaskEMCALTender(distBC, recalibClus, recalcClusPos, nonLinearCorr, remExoticCell, remExoticClus,
						     fidRegion, calibEnergy, calibTime, remBC, nonLinFunct, reclusterize, seedthresh,
						     cellthresh, clusterizer, trackMatch, updateCellOnly, timeMin, timeMax, timeCut, pass);
      if (0&&dType != kAod) {
	AliTender *alitender = dynamic_cast<AliTender*>(tender);
	if (!alitender) {
	  ::Error("Incorrect tender type!");
	  return;
	} 
	if (rType == kLocal) 
	  alitender->SetDefaultCDBStorage("local://$ALICE_ROOT/OCDB"); 
	else if (rType == kPdsf || rType == kCarver) 
	  alitender->SetDefaultCDBStorage("local://$ALICE_OCDB"); 
      }
    }

    AliAnalysisTaskEMCALClusterizeFast *taskClust = AddTaskClusterizerFast("ClusterizerFast", "", "", AliEMCALRecParam::kClusterizerv2, 0.05, 0.1, -50e-6, 50e-6, 25e-6, 
									   kTRUE, kFALSE, AliAnalysisTaskEMCALClusterizeFast::kFEEData);
    AliAnalysisTaskSE *taskPrep = AddTaskEmcalClusterMaker(AliEMCALRecoUtils::kPi0MCv3, kFALSE, 0, clusName, 0.15, kFALSE);

    if (dohadcorr) {
      strcpy(corrClusName, "CaloClustersCorr");
      TString emcalTracksName(Form("EmcalTracks_%s",tracksName));
      TString emcalClusName(Form("EmcalClusters_%s",clusName));
      
      // EmcalParticles maker
      AliEmcalParticleMaker *ePartTask = AddTaskEmcalParticleMaker(tracksName, clusName, emcalTracksName, emcalClusName, "AliEmcalParticleMaker");
      
      // Cluster-track matcher task
      AliEmcalClusTrackMatcherTask *matcherTask = AddTaskEmcalClusTrackMatcher(emcalTracksName, emcalClusName, 0.1, kFALSE, kTRUE);
      
      // Hadronic correction task
      AliHadCorrTask *hCorrTask = AddTaskHadCorr(emcalTracksName, emcalClusName, corrClusName, hadCorr, 0.15, 0.030, 0.015, 0, kTRUE, kTRUE, "AnalysisResults.root");
      if (pdsftrain)
	RequestMemory(hCorrTask,200*1024);
    }
    else {
      strcpy(corrClusName, clusName);
    }
  }
 
  // QA tasks
  if (1 && strcmp(mcTracksName,"")!=0) {
    AliAnalysisTaskSAQA *mcqatask = AddTaskSAQA(mcTracksName, "", "", "", "", 0.2, jetPtCut, 0, partLevPtCut, partLevPtCut, "TPC");
    mcqatask->SetHistoBins(300,0,150);
    mcqatask->SetMC(kTRUE);
    mcqatask->SetParticleLevel(kTRUE);
    AliParticleContainer *partCont = mcqatask->GetParticleContainer(0);
    partCont->SelectPhysicalPrimaries(kTRUE);
    partCont->SelectHIJING(selectHijingPart);
    if (pdsftrain)
      RequestMemory(mcqatask,400*1024);
  }

  if (doemcal) {
    AliAnalysisTaskSAQA *qatask1 = 0;
    qatask1 = AddTaskSAQA("", clusName, "", "", "", 0.2, jetPtCut, 0, trackPtCut, clusPtCut, "TPC");
    qatask1->SetMC(kTRUE);
    qatask1->SetHistoBins(300,0,150);
    if (pdsftrain)
      RequestMemory(qatask,100*1024);

    AliAnalysisTaskSAQA *qatask = AddTaskSAQA(tracksName, corrClusName, "", "", "", 0.2, jetPtCut, 0, trackPtCut, clusPtCut, "TPC");
    qatask->SetMC(kTRUE);
    qatask->SetHistoBins(300,0,150);
    AliParticleContainer *partCont = qatask->GetParticleContainer(0);
    partCont->SelectHIJING(selectHijingPart);
    if (pdsftrain)
      RequestMemory(qatask,200*1024);
  }
  else {
    AliAnalysisTaskSAQA *qatask = AddTaskSAQA(tracksName, "", "", "", "", 0.2, jetPtCut, 0, trackPtCut, clusPtCut, "TPC");
    qatask->SetMC(kTRUE);
    qatask->SetHistoBins(300,0,150);
    AliParticleContainer *partCont = qatask->GetParticleContainer(0);
    partCont->SelectHIJING(selectHijingPart);
    if (pdsftrain)
      RequestMemory(qatask,200*1024);
  }

  if (selectHijingPart && 0) {
    AliAnalysisTaskSAQA *qatask = AddTaskSAQA(tracksName, "", "", "", "", 0.2, jetPtCut, 0, trackPtCut, clusPtCut, "TPC", "AliAnalysisTaskSAQA_WithEmb");
    qatask->SetMC(kTRUE);
    qatask->SetHistoBins(300,0,150);
    if (pdsftrain)
      RequestMemory(qatask,200*1024);
  }

  AliEmcalJetTask *chKtJetTask = 0;
  char *chKtJetsName = 0;

  if (doemb) {
    chKtJetTask = AddTaskEmcalJet(tracksName, "", 0, 0.2, 1, trackPtCut, clusPtCut);
    chKtJetsName = chKtJetTask->GetName();
  }
 
  if (dotpconly) {
    // Jet finder task
    AliEmcalJetTask *chJetTaskch = AddTaskEmcalJet(tracksName, "", 1, jetRadius, 1, trackPtCut, clusPtCut, 0.005, 1, "Jet", 1, kFALSE, kFALSE);
    const char *chJetsName = chJetTaskch->GetName();
    AliEmcalJetTask *chMcJetTaskch = AddTaskEmcalJet(mcTracksName, "", 1, jetRadius, 1, partLevPtCut, partLevPtCut, 0.005, 1, "Jet", 1, kTRUE, kFALSE);
    const char *chMcJetsName = chMcJetTaskch->GetName();

    TString chRhoName("");

    if (doemb) {
      chRhoName = "ChargedRho";
      // Rho task charged jet
      AliAnalysisTaskRho *chargedRhoTask = AddTaskRho(chKtJetsName, tracksName, "", chRhoName, 
						      0.2, "TPC", 0.01, 0, 0, 2, kTRUE);
      if (pdsftrain)
	RequestMemory(chargedRhoTask,150*1024);
    }
    
    // Analysis task
    AliJetResponseMaker *taskCharged = AddTaskJetResponseMaker(tracksName, "", chJetsName, chRhoName, jetRadius, mcTracksName, "", chMcJetsName, "", jetRadius,
							      jetPtCut, jetAreaCut, jetBias, biasType, matching, matchingLevel, matchingLevel, "TPC");
    taskCharged->SetHistoType(histoType);
    taskCharged->GetParticleContainer(1)->SelectPhysicalPrimaries(kTRUE);
    if (pdsftrain) 
      RequestMemory(taskCharged,150*1024);
  }
  
  if (doemcal) {
    
    AliEmcalJetTask *jetTask = AddTaskEmcalJet(tracksName, corrClusName, 1, jetRadius, 0, trackPtCut, clusPtCut, 0.005, 1, "Jet", 1, kFALSE, kFALSE);
    const char *jetsName = jetTask->GetName();
    AliEmcalJetTask *mcJetTask = AddTaskEmcalJet(mcTracksName, "", 1, jetRadius, 0, partLevPtCut, partLevPtCut, 0.005, 1, "Jet", 1, kTRUE, kFALSE);
    const char *mcJetsName = mcJetTask->GetName();

    TString chRhoSmallName("");
    TString rhoNameMeth2("");

    if (doemb) {
      // 300 MeV
      TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
      sfunc->SetParameter(2,1.80183);
      sfunc->SetParameter(1,-0.0134124);
      sfunc->SetParameter(0,0.000146031);

      chRhoSmallName = "ChargedRhoSmall";
      rhoNameMeth2 = chRhoSmallName;
      rhoNameMeth2 += "_Scaled";

      AliAnalysisTaskRho *rhoTaskSmall = AddTaskRho(chKtJetsName, tracksName, "", chRhoSmallName, 
						    0.2, "USER", 0.01, 0, 0, 0, kTRUE, chRhoSmallName);
      rhoTaskSmall->SetJetEtaLimits(kRhoMinEta,kRhoMaxEta);
      rhoTaskSmall->SetJetPhiLimits(kRhoMinPhi,kRhoMaxPhi);
      rhoTaskSmall->SetScaleFunction(sfunc);
      rhoTaskSmall->SetClusName(corrClusName);

      if (pdsftrain)
	RequestMemory(rhoTaskSmall,150*1024);
    }

    // Analysis task
    AliJetResponseMaker *taskFull = AddTaskJetResponseMaker(tracksName, corrClusName, jetsName, rhoNameMeth2, jetRadius, mcTracksName, "", mcJetsName, "", jetRadius,
							      jetPtCut, jetAreaCut, jetBias, biasType, matching, matchingLevel, matchingLevel, "EMCAL");
    taskFull->SetHistoType(histoType);
    taskFull->SetCaloCellsName(cellsName);
    taskFull->GetParticleContainer(1)->SelectPhysicalPrimaries(kTRUE);
    if (pdsftrain) 
      RequestMemory(taskFull,150*1024);
  }

  UInt_t physSel = 0;
  //UInt_t physSel = AliEmcalPhysicsSelection::kEmcalOk;
  //UInt_t physSel = AliVEvent::kINT7;
  //UInt_t physSel = AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  //UInt_t physSel = AliVEvent::kEMCEGA;
  //UInt_t physSel = AliVEvent::kEMCEJE;
  //UInt_t physSel = AliEmcalPhysicsSelection::kEmcalHT;

  if (1) {
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      ::Error("AddTaskJetResp", "No analysis manager to connect to.");
      return NULL;
    }  

    TObjArray *toptasks = mgr->GetTasks();
    for (Int_t i = 0; i < toptasks->GetEntries(); ++i) {
      AliAnalysisTaskSE *task = dynamic_cast<AliAnalysisTaskSE*>(toptasks->At(i));
      if (!task)
        continue;

      if (!task->InheritsFrom("AliPhysicsSelectionTask")) {
        task->SelectCollisionCandidates(physSel);
      }

      if (task->InheritsFrom("AliAnalysisTaskEmcal")) {
	AliAnalysisTaskEmcal *taskEmcal = dynamic_cast<AliAnalysisTaskEmcal*>(task);
	if (forcePP) taskEmcal->SetForceBeamType(AliAnalysisTaskEmcal::kpp);
	taskEmcal->SetVzRange(-10,10);
	taskEmcal->SetIsPythia(isPythia);
      }
    }
  }
}
