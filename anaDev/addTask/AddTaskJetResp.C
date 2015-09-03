// AddTaskJetResp.C

void AddTaskJetResp(const char *datatype = "AOD", const char *runtype = "local",
		    Bool_t dohf = kTRUE, Bool_t dotpconly = kFALSE, Bool_t doemcal = kFALSE,
		    Bool_t dohadcorr = kFALSE, Bool_t doTrackingQA = kFALSE)
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

  const char *clusName      = "EmcCaloClusters";
  const char *mcTracksName  = "mcparticles";
  //const char *tracksName    = "HybridTracks";
  const char *tracksName    = "tracks";
  Double_t jetRadius = 0.4;

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
  Int_t       biasType            = 0;   //  0 = charged, 1 = neutral, 2 = both
  Float_t     jetPtCut            = 1;
  Float_t     jetAreaCut          = 0.;
  Float_t     hadCorr             = 2;
  Int_t       histoType           = 1;  // 0 = TH2, 1 = THnSparse
  Float_t     propDist            = 440.;
  Bool_t      selectHijingPart    = kFALSE;
  Bool_t      isPythia            = kTRUE;
  Float_t     trackingEff         = 1.;
  Bool_t      forcePP             = kTRUE;
  Double_t    kGhostArea          = 0.01;
  Int_t eFlavourJetMatchingType   = AliAnalysisTaskDmesonJetCorrelations::kGeometricalMatching;    // kGeometricalMatching, kConstituentMatching, kJetLoop

  AliAnalysisTaskSEDmesonsFilterCJ::ECandidateType kDmesonType = AliAnalysisTaskSEDmesonsFilterCJ::kD0toKpi;
  AliAnalysisTaskDmesonJetCorrelations::ECandidateType kDmesonCorrType = AliAnalysisTaskDmesonJetCorrelations::kD0toKpi; //kDstartoKpipi  
  AliEmcalJet::EFlavourTag kFlavourCut = AliEmcalJet::kD0;
  TString sDmesonCandName = "DSBcandidates";
  TString sTracksDcandidatesName(sDmesonCandName);
  sTracksDcandidatesName += "AndTracksD0MCrec";
  TString mcTracksDMesonname = "mcparticlesD0";

  TF1* trackingEff_function = 0;

  Double_t kRhoMinEta = -0.7+jetRadius;
  Double_t kRhoMaxEta = 0.7-jetRadius;
  Double_t kRhoMinPhi = 30*TMath::DegToRad()+jetRadius;
  Double_t kRhoMaxPhi = 230*TMath::DegToRad()-jetRadius;

  UInt_t   matching = 1; //1=geometrical, 2=MClabel
  Double_t matchingLevel = 0.99;

  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalEsdTrackFilter.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C"); 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C"); 
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAJF.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskTrackingQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskMCHFParticleSelector.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskSEDmesonsFilterCJ.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskDmesonJetCorr.C");

  if (dType == kEsd && strcmp(mcTracksName,"")!=0) {
    // MC particle selector
    AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(mcTracksName, kFALSE, kFALSE, 1);
    mcPartTask->SetOnlyPhysPrim(kTRUE);
    mcPartTask->SetOnlyHIJING(kFALSE);
  }

  if (0 && (dType == kEsd || dType == kAod)) {
    if (dType == kEsd) {
      // Hybrid tracks maker for ESD
      AliEmcalEsdTrackFilterTask *hybTask = AddTaskEmcalEsdTrackFilter(tracksName, "Hybrid_LHC10f7a");
      //hybTask->SetTrackEfficiency(trackingEff);
      hybTask->SetTrackEfficiency(trackingEff_function);
      hybTask->SetDoPropagation(dohadcorr);
      hybTask->SetDist(propDist);
      hybTask->SetMC(kTRUE);
    }
    else if (dType == kAod) {
      // Hybrid tracks maker for AOD
      AliEmcalAodTrackFilterTask *hybTask = AddTaskEmcalAodTrackFilter(tracksName, "tracks", "LHC10f7a");
      //hybTask->SetTrackEfficiency(trackingEff);
      hybTask->SetTrackEfficiency(trackingEff_function);
      hybTask->SetAttemptProp(dohadcorr);
      hybTask->SetDist(propDist);
      hybTask->SetMC(kTRUE);
    }
  }

  if (doTrackingQA) {
    AliEmcalTrackingQATask* trackingQAtask = AddTaskTrackingQA(tracksName, mcTracksName, selectHijingPart);
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
      /*
      if (dType != kAod) {
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
      */
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
    }
    else {
      strcpy(corrClusName, clusName);
    }
  }
 
  // QA tasks
  if (strcmp(mcTracksName,"")!=0) {
    AliAnalysisTaskSAQA *mcqatask = AddTaskSAQA(mcTracksName, "", "", "", "", jetRadius, jetPtCut, 0, partLevPtCut, partLevPtCut, "TPC");
    mcqatask->SetHistoBins(300,0,150);
    mcqatask->SetMC(kTRUE);
    mcqatask->SetParticleLevel(kTRUE);
    AliParticleContainer *partCont = mcqatask->GetParticleContainer(0);
    partCont->SelectPhysicalPrimaries(kTRUE);
    partCont->SelectHIJING(selectHijingPart);
  }

  if (doemcal) {
    AliAnalysisTaskSAQA *qatask1 = 0;
    qatask1 = AddTaskSAQA("", clusName, "", "", "", jetRadius, jetPtCut, 0, trackPtCut, clusPtCut, "TPC");
    qatask1->SetMC(kTRUE);
    qatask1->SetHistoBins(300,0,150);

    if (dohadcorr) {
      AliAnalysisTaskSAQA *qatask = AddTaskSAQA(tracksName, corrClusName, "", "", "", jetRadius, jetPtCut, 0, trackPtCut, clusPtCut, "TPC");
      qatask->SetMC(kTRUE);
      qatask->SetHistoBins(300,0,150);
      AliParticleContainer *partCont = qatask->GetParticleContainer(0);
      partCont->SelectHIJING(selectHijingPart);
    }
  }
  else {
    AliAnalysisTaskSAQA *qatask = AddTaskSAQA(tracksName, "", "", "", "", jetRadius, jetPtCut, 0, trackPtCut, clusPtCut, "TPC");
    qatask->SetMC(kTRUE);
    qatask->SetHistoBins(300,0,150);
    AliParticleContainer *partCont = qatask->GetParticleContainer(0);
    partCont->SelectHIJING(selectHijingPart);
    partCont->SetClassName("AliAODTrack");
    partCont->SetFilterHybridTracks(kTRUE);
    qatask->SetAODfilterBits(256, 512);
  }
 
  if (dotpconly) {
    // Jet finder task
    AliEmcalJetTask *chJetTaskch = AddTaskEmcalJet(tracksName, "", 1, jetRadius, 1, trackPtCut, clusPtCut, kGhostArea, 1, "Jet", 1, kFALSE, kFALSE);
    const char *chJetsName = chJetTaskch->GetName();
    AliEmcalJetTask *chMcJetTaskch = AddTaskEmcalJet(mcTracksName, "", 1, jetRadius, 1, partLevPtCut, partLevPtCut, kGhostArea, 1, "Jet", 1, kTRUE, kFALSE);
    const char *chMcJetsName = chMcJetTaskch->GetName();
    
    // Analysis task
    AliJetResponseMaker *taskCharged = AddTaskJetResponseMaker(tracksName, "", chJetsName, "", jetRadius, mcTracksName, "", chMcJetsName, "", jetRadius,
							      jetPtCut, jetAreaCut, jetBias, biasType, matching, matchingLevel, matchingLevel, "TPC");
    taskCharged->SetHistoType(histoType);
    taskCharged->GetParticleContainer(1)->SelectPhysicalPrimaries(kTRUE);
  }
  
  if (doemcal) {
    AliEmcalJetTask *jetTask = AddTaskEmcalJet(tracksName, corrClusName, 1, jetRadius, 0, trackPtCut, clusPtCut, kGhostArea, 1, "Jet", 1, kFALSE, kFALSE);
    const char *jetsName = jetTask->GetName();
    AliEmcalJetTask *mcJetTask = AddTaskEmcalJet(mcTracksName, "", 1, jetRadius, 0, partLevPtCut, partLevPtCut, kGhostArea, 1, "Jet", 1, kTRUE, kFALSE);
    const char *mcJetsName = mcJetTask->GetName();

    // Analysis task
    AliJetResponseMaker *taskFull = AddTaskJetResponseMaker(tracksName, corrClusName, jetsName, "", jetRadius, mcTracksName, "", mcJetsName, "", jetRadius,
                                                            jetPtCut, jetAreaCut, jetBias, biasType, matching, matchingLevel, matchingLevel, "EMCAL");
    taskFull->SetHistoType(histoType);
    taskFull->SetCaloCellsName(cellsName);
    taskFull->GetParticleContainer(1)->SelectPhysicalPrimaries(kTRUE);
  }

  // HF-jet analysis
  if (dohf) {

    AliAnalysisTaskSEDmesonsFilterCJ* pDMesonFilterRec = AddTaskSEDmesonsFilterCJ(kDmesonType,
										   "",
										   kTRUE,  //   Bool_t theMCon
										   kTRUE,   //   Bool_t reco
										   "");
    pDMesonFilterRec->SetCombineDmesons(kTRUE);
    //pDMesonFilterRec->SetRejectDfromB(kFALSE);
    //pDMesonFilterRec->SetKeepOnlyDfromB(kTRUE);
    AliParticleContainer* trackContDMeson = pDMesonFilterRec->AddParticleContainer(tracksName);
    trackContDMeson->SetClassName("AliAODTrack");
    trackContDMeson->SetFilterHybridTracks(kTRUE);

    if (1) {
      // QA task
      AliAnalysisTaskSAQA *pQADcandidatesTask = AddTaskSAQA(sTracksDcandidatesName, "", "", "", "", jetRadius, 1, 0, trackPtCut, clusPtCut, "TPC");
      pQADcandidatesTask->SetHistoBins(200, 0, 30);
    }
    
    AliEmcalJetTask *pChJetDMesonTask = AddTaskEmcalJet(sTracksDcandidatesName, "", 1, jetRadius, 1, trackPtCut, clusPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, 1);
    TString sChJetsDMesonName = pChJetDMesonTask->GetName();

    AliAnalysisTaskDmesonJetCorrelations* pDMesonJetCorrRec = AddTaskDmesonJetCorr(kDmesonCorrType, "", 
										   sTracksDcandidatesName, "", sChJetsDMesonName, "",
										   jetRadius, jetPtCut, jetAreaCut, "TPC", 0, sDmesonCandName,
										   "AliAnalysisTaskDmesonJetCorrelations", "MCrec");
    pDMesonJetCorrRec->SetMaxR(jetRadius);
    pDMesonJetCorrRec->SetMatchingType(eFlavourJetMatchingType);
    pDMesonJetCorrRec->SetPlotOnlyAcceptedJets(kTRUE);
    pDMesonJetCorrRec->SetShowDeltaEta(kTRUE);
    pDMesonJetCorrRec->SetShowDeltaPhi(kTRUE);
    if (kDmesonCorrType == AliAnalysisTaskDmesonJetCorrelations::kDstartoKpipi) pDMesonJetCorrRec->SetShow2ProngInvMass(kTRUE);
    pDMesonJetCorrRec->SetShowInvMass(kTRUE);
    pDMesonJetCorrRec->SetShowJetConstituents(kTRUE);
    return;
    // MC particle selector    
    AliMCHFParticleSelector *mcPartTask = AddTaskMCHFParticleSelector(mcTracksDMesonname, kFALSE, kTRUE, 0.);
    mcPartTask->SetOnlyPhysPrim(kTRUE);
    if (kDmesonType == AliAnalysisTaskSEDmesonsFilterCJ::kDstartoKpipi) {
      mcPartTask->SelectCharmtoDStartoKpipi();
    }
    else {
      mcPartTask->SelectCharmtoD0toKpi();
    }
    //mcPartTask->SetRejectDfromB(kFALSE);
    //mcPartTask->SetKeepOnlyDfromB(kTRUE);
    //mcPartTask->SelectCharmtoD0toKpi();

    AliEmcalJetTask *chMcJetTaskchDMeson = AddTaskEmcalJet(mcTracksDMesonname, "", 1, jetRadius, 1, partLevPtCut, partLevPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, 0);
    const char *chMcJetsDMesonName = chMcJetTaskchDMeson->GetName();
    
    AliAnalysisTaskSEDmesonsFilterCJ* pDMesonFilterGen = AddTaskSEDmesonsFilterCJ(kDmesonType,
										  "",
										  kTRUE,  //   Bool_t theMCon
										  kFALSE,   //   Bool_t reco
										  "");
    //pDStarMesonFilterTruthTask->SetRejectDfromB(kFALSE);
    //pDStarMesonFilterTruthTask->SetKeepOnlyDfromB(kTRUE);

    AliAnalysisTaskDmesonJetCorrelations* pDMesonJetCorrGen = AddTaskDmesonJetCorr(kDmesonCorrType, "", 
										   mcTracksDMesonname, "", chMcJetsDMesonName, "",
										   jetRadius, jetPtCut, jetAreaCut, "TPC", 0, "Dcandidates",
										   "AliAnalysisTaskDmesonJetCorrelations", "MC");
    pDMesonJetCorrGen->SetMaxR(jetRadius);
    pDMesonJetCorrGen->SetMatchingType(eFlavourJetMatchingType);
    pDMesonJetCorrGen->SetPlotOnlyAcceptedJets(kTRUE);
    pDMesonJetCorrGen->SetShowDeltaEta(kTRUE);
    pDMesonJetCorrGen->SetShowDeltaPhi(kTRUE);
    if (kDmesonCorrType == AliAnalysisTaskDmesonJetCorrelations::kDstartoKpipi) pDMesonJetCorrGen->SetShow2ProngInvMass(kTRUE);
    pDMesonJetCorrGen->SetShowInvMass(kTRUE);
    pDMesonJetCorrGen->SetParticleLevel(kTRUE);
    pDMesonJetCorrGen->SetShowJetConstituents(kTRUE);

    // Analysis task
    AliJetResponseMaker *taskCharged = AddTaskJetResponseMaker(sTracksDcandidatesName, "", sChJetsDMesonName, "", jetRadius, mcTracksDMesonname, "", chMcJetsDMesonName, "", jetRadius,
							      jetPtCut, jetAreaCut, jetBias, biasType, matching, matchingLevel, matchingLevel, "TPC");
    taskCharged->SetHistoType(histoType);
    //taskCharged->GetParticleContainer(1)->SelectPhysicalPrimaries(kTRUE);
    taskCharged->SetFlavourZAxis(kTRUE);
    taskCharged->SetFlavourPtAxis(kTRUE);
    taskCharged->GetJetContainer(0)->SetFlavourCut(kFlavourCut);
    taskCharged->GetJetContainer(1)->SetFlavourCut(kFlavourCut);
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
        taskEmcal->SetNeedEmcalGeom(doemcal);
      }
    }
  }
}
