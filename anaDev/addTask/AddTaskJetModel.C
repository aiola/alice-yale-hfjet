// $Id$

void AddTaskJetModel(Bool_t pdsftrain = kTRUE, const char *datatype = "AOD", const char *runtype = "carver",
		     Bool_t doFullJets = kTRUE, Bool_t doChargedJets = kFALSE, Bool_t doEmb = kTRUE, Bool_t doHadCorr = kTRUE, Bool_t doTender = kTRUE)
{
  // avoid automatic setting of physics selection
  //gSystem->Setenv("ETRAIN_DEFPSEL",0);

  TString trainpath(gSystem->Getenv("ETRAIN_PATH"));

  enum eDataType { kAod, kEsd, kSesd };
  enum eRunType  { kLocal, kGrid, kPdsf, kCarver };

  eRunType rType;
  if (!strcmp(runtype, "grid")) 
    rType = kGrid;
  else if (!strcmp(runtype, "local")) 
    rType = kLocal;
  else if (!strcmp(runtype, "pdsf")) 
    rType = kPdsf;
  else if (!strcmp(runtype, "carver")) 
    rType = kCarver;
  else {
    cout << "Incorrect run option, check first argument of run macro" << endl;
    cout << "runtype = local, grid, pdsf, carver" << endl;
    return;
  }

  TString period(gSystem->Getenv("ETRAIN_DATASET"));
  period.ToLower();
  if (pdsftrain) {
    if (!period.BeginsWith("lhc12a15e") && !period.BeginsWith("lhc12a15a")) {
      cout << "AddTaskJetModel.C ignored because dataset is not lhc12a15e nor lhc12a15a." << endl;
      return;
    }
  }

  eDataType dType;
  if (!strcmp(datatype, "ESD"))
    dType = kEsd;
  else if (!strcmp(datatype, "sESD"))
    dType = kSesd;
  else if (!strcmp(datatype, "AOD"))
    dType = kAod;
  else {
    cout << "Incorrect data type option, check third argument of run macro" << endl;
    cout << "datatype = AOD, ESD or sESD (skimmed ESD)" << endl;
    return;
  }

  Double_t minJetPt = 0;
  Int_t minMClabel = 0;

  TString tracksName("PicoTracks");
  TString mcTracksName("MCSelectedParticles");

  TString clusName(dType == kAod ? "caloClusters" : "CaloClusters");
  TString corrClusName;

  Double_t jetRadius     = 0.2;
  Double_t trackPtCut    = 0.15;
  Double_t clusPtCut     = 0.30;
  Double_t partLevPtCut  = 0.;
  Double_t jetAreaCut    = 0.;
  Double_t jetPtCut      = 1;
  Double_t ghostArea     = 0.005;
  Double_t hadCorrF      = 2;
  Int_t    histoType     = 1;
  Double_t propDist      = 440.;

  TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
  
  if (hadCorrF > 1.9) {
    // 300 MeV, had corr = 2
    sfunc->SetParameter(2,1.8078);
    sfunc->SetParameter(1,-0.0115442);
    sfunc->SetParameter(0,0.000109722);
  }
  else if (hadCorrF > 1.6) {
    // 300 MeV, had corr = 1.7
    sfunc->SetParameter(2, 1.79148);
    sfunc->SetParameter(1, -1.13795e-02);
    sfunc->SetParameter(0, 1.14528e-04);
  }
  else if (hadCorrF > 1.2) {
    // 300 MeV, had corr = 1.3
    sfunc->SetParameter(2,1.90655);
    sfunc->SetParameter(1,-0.0159274);
    sfunc->SetParameter(0,0.0001706891);
  }
  else {
    // 300 MeV, had corr = 0
    sfunc->SetParameter(2,2.01265);
    sfunc->SetParameter(1,-0.0128864);
    sfunc->SetParameter(0,0.000100894);
  }

  UInt_t   matching = 1; //1=geometrical, 2=MClabel, 3=SameCollections
  Double_t matchingLevel1 = 0.99;
  Double_t matchingLevel2 = 0.99;

  Float_t kJetLeadingTrackBias = 5;
  Int_t kNcent = 0;
  
  Int_t kNjetResp = 1;
  
  if (kJetLeadingTrackBias > 1)
    kNjetResp = 2;
  
  if (kJetLeadingTrackBias > 5)
    kNjetResp = 3;
  
  kNjetResp *= kNcent!=0 ? kNcent : 1;

  UInt_t physSel = 0;
  //UInt_t physSel = AliVEvent::kAny;
  //UInt_t physSel = AliVEvent::kINT7;
  //UInt_t physSel = AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  //UInt_t physSel = AliVEvent::kEMCEGA;
  //UInt_t physSel = AliVEvent::kEMCEJE;
  //UInt_t physSel = AliEmcalPhysicsSelection::kEmcalHT;
  //UInt_t physSel = AliEmcalPhysicsSelection::kEmcalOk;

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskJetModel", "No analysis manager to connect to.");
    return;
  }  

  // Loading macros
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskSAQA.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskRho.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskDeltaPt.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskSAJF.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetEmbedding.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskJetRespPtHard.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalAodTrackFilter.C");

  if (dType == kEsd || dType == kAod) {
    const char *inputTracks = "HybridTracks";

    if (dType == kEsd) {
      // Hybrid tracks maker for ESD
      AliEmcalEsdTpcTrackTask *hybTask = AddTaskEmcalEsdTpcTrack(inputTracks, "Hybrid_LHC11h", kTRUE);
      hybTask->SelectCollisionCandidates(physSel);
      
      // Track propagator
      AliEmcalTrackPropagatorTask *propTask = AddTaskEmcalTrackPropagator(inputTracks, propDist);
      propTask->SelectCollisionCandidates(physSel);
    }

    if (dType == kAod) {
      // Hybrid tracks maker for AOD
      AddTaskEmcalAodTrackFilter *hybTask = AddTaskEmcalAodTrackFilter(inputTracks, "tracks", "lhc11h");
      hybTask->SelectCollisionCandidates(physSel);
    }

    // PicoTracks maker
    AliEmcalPicoTrackMaker *pTrackTask = AddTaskEmcalPicoTrackMaker(tracksName, inputTracks);
    pTrackTask->SelectCollisionCandidates(physSel);
    pTrackTask->SetTrackEfficiency(0.95);
    pTrackTask->SetUseNegativeLabels(kTRUE);
  }
  
  if (1) {
    // MC particle selector
    mcTracksName = "MCSelectedParticles";
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
    AliEmcalMCTrackSelector *mcPartTask = AddTaskMCTrackSelector(mcTracksName, kFALSE, kFALSE);
    mcPartTask->SelectCollisionCandidates(physSel);
  }

  // Tender Supplies
  if (doTender) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEMCALTender.C");
    AliAnalysisTaskSE *tender = AddTaskEMCALTender(kTRUE, kTRUE, kTRUE, kTRUE, kTRUE, kFALSE, kTRUE, kTRUE, kTRUE, AliEMCALRecoUtils::kBeamTestCorrected,
						   kFALSE, 0.1, 0.05, AliEMCALRecParam::kClusterizerv2, kFALSE, kFALSE, -1, 1e6, 1e6, "pass2");
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
    tender->SelectCollisionCandidates(physSel);
    if (pdsftrain)
      RequestMemory(tender, 250 * 1024);
  }
  
  if (doHadCorr) {
    corrClusName = "CaloClustersCorr";
    
    TString emcalTracksName("EmcalTracks");
    TString emcalClusName("EmcalClusters");
    
    if (!mgr->GetTask(Form("HadCorr_%s",corrClusName.Data()))) {
      // EmcalParticles maker
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalParticleMaker.C");
      AliEmcalParticleMaker *ePartTask = AddTaskEmcalParticleMaker(tracksName, clusName, emcalTracksName, emcalClusName, "AliEmcalParticleMaker");
      ePartTask->SelectCollisionCandidates(physSel);
      
      // Cluster-track matcher task
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C");
      AliEmcalClusTrackMatcherTask *matcherTask = AddTaskEmcalClusTrackMatcher(emcalTracksName, emcalClusName, 0.1, kFALSE);
      matcherTask->SelectCollisionCandidates(physSel);
      
      // Hadronic correction task
      gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskHadCorr.C");
      AliHadCorrTask *hCorrTask = AddTaskHadCorr(emcalTracksName, emcalClusName, corrClusName, 
						 hadCorrF, 0.15, 0.030, 0.015, 0, kTRUE, kTRUE);
      hCorrTask->SelectCollisionCandidates(physSel);
      hCorrTask->SetHistoBins(100,0,100);
      if (pdsftrain)
	RequestMemory(hCorrTask,400*1024);
    }
  }
  else {
    corrClusName = clusName;
  }
  
  if (doEmb) {
    TF1 *boltzman = new TF1("boltzman","exp(-x/[0])",0,250);
    boltzman->SetNpx(1000);
    boltzman->SetParameter(0,0.550);

    TF1 *density = new TF1("density","x<[0] ? 0 : (x>[1] ? 0 : 1)",0,6000);
    density->SetNpx(24000);
    density->SetParameter(0,3100);
    density->SetParameter(1,5150);
    
    AliJetEmbeddingTask *embTask = AddTaskJetEmbedding(tracksName,"","JetEmbeddingTask",0,2);
    //embTask->SetNTracks(0);
    embTask->SetMarkMC(0);
    embTask->SetDensitySpectrum(density);
    embTask->SetPtSpectrum(boltzman);
    embTask->SelectCollisionCandidates(physSel);
    if (pdsftrain)
      RequestMemory(embTask, 1000 * 1024);
  }

  if (1) {
    AliAnalysisTaskSAQA *mcqatask = AddTaskSAQA(mcTracksName, "", "", "", "", jetRadius, jetPtCut, 0, 0, 0, "TPC");
    mcqatask->SetHistoBins(200,0,100);
    mcqatask->SetParticleLevel(kTRUE);
    mcqatask->SetMC(kTRUE);
    mcqatask->SelectCollisionCandidates(physSel);
    if (pdsftrain)
      RequestMemory(mcqatask, 100 * 1024);
  }

  if (0) {
    // QA task
    AliAnalysisTaskSAQA *qaTask = AddTaskSAQA(tracksName, corrClusName, "", "", "",
					      0, 0, 0, trackPtCut, clusPtCut, "TPC", "AliAnalysisTaskSAQA");
    qaTask->SetHistoBins(100, 0, 100);
    qaTask->SelectCollisionCandidates(physSel);
    qaTask->SetMC(kTRUE);
    qaTask->SetIsEmbedded(doEmb);
    if (pdsftrain)
      RequestMemory(qaTask, 400 * 1024);
  }

  gROOT->LoadMacro("$ETRAIN_ROOT/saiola/AddTaskMCClustersQA.C");

  if (0) {
    AliMCClustersQATask *qaTask = AddTaskMCClustersQA(corrClusName, mcTracksName);
    qaTask->SelectCollisionCandidates(physSel);
    qaTask->SetIsPythia(kTRUE);
    if (pdsftrain)
      RequestMemory(qaTask, 800 * 1024);
  }

  if (0) {
    AliMCClustersQATask *qaTask = AddTaskMCClustersQA(clusName, mcTracksName);
    qaTask->SelectCollisionCandidates(physSel);
    qaTask->SetMakeParticleHisto(kFALSE);
    if (pdsftrain)
      RequestMemory(qaTask, 400 * 1024);
  }

  if (doEmb) {
    if (1) {
      // QA task
      AliAnalysisTaskSAQA *qaTask = AddTaskSAQA(tracksName, "", "", "", "",
						0, 0, 0, trackPtCut, 0, "TPC", "AliAnalysisTaskSAQA_Background");
      qaTask->GetParticleContainer(0)->SetMCTrackBitMap(TObject::kBitMask);
      qaTask->GetParticleContainer(0)->SetTrackBitMap(0);
      qaTask->SetHistoBins(200, 0, 100);
      qaTask->SelectCollisionCandidates(physSel);
      if (pdsftrain)
	RequestMemory(qaTask, 400 * 1024);
    }

    if (1) {
      // QA task
      AliAnalysisTaskSAQA *qaTask = AddTaskSAQA(tracksName, corrClusName, "", "", "",
						0, 0, 0, trackPtCut, clusPtCut, "TPC", "AliAnalysisTaskSAQA_Signal");
      qaTask->GetParticleContainer(0)->SetMCTrackBitMap(0);
      qaTask->GetParticleContainer(0)->SetTrackBitMap(TObject::kBitMask);
      qaTask->GetClusterContainer(0)->SetMCClusterBitMap(0);
      qaTask->GetClusterContainer(0)->SetClusterBitMap(TObject::kBitMask);
      qaTask->SetMC(kTRUE);
      qaTask->SetHistoBins(200, 0, 100);
      qaTask->SelectCollisionCandidates(physSel);
      if (pdsftrain)
	RequestMemory(qaTask, 400 * 1024);
    }
  }

  TString chRhoName;
  TString chKtJetsName;

  if (doEmb) {
    AliEmcalJetTask *chKtJetTask = AddTaskEmcalJet(tracksName, "", 0, 0.2, 1, trackPtCut, clusPtCut, ghostArea);
    chKtJetTask->SetMinMCLabel(minMClabel);
    //chKtJetTask->SelectConstituents(0, TObject::kBitMask);
    chKtJetTask->SelectCollisionCandidates(physSel);
    if (pdsftrain)
      RequestMemory(chKtJetTask, 500 * 1024);
    
    chKtJetsName = chKtJetTask->GetName();
    chRhoName = "TpcRho_ExLJ";

    // Rho task charged jet
    AliAnalysisTaskRho *chRhoTask = AddTaskRho(chKtJetsName, tracksName, "", chRhoName, 
					       0.2, "TPC", 0.01, 0, 0, 2, kTRUE, "Rho_ExLJ");
    chRhoTask->GetParticleContainer(0)->SetMinMCLabel(minMClabel);
    chRhoTask->GetParticleContainer(0)->SetMCTrackBitMap(TObject::kBitMask);
    chRhoTask->SetHistoBins(250,0,250);
    chRhoTask->SelectCollisionCandidates(physSel);

    if (1) {
      AliAnalysisTaskDeltaPt *DeltaPtChTask  = AddTaskDeltaPt(tracksName, "", "", 
							      "", "", "", 
							      "", "", chRhoName, 
							      jetRadius, jetAreaCut, 
							      trackPtCut, clusPtCut, "TPC","AliAnalysisTaskDeltaPt_Background");
      DeltaPtChTask->GetParticleContainer(0)->SetMinMCLabel(minMClabel);
      DeltaPtChTask->GetParticleContainer(0)->SetMCTrackBitMap(TObject::kBitMask);
      DeltaPtChTask->SetMCJetPtThreshold(1);
      DeltaPtChTask->SelectCollisionCandidates(physSel);
      if (pdsftrain)
	RequestMemory(DeltaPtChTask, 100 * 1024);
    }
  }

  if (doChargedJets) {
    AliEmcalJetTask *chDetLevJetTask = AddTaskEmcalJet(tracksName, "", 1, jetRadius, 1, trackPtCut, clusPtCut, ghostArea);
    chDetLevJetTask->SetMinMCLabel(minMClabel);
    chDetLevJetTask->SelectCollisionCandidates(physSel);
    if (pdsftrain)
      RequestMemory(chDetLevJetTask, 200 * 1024);

    const char *chDetLevJetsName = chDetLevJetTask->GetName();
    
    if (doEmb) {
      if (1) {
	AliAnalysisTaskSAJF *SpectraChTask = AddTaskSAJF(tracksName, "", chDetLevJetsName, chRhoName, 
							 jetRadius, jetPtCut, jetAreaCut, 
							 "TPC");
	SpectraChTask->SetMinMCLabel(minMClabel);
	SpectraChTask->SetNLeadingJets(1);
	SpectraChTask->SetIsEmbedded(doEmb);
	SpectraChTask->SetHistoBins(400,0,200);
	SpectraChTask->SelectCollisionCandidates(physSel);
	if (pdsftrain)
	  RequestMemory(SpectraChTask, 400 * 1024);
      }
      
      if (0) {
	AliAnalysisTaskDeltaPt *DeltaPtChTask  = AddTaskDeltaPt(tracksName, "", "", 
								"", "", "", 
								"", "", chRhoName, 
								jetRadius, jetAreaCut, 
								trackPtCut, clusPtCut, "TPC", "AliAnalysisTaskDeltaPt");
	DeltaPtChTask->SetMinMCLabel(minMClabel);
	DeltaPtChTask->SetMCJetPtThreshold(1);
	DeltaPtChTask->SelectCollisionCandidates(physSel);
	if (pdsftrain)
	  RequestMemory(DeltaPtChTask, 100 * 1024);
      }
    }

    AliEmcalJetTask *chPartLevJetTask = AddTaskEmcalJet(mcTracksName, "", 1, jetRadius, 1, partLevPtCut, partLevPtCut, ghostArea);
    chPartLevJetTask->SelectCollisionCandidates(physSel);
    if (pdsftrain)
      RequestMemory(chPartLevJetTask, 200 * 1024);

    const char *chPartLevJetsName       = chPartLevJetTask->GetName();

    if (0) {
      AliAnalysisTaskSAJF *SpectraPartLevChTask = AddTaskSAJF(mcTracksName, "", chPartLevJetsName, "", 
							      jetRadius, jetPtCut, jetAreaCut, 
							      "TPC"); 
      SpectraPartLevChTask->SetNLeadingJets(1);
      SpectraPartLevChTask->SetMaxClusterPt(1000);
      SpectraPartLevChTask->SetMaxTrackPt(1000);
      SpectraPartLevChTask->SelectCollisionCandidates(physSel);
      if (pdsftrain)
	RequestMemory(SpectraPartLevChTask, 250 * 1024);
    }

    if (1) { 
      AliJetResponseMaker *rmTaskCh = AddTaskJetRespPtHard(tracksName, "", chDetLevJetsName, chRhoName, jetRadius,
							   mcTracksName, "", chPartLevJetsName, "", jetRadius,
							   jetPtCut, jetAreaCut, kJetLeadingTrackBias, 0, matching, matchingLevel1, matchingLevel1, 
							   "TPC", -999, -999, kNcent, "AliJetResponseMaker");
      for (Int_t i = 0; i < kNjetResp; i++) {
	rmTaskCh[i].SelectCollisionCandidates(physSel);
	rmTaskCh[i].SetHistoType(histoType);
	rmTaskCh[i].SetIsParticleLevel(kFALSE,0);
	rmTaskCh[i].SetIsParticleLevel(kTRUE,1);
	rmTaskCh[i].SetIsEmbedded(doEmb);
	rmTaskCh[i].SetMinMCLabel(minMClabel);
	rmTaskCh[i].SetMCLabelShift(minMClabel);
	rmTaskCh[i].SetHistoBins(500,0,250);
	if (pdsftrain)
	  RequestMemory(&rmTaskCh[i],500*1024);
      }
    }

    if (0) { // MC label
      AliJetResponseMaker *rmTaskCh = AddTaskJetRespPtHard(tracksName, "", chDetLevJetsName, chRhoName, jetRadius,
							   mcTracksName, "", chPartLevJetsName, "", jetRadius,
							   jetPtCut, jetAreaCut, 0, kJetLeadingTrackBias, 2, 0.99, 0.99, 
							   "TPC", -999, -999, kNcent, "AliJetResponseMaker_MCLabel");
      for (Int_t i = 0; i < kNjetResp; i++) {
	rmTaskCh[i].SelectCollisionCandidates(physSel);
	rmTaskCh[i].SetHistoType(histoType);
	rmTaskCh[i].SetIsParticleLevel(kFALSE,0);
	rmTaskCh[i].SetIsParticleLevel(kTRUE,1);
	rmTaskCh[i].SetIsEmbedded(doEmb);
	rmTaskCh[i].SetMinMCLabel(minMClabel);
	rmTaskCh[i].SetMCLabelShift(minMClabel);
	rmTaskCh[i].SetHistoBins(500,0,250);
	rmTaskCh[i].SetNEFAxis(0);
	rmTaskCh[i].SetZAxis(1);
	if (pdsftrain)
	  RequestMemory(&rmTaskCh[i],500*1024);
      }
    }
    
    if (doEmb) {
      AliEmcalJetTask *chDetLevMConlyJetTask = AddTaskEmcalJet(tracksName, "", 1, jetRadius, 1, trackPtCut, clusPtCut, ghostArea, "JetMConly");
      chDetLevMConlyJetTask->SetMinMCLabel(minMClabel);
      chDetLevMConlyJetTask->SelectConstituents(TObject::kBitMask, 0);
      chDetLevMConlyJetTask->SelectCollisionCandidates(physSel);
      if (pdsftrain)
	RequestMemory(chDetLevMConlyJetTask, 200 * 1024);

      const char *chDetLevMConlyJetsName  = chDetLevMConlyJetTask->GetName();
  
      if (0) {
	AliAnalysisTaskSAJF *SpectraChMConlyTask = AddTaskSAJF(tracksName, "", chDetLevMConlyJetsName, chRhoName, 
							       jetRadius, jetPtCut, jetAreaCut, 
							       "TPC");
	SpectraChTask->SetMinMCLabel(minMClabel);
	SpectraChMConlyTask->SetNLeadingJets(1);
	SpectraChMConlyTask->SetTrackBitMap(TObject::kBitMask);
	SpectraChMConlyTask->SelectCollisionCandidates(physSel);
	if (pdsftrain)
	  RequestMemory(SpectraChMConlyTask, 250 * 1024);
      }
      
      if (1) {
	AliJetResponseMaker *rmTaskCh = AddTaskJetRespPtHard(tracksName, "", chDetLevMConlyJetsName, "", jetRadius,
							     mcTracksName, "", chPartLevJetsName, "", jetRadius,
							     jetPtCut, jetAreaCut, kJetLeadingTrackBias, 0, matching, matchingLevel1, matchingLevel2, 
							     "TPC", -999, -999, kNcent);
	for (Int_t i = 0; i < kNjetResp; i++) {
	rmTaskCh[i].SelectCollisionCandidates(physSel);
	rmTaskCh[i].SetHistoType(histoType);
	rmTaskCh[i].SetIsParticleLevel(kFALSE,0);
	rmTaskCh[i].SetIsParticleLevel(kTRUE,1);
	rmTaskCh[i].SetIsEmbedded(kFALSE);
	rmTaskCh[i].SetMinMCLabel(minMClabel);
	rmTaskCh[i].SetMCLabelShift(minMClabel);
	rmTaskCh[i].SetHistoBins(500,0,250);
	  if (pdsftrain)
	    RequestMemory(&rmTaskCh[i],200*1024);
	}
      }

      if (1) {
	AliJetResponseMaker *rmTaskCh = AddTaskJetRespPtHard(tracksName, "", chDetLevJetsName, chRhoName, jetRadius,
							     tracksName, "", chDetLevMConlyJetsName, "", jetRadius,
							     jetPtCut, jetAreaCut, kJetLeadingTrackBias, 0, matching, matchingLevel1, matchingLevel2, 
							     "TPC", -999, -999, kNcent);
	for (Int_t i = 0; i < kNjetResp; i++) {
	  rmTaskCh[i].SelectCollisionCandidates(physSel);
	  rmTaskCh[i].SetHistoType(histoType);
	  rmTaskCh[i].SetIsParticleLevel(kFALSE,0);
	  rmTaskCh[i].SetIsParticleLevel(kFALSE,1);
	  rmTaskCh[i].SetIsEmbedded(kTRUE);
	  rmTaskCh[i].SetMinMCLabel(minMClabel);
	  rmTaskCh[i].SetMCLabelShift(minMClabel);
	  rmTaskCh[i].SetHistoBins(500,0,250);
	  if (pdsftrain)
	    RequestMemory(&rmTaskCh[i],200*1024);
	}
      }
    }
  }

  // EMCal analysis
  if (doFullJets) {
    Double_t kEmcalMinEta = -0.7;
    Double_t kEmcalMaxEta = 0.7;
    Double_t kEmcalMinPhi = 80*TMath::DegToRad();
    Double_t kEmcalMaxPhi = 180*TMath::DegToRad();

    TString fuRhoName = chRhoName;
    
    if (0) {
      gROOT->LoadMacro("$ALICE_ROOT/PWGJE/EMCALJetTasks/macros/AddTaskScale.C");
      AliAnalysisTaskScale *scaleTask = AddTaskScale(tracksName, corrClusName, trackPtCut, clusPtCut);
      scaleTask->SetScaleFunction(sfunc);
      scaleTask->SetMinMCLabel(minMClabel);
      scaleTask->SelectCollisionCandidates(physSel);

      Double_t kRhoMinEta     = -0.5;
      Double_t kRhoMaxEta     = 0.5;
      Double_t kRhoMinPhi     = 30*TMath::DegToRad()+0.2;
      Double_t kRhoMaxPhi     = 230*TMath::DegToRad()-0.2;

      fuRhoName = "TpcRho_Small_Scaled";
      AliAnalysisTaskRho *rhoTaskSmall = AddTaskRho(chKtJetsName, tracksName, corrClusName, "TpcRho_Small", 
						    0.2, user, 0.01, 0, 0, 0, kFALSE, "Rho_Small");
      rhoTaskSmall->SetMinMCLabel(minMClabel);
      rhoTaskSmall->SetScaleFunction(sfunc);
      rhoTaskSmall->SetJetEtaLimits(kRhoMinEta,kRhoMaxEta);
      rhoTaskSmall->SetJetPhiLimits(kRhoMinPhi,kRhoMaxPhi);
      rhoTaskSmall->SetHistoBins(50,0,250);
      rhoTaskSmall->SelectCollisionCandidates(physSel);
    }

    AliEmcalJetTask *fuDetLevJetTask = AddTaskEmcalJet(tracksName, corrClusName, 1, jetRadius, 0, trackPtCut, clusPtCut, ghostArea);
    fuDetLevJetTask->SetMinMCLabel(minMClabel);
    fuDetLevJetTask->SelectCollisionCandidates(physSel);
    if (pdsftrain)
      RequestMemory(fuDetLevJetTask, 200 * 1024);

    const char *fuDetLevJetsName = fuDetLevJetTask->GetName();

    if (doEmb) {
      if (1) {
	AliAnalysisTaskSAJF *SpectraFuTask = AddTaskSAJF(tracksName, corrClusName, fuDetLevJetsName, fuRhoName, 
							 jetRadius, jetPtCut, jetAreaCut, 
							 "EMCAL"); 
	SpectraFuTask->SetMinMCLabel(minMClabel);
	SpectraFuTask->SetNLeadingJets(1);
	SpectraFuTask->SetIsEmbedded(doEmb);
	SpectraFuTask->SelectCollisionCandidates(physSel);
	SpectraFuTask->SetHistoBins(400,0,200);
	if (pdsftrain)
	  RequestMemory(SpectraFuTask, 400 * 1024);
      }

      if (1) {
	AliAnalysisTaskDeltaPt *DeltaPtFuTask = AddTaskDeltaPt(tracksName, corrClusName, "", 
							       "", "", "", 
							       "", "", fuRhoName, 
							       jetRadius, jetAreaCut, 
							       trackPtCut, clusPtCut, "EMCAL", "AliAnalysisTaskDeltaPt"); 
	DeltaPtFuTask->SetMinMCLabel(minMClabel);
	DeltaPtFuTask->SelectCollisionCandidates(physSel);
      }
    }

    AliEmcalJetTask *fuPartLevJetTask = AddTaskEmcalJet(mcTracksName, "", 1, jetRadius, 0, partLevPtCut, partLevPtCut);
    fuPartLevJetTask->SetMinMCLabel(minMClabel);
    fuPartLevJetTask->SelectCollisionCandidates(physSel);
    if (pdsftrain)
      RequestMemory(fuPartLevJetTask, 200 * 1024);

    const char *fuPartLevJetsName        = fuPartLevJetTask->GetName();

    if (0) {
      AliAnalysisTaskSAJF *SpectraMCFuTask = AddTaskSAJF(mcTracksName, "", fuPartLevJetsName, "", 
							 jetRadius, jetPtCut, jetAreaCut, 
							 "EMCAL");
      SpectraMCFuTask->SetNLeadingJets(1);
      SpectraMCFuTask->SetMaxClusterPt(1000);
      SpectraMCFuTask->SetMaxTrackPt(1000);
      SpectraMCFuTask->SelectCollisionCandidates(physSel);
    }

    if (1) {
      AliJetResponseMaker *rmTaskFu = AddTaskJetRespPtHard(tracksName, corrClusName, fuDetLevJetsName, fuRhoName, jetRadius,
							   mcTracksName, "", fuPartLevJetsName, "", jetRadius,
							   jetPtCut, jetAreaCut, kJetLeadingTrackBias, 0, matching, matchingLevel1, matchingLevel2, 
							   "EMCAL", -999, -999, kNcent, "AliJetResponseMaker");
      for (Int_t i = 0; i < kNjetResp; i++) {
	rmTaskFu[i].SelectCollisionCandidates(physSel);
	rmTaskFu[i].SetHistoType(histoType);
	rmTaskFu[i].SetIsParticleLevel(kFALSE,0);
	rmTaskFu[i].SetIsParticleLevel(kTRUE,1);
	rmTaskFu[i].SetIsEmbedded(doEmb);
	rmTaskFu[i].SetMinMCLabel(minMClabel);
	rmTaskFu[i].SetMCLabelShift(minMClabel);
	rmTaskFu[i].SetHistoBins(500,0,250);
	if (pdsftrain)
	  RequestMemory(&rmTaskFu[i],500*1024);
      }
    }

    if (0) { // MC label
      AliJetResponseMaker *rmTaskFu = AddTaskJetRespPtHard(tracksName, corrClusName, fuDetLevJetsName, fuRhoName, jetRadius,
							   mcTracksName, "", fuPartLevJetsName, "", jetRadius,
							   jetPtCut, jetAreaCut, kJetLeadingTrackBias, 0, 2, 0.99, 0.99, 
							   "EMCAL", -999, -999, kNcent, "AliJetResponseMaker_MCLabel");
      for (Int_t i = 0; i < kNjetResp; i++) {
	rmTaskFu[i].SelectCollisionCandidates(physSel);
	rmTaskFu[i].SetHistoType(histoType);
	rmTaskFu[i].SetIsParticleLevel(kFALSE,0);
	rmTaskFu[i].SetIsParticleLevel(kTRUE,1);
	rmTaskFu[i].SetIsEmbedded(doEmb);
	rmTaskFu[i].SetMinMCLabel(minMClabel);
	rmTaskFu[i].SetMCLabelShift(minMClabel);
	rmTaskFu[i].SetHistoBins(500,0,250);
	if (pdsftrain)
	  RequestMemory(&rmTaskFu[i],500*1024);
      }
    }
    
    if (doEmb) {
      AliEmcalJetTask *fuDetLevMConlyJetTask = AddTaskEmcalJet(tracksName, corrClusName, 1, jetRadius, 0, trackPtCut, clusPtCut, ghostArea, "JetMConly");
      fuDetLevMConlyJetTask->SetMinMCLabel(minMClabel);
      fuDetLevMConlyJetTask->SelectConstituents(TObject::kBitMask, 0);
      fuDetLevMConlyJetTask->SelectCollisionCandidates(physSel);
      if (pdsftrain)
	RequestMemory(fuDetLevMConlyJetTask, 200 * 1024);

      const char *fuDetLevMConlyJetsName   = fuDetLevMConlyJetTask->GetName();

      if (0) {
	AliAnalysisTaskSAJF *SpectraMConlyFuTask = AddTaskSAJF(tracksName, corrClusMCOnlyName, fuDetLevMConlyJetsName, "", 
							       jetRadius, jetPtCut, jetAreaCut, 
							       "EMCAL"); 
	SpectraMConlyFuTask->SetMinMCLabel(minMClabel);
	SpectraMConlyFuTask->SetNLeadingJets(1);
	SpectraMConlyFuTask->SetParticleBitMap(TObject::kBitMask);
	SpectraMConlyFuTask->SelectCollisionCandidates(physSel);
      }
	
      if (1) {
	AliJetResponseMaker *rmTaskFu = AddTaskJetRespPtHard(tracksName, corrClusName, fuDetLevMConlyJetsName, "", jetRadius,
							     mcTracksName, "", fuPartLevJetsName, "", jetRadius,
							     jetPtCut, jetAreaCut, kJetLeadingTrackBias, 0, matching, matchingLevel1, matchingLevel2, 
							     "EMCAL", -999, -999, kNcent);
	for (Int_t i = 0; i < kNjetResp; i++) {
	  rmTaskFu[i].SelectCollisionCandidates(physSel);
	  rmTaskFu[i].SetHistoType(histoType);
	  rmTaskFu[i].SetIsParticleLevel(kFALSE,0);
	  rmTaskFu[i].SetIsParticleLevel(kTRUE,1);
	  rmTaskFu[i].SetIsEmbedded(kFALSE);
	  rmTaskFu[i].SetMinMCLabel(minMClabel);
	  rmTaskFu[i].SetMCLabelShift(minMClabel);
	  rmTaskFu[i].SetHistoBins(500,0,250);
	  if (pdsftrain)
	    RequestMemory(&rmTaskFu[i],150*1024);
	}
      }

      if (1) {
	AliJetResponseMaker *rmTaskFu = AddTaskJetRespPtHard(tracksName, corrClusName, fuDetLevJetsName, fuRhoName, jetRadius,
							     tracksName, corrClusName, fuDetLevMConlyJetsName, "", jetRadius,
							     jetPtCut, jetAreaCut, kJetLeadingTrackBias, 0, matching, matchingLevel1, matchingLevel2, 
							     "EMCAL", -999, -999, kNcent);
	for (Int_t i = 0; i < kNjetResp; i++) {
	  rmTaskFu[i].SelectCollisionCandidates(physSel);
	  rmTaskFu[i].SetHistoType(histoType);
	  rmTaskFu[i].SetIsParticleLevel(kFALSE,0);
	  rmTaskFu[i].SetIsParticleLevel(kFALSE,1);
	  rmTaskFu[i].SetIsEmbedded(kTRUE);
	  rmTaskFu[i].SetMinMCLabel(minMClabel);
	  rmTaskFu[i].SetMCLabelShift(minMClabel);
	  rmTaskFu[i].SetHistoBins(500,0,250);
	  if (pdsftrain)
	    RequestMemory(&rmTaskFu[i],150*1024);
	}
      }

    }
  }

  if (1) {
    TObjArray *toptasks = mgr->GetTasks();
    for (Int_t i = 0; i < toptasks->GetEntries(); ++i) {
      AliAnalysisTaskEmcal *taskEmcal = dynamic_cast<AliAnalysisTaskEmcal*>(toptasks->At(i));
      if (taskEmcal && taskEmcal->InheritsFrom("AliAnalysisTaskEmcal")) {
	taskEmcal->SetForceBeamType(AliAnalysisTaskEmcal::kpp);
	taskEmcal->SetVzRange(-10,10);
	taskEmcal->SetIsPythia(kTRUE);
      }
    }
  }
}
