// AddTaskJetSim.C

void AddTaskJetSim()
{
  const char *mcTracksName        = "MCParticles";
  Double_t    jetRadius           = 0.6;
  Float_t     partLevPtCut        = 0.;
  Float_t     jetBias             = 0;
  Int_t       biasType            = 0;   //  0 = charged, 1 = neutral, 2 = both
  Float_t     jetPtCut            = 1;
  Float_t     jetAreaCut          = 0.;
  Int_t       histoType           = 1;  // 0 = TH2, 1 = THnSparse
  Bool_t      forcePP             = kTRUE;
  Double_t    kGhostArea          = 1;
  Int_t eFlavourJetMatchingType   = AliAnalysisTaskDmesonJetCorrelations::kCandidateConstituentMatching;    // kGeometricalMatching, kDaughterConstituentMatching, kCandidateConstituentMatching, kJetLoop

  AliAnalysisTaskDmesonJetCorrelations::ECandidateType kDmesonCorrType = AliAnalysisTaskDmesonJetCorrelations::kD0toKpi; //kDstartoKpipi  
  AliEmcalJet::EFlavourTag kFlavourCut = AliEmcalJet::kD0;
  TString mcTracksDMesonname = "mcparticlesD0";

  UInt_t   matching = 1; //1=geometrical, 2=MClabel
  Double_t matchingLevel = 0.99;
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMCTrackSelector.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskJetResponseMaker.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAQA.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAJF.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskMCHFParticleSelector.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskDmesonJetCorr.C");
  
  // HF-jet analysis
  if (1) {
    // MC particle selector    
    AliMCHFParticleSelector *mcPartTask = AddTaskMCHFParticleSelector(mcTracksDMesonname, kFALSE, kTRUE, 1.);
    mcPartTask->SetOnlyPhysPrim(kTRUE);
    if (kDmesonCorrType == AliAnalysisTaskDmesonJetCorrelations::kDstartoKpipi) {
      mcPartTask->SelectCharmtoDStartoKpipi();
    }
    else {
      mcPartTask->SelectCharmtoD0toKpi();
    }
    //mcPartTask->SetRejectDfromB(kFALSE);
    //mcPartTask->SetKeepOnlyDfromB(kTRUE);
    //mcPartTask->SelectCharmtoD0toKpi();

    AliEmcalJetTask *chMcJetTaskchDMeson = AddTaskEmcalJet(mcTracksDMesonname, "", 1, jetRadius, 0, partLevPtCut, partLevPtCut, kGhostArea, 1, "Jet", 0., kFALSE, kFALSE, 0);
    const char *chMcJetsDMesonName = chMcJetTaskchDMeson->GetName();
    
    AliAnalysisTaskDmesonJetCorrelations* pDMesonJetCorrGen = AddTaskDmesonJetCorr(kDmesonCorrType, "", 
										   mcTracksDMesonname, "", chMcJetsDMesonName, "",
										   jetRadius, jetPtCut, jetAreaCut, "TPC", 0, "", kFALSE,
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
  }

  if (1) {
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      ::Error("AddTaskJetResp", "No analysis manager to connect to.");
      return NULL;
    }  

    TObjArray *toptasks = mgr->GetTasks();
    for (Int_t i = 0; i < toptasks->GetEntries(); ++i) {
      AliAnalysisTaskSE *task = dynamic_cast<AliAnalysisTaskSE*>(toptasks->At(i));
      if (!task) continue;

      if (task->InheritsFrom("AliAnalysisTaskEmcal")) {
	AliAnalysisTaskEmcal *taskEmcal = dynamic_cast<AliAnalysisTaskEmcal*>(task);
	if (forcePP) taskEmcal->SetForceBeamType(AliAnalysisTaskEmcal::kpp);
	taskEmcal->SetVzRange(-10,10);
      }
    }
  }
}
