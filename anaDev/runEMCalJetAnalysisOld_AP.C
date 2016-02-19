void runEMCalJetAnalysisOld_AP(
			       const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
			       const char   *cLocalFiles    = "fileLists/files_LHC13b4_plus_AOD152_5.txt",   // set the local list file
			       UInt_t        iNumFiles      = 100,                                     // number of files analyzed locally
			       UInt_t        iNumEvents     = 20000,                                    // number of events to be analyzed
			       const char   *cRunPeriod     = "LHC13b4",                                // set the run period
			       const char   *cTaskName      = "JetAna"
			       )
{
  // AliEmcalPhysicsSelection::kEmcalOk, AliEmcalPhysicsSelection::kEmcalH,
  // AliVEvent::kINT7, AliVEvent::kMB, AliVEvent::kCentral, AliVEvent::kSemiCentral,
  // AliVEvent::kEMCEGA, AliVEvent::kEMCEJE
  UInt_t kPhysSel = AliVEvent::kAnyINT;
  
  enum eDataType { kAod, kEsd };

  eDataType iDataType;
  if (!strcmp(cDataType, "ESD")) {
    iDataType = kEsd;
  }
  else if (!strcmp(cDataType, "AOD")) {
    iDataType = kAod;
  }
  else {
    Printf("Incorrect data type option, check third argument of run macro.");
    Printf("datatype = AOD or ESD");
    return;
  }

  Bool_t doNotStart = kTRUE;
  TString OCDBpath = "local:///Volumes/DATA/ALICE/OCDB/2012";
    
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/runEMCalJetAnalysisOld.C");
  AliAnalysisManager* pMgr = runEMCalJetAnalysisOld(cDataType, cLocalFiles, iNumFiles, iNumEvents, cRunPeriod, kPhysSel, cTaskName, kTRUE, kFALSE, OCDBpath, doNotStart);

  if (doNotStart) {
    //pMgr->SetUseProgressBar(kFALSE, 10);
    //pMgr->SetDebugLevel(2);
  
    //pMgr->AddClassDebug("AliEmcalJetTask", AliLog::kDebug+100);

    TChain* pChain = 0;
    if (iDataType == kAod) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      pChain = CreateAODChain(cLocalFiles, iNumFiles, 0, kFALSE);
    }
    else {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      pChain = CreateESDChain(cLocalFiles, iNumFiles, 0, kFALSE);
    }

    // start analysis
    Printf("Starting Analysis...");
    pMgr->StartAnalysis("local", pChain, iNumEvents);
  }
}
