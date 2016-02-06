void runEMCalJetAnalysisNew_AP(
			       const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
			       const char   *cLocalFiles    = "fileLists/files_LHC12c_pass2_AOD.txt",   // set the local list file
			       UInt_t        iNumFiles      = 100,                                     // number of files analyzed locally
			       UInt_t        iNumEvents     = 10000,                                    // number of events to be analyzed
			       const char   *cRunPeriod     = "LHC12c",                                // set the run period
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

  enum EAnaType {
    kLocal = 0,
    kGrid = 2,
  };

  EAnaType anaType = kLocal;

  TString OCDBpath;
  
  if (anaType == kLocal) {
    OCDBpath = "local:///Volumes/DATA/ALICE/OCDB/2012";
  }
  else {
    OCDBpath = "raw://";
  }
    
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/runEMCalJetAnalysisNew.C");
  AliAnalysisManager* pMgr = runEMCalJetAnalysisNew(cDataType, cLocalFiles, iNumFiles, iNumEvents, cRunPeriod, kPhysSel, cTaskName, kTRUE, kTRUE, OCDBpath, anaType, "terminate");

  if (anaType == kLocal) {
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
