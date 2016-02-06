void runEMCalJetAnalysisOld_AP(
			       const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
			       const char   *cLocalFiles    = "fileLists/files_LHC11h_2_AOD145.txt",   // set the local list file
			       UInt_t        iNumFiles      = 100,                                     // number of files analyzed locally
			       UInt_t        iNumEvents     = 100,                                    // number of events to be analyzed
			       const char   *cRunPeriod     = "LHC11h",                                // set the run period
			       const char   *cTaskName      = "JetAna"
			       )
{
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
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/runEMCalJetAnalysisOld.C");
  AliAnalysisManager* pMgr = runEMCalJetAnalysisOld(cDataType, cLocalFiles, iNumFiles, iNumEvents, cRunPeriod, cTaskName, kTRUE, kTRUE, "local:///Volumes/DATA/ALICE/OCDB/2011", kTRUE);

  //pMgr->SetUseProgressBar(kFALSE, 250);
  //pMgr->SetDebugLevel(1);
  
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
