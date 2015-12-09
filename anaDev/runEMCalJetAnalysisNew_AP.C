void runEMCalJetAnalysisNew_AP(
			       const char   *cDataType      = "AOD",                                   // set the analysis type, AOD or ESD
			       const char   *cLocalFiles    = "fileLists/files_LHC11h_2_AOD145.txt",   // set the local list file
			       UInt_t        iNumFiles      = 100,                                     // number of files analyzed locally
			       UInt_t        iNumEvents     = 5000,                                    // number of events to be analyzed
			       const char   *cRunPeriod     = "LHC11h",                                // set the run period
			       const char   *cTaskName      = "JetAna"
			       )
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/runEMCalJetAnalysisNew.C");
  runEMCalJetAnalysisNew(cDataType, cLocalFiles, iNumFiles, iNumEvents, cRunPeriod, cTaskName);
}
