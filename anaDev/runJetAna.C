// runJetAna.C

class AliAnalysisGrid;
class AliESDInputHandler;
class AliVEvent;
class AliEmcalPhysicsSelectionTask;
class AliAnalysisManager;
class AliAODInputHandler;
class AliAODOutputHandler;
class AliEmcalCompatTask;
class AliPicoTrackFixer;
class AliCentralitySelectionTask;
class AliEmcalSetupTask;

void LoadLibs();
void LoadMacros();
AliAnalysisGrid* CreateAlienHandler(const char *cTaskname, const char *cGridMode, const TArrayI &oRunNumbers, 
				    const char *cRunPeriod = "LHC110b", const char *cDataPattern = "*ESDs/pass2/AOD137/*AOD.root");

//______________________________________________________________________________
void runJetAna( 
	 const char   *cDataType     = "AOD",                                 // set the analysis type, AOD or ESD
         const char   *cRunType      = "local",                               // local or grid
         const char   *cGridMode     = "test",                                // set the grid run mode (can be "full", "test", "offline", "submit" or "terminate")
	 const char   *cLocalFiles   = "fileLists/files_LHC10b_AOD137.txt",   // set the local list file
	 UInt_t        iNumFiles     = 100,                                   // number of files analyzed locally
	 UInt_t        iNumEvents    = 50000,                                 // number of events to be analyzed
	 const char   *cRunPeriod    = "LHC110b",                             // set the run period
         const char   *cTaskName     = "JetAna",                              // sets name of grid generated macros
         Bool_t       *bDoEmcal      = kFALSE,
         Bool_t       *bDoHF         = kTRUE
         )
{
  //gSystem->SetFPEMask(TSystem::kInvalid | TSystem::kDivByZero | TSystem::kOverflow | TSystem::kUnderflow);

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
    Printf("runtype = local, grid, pdsf, carver");
    return;
  }

  eDataType iDataType;
  if (!strcmp(cDataType, "ESD"))
    iDataType = kEsd;
  else if (!strcmp(cDataType, "AOD"))
    iDataType = kAod;
  else {
    Printf("Incorrect data type option, check third argument of run macro.");
    Printf("datatype = AOD or ESD");
    return;
  }

  Printf("%s %s analysis chosen.", cRunType, cDataType);

  TString sLocalFiles(cLocalFiles);
  if (iRunType == kLocal || iRunType == kPdsf || iRunType == kCarver) {
    if (sLocalFiles == "") {
      Printf("You need to provide the list of local files!");
      return;
    }
    Printf("Setting local analysis for %d files from list %s, max events = %d", iNumFiles, sLocalFiles.Data(), iNumEvents);
  }

  LoadLibs();
  LoadMacros();

  UInt_t kPrePhysSel = AliVEvent::kMB;
  //UInt_t kPrePhysSel = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  //UInt_t kPrePhysSel = AliVEvent::kEMCEGA;
  //UInt_t kPrePhysSel = AliVEvent::kEMCEJE;

  UInt_t kPhysSel = AliEmcalPhysicsSelection::kEmcalOk;
  //UInt_t kPhysSel = AliVEvent::kINT7;
  //UInt_t kPhysSel = AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  //UInt_t kPhysSel = AliVEvent::kEMCEGA;
  //UInt_t kPhysSel = AliVEvent::kEMCEJE;
  //UInt_t kPhysSel = AliEmcalPhysicsSelection::kEmcalHT;

  // Analysis manager
  AliAnalysisManager* pMgr = new AliAnalysisManager(cTaskName);

  if (iDataType == kAod) {
    AliAODInputHandler* pAODHandler = AddAODHandler();
  }
  else {   //ESD or skimmed ESD
    AliESDInputHandler* pESDHandler = AddESDHandler();
  }

  if (0) {
    AliAODHandler* pAODOutHandler = AddAODOutputHandler();
  }

  // Physics selection task
  if (1) {
    AliEmcalPhysicsSelectionTask *pPhysSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kTRUE,
                                                                              kPrePhysSel,
                                                                              5, 5, 10, kTRUE);
  }

  // Centrality task
  if (0 && iDataType == kEsd) {
    AliCentralitySelectionTask *pCentralityTask = AddTaskCentrality(kTRUE);
    pCentralityTask->SelectCollisionCandidates(kPhysSel);
  }
    
  // PID response
  if (bDoHF) {
    AliAnalysisTaskPIDResponse* pPIDtask = (AliAnalysisTaskPIDResponse*)AddTaskPIDResponse(kFALSE);
    pPIDtask->SelectCollisionCandidates(kPhysSel);
  }
  
  // Setup task
  if (bDoEmcal) {
    AliEmcalSetupTask *pSetupTask = AddTaskEmcalSetup();
  }

  // User task
  AddTaskJetAna(cDataType, cRunType, kPhysSel);
	
  if (!pMgr->InitAnalysis()) return;
  pMgr->PrintStatus();

  if (iRunType == kGrid) {
    // need some cleaning up here

    /*
    TString dataPattern;
    if (dType == kAod) {
      dataPattern = "*ESDs/pass2/AOD115/*AOD.root";
    }
    else {   // ESD
      dataPattern = "*ESDs/pass2/*ESDs.root";
    }

    TString runPeriod(runperiod);
    if (runPeriod == "LHC11h")
      runPeriod += "_2";
    AliAnalysisGrid *plugin = CreateAlienHandler(taskName.Data(), gridmode, runNumbers, runPeriod.Data(), dataPattern.Data()); 
    mgr->SetGridHandler(plugin);

    // start analysis
    cout << "Starting Analysis...";
    mgr->SetDebugLevel(2);
    mgr->StartAnalysis("grid");
    */
  }
  else {  // local, carver, pdsf

    TChain* pChain = 0;
    if (iDataType == kAod) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      if (bDoHF) {
        pChain = CreateAODChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE, "AliAOD.VertexingHF.root");
      }
      else {
        pChain = CreateAODChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE, "");
      }
    }
    else { 
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      pChain = CreateESDChain(sLocalFiles.Data(), iNumFiles, 0, kFALSE);
    }
    
    // start analysis
    Printf("Starting Analysis...");
    pMgr->SetUseProgressBar(1, 250);
    //pMgr->SetDebugLevel(2);

    // To have more debug info
    //pMgr->AddClassDebug("AliAnalysisTaskSE", AliLog::kDebug+100);
    //pMgr->AddClassDebug("AliAnalysisTaskSEDmesonsFilterCJ", AliLog::kDebug+100);
    //pMgr->AddClassDebug("AliAnalysisTaskDmesonJetCorrelations", AliLog::kDebug+100);
    //pMgr->AddClassDebug("AliEmcalJetTask", AliLog::kDebug+100);

    TFile *pOutFile = new TFile("train.root","RECREATE");
    pOutFile->cd();
    pMgr->Write();
    pOutFile->Close();
    delete pOutFile;

    pMgr->StartAnalysis("local", pChain, iNumEvents);
  }
}

//______________________________________________________________________________
void LoadLibs()
{
  // load fastjet libraries 3.x
  gSystem->Load("libCGAL");
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical");
  gSystem->Load("$FASTJET/lib/libfastjetplugins");
  gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");
}

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *cTaskName, const char *cGridMode, const TArrayI &oRunNumbers, 
				    const char *cRunPeriod, const char *cDataPattern)
{
  // need some clean up here

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetRunMode(gridmode);
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-05-1");
  plugin->SetAliROOTVersion("v5-04-41a-AN");
  plugin->AddExternalPackage("fastjet::v2.4.2");

  // Declare input data to be processed.
  
  // Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  TString datadir("/alice/data/");
  TString runP(runPeriod);
  if (runP.Contains("10"))
    datadir += "2010/";
  else if (runP.Contains("11"))
    datadir += "2011/";
  else if (runP.Contains("12"))
    datadir += "2012/";
  else {
    AliError(Form("Run period %s not recognized!", runPeriod));
    return 0;
  }
  datadir += runP;
  plugin->SetGridDataDir(datadir);

  // Set data search pattern
  plugin->SetDataPattern(dataPattern);
  plugin->SetRunPrefix("000");   // real data
  // ...then add run numbers to be considered
  for (Int_t i = 0; i < runNumbers.GetSize(); i++) {
    plugin->AddRunNumber(runNumbers[i]);
  }
  plugin->SetNrunsPerMaster(3);
  plugin->SetOutputToRunNo();
  // comment out the next line when using the "terminate" option, unless
  // you want separate merged files for each run
  plugin->SetMergeViaJDL();
  
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir(taskname);
  
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("out"); // In this case will be $HOME/taskname/out
  
  plugin->AddIncludePath("-I$ALICE_ROOT/include");
  plugin->AddIncludePath("-I$ALICE_ROOT/PWGJE/EMCALJetTasks/");
  plugin->AddIncludePath("-I$ALICE_ROOT/JETAN/fastjet -I$ALICE_ROOT/JETAN/fastjet/fastjet");
  plugin->AddIncludePath("-I$PWD/.");

  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  //plugin->SetAnalysisSource("");
  
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalRootLibs("libGui.so libXMLParser.so libMinuit2.so libProof.so");
  plugin->SetAdditionalLibs("libCDB.so libRAWDatabase.so libSTEER.so libEMCALUtils.so libEMCALraw.so "
                            "libEMCALbase.so libEMCALrec.so libPHOSUtils.so libPWGCaloTrackCorrBase.so "
                            "libPWGGACaloTrackCorrelations.so libPWGGACaloTasks.so libTRDbase.so libVZERObase.so "
                            "libVZEROrec.so libTENDER.so libTENDERSupplies.so libJETAN.so libCGAL.so libfastjet.so "
                            "libsiscone.so libSISConePlugin.so libFASTJETAN.so libPWGTools.so PWGEMCAL.par "
                            "libPWGGAEMCALTasks.so PWGJEEMCALJetTasks.par");

  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  // To only save certain files, use SetDefaultOutputs(kFALSE), and then
  // SetOutputFiles("list.root other.filename") to choose which files to save
  plugin->SetDefaultOutputs();
  //plugin->SetOutputFiles("list.root");
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro(Form("%s.C",taskname));
  
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(5);
  
  // Optionally modify the executable name (default analysis.sh)
  plugin->SetExecutable(Form("%s.sh",taskname));
  
  // set number of test files to use in "test" mode
  plugin->SetNtestFiles(1);
  
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);
  
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(50000);
  
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  plugin->SetJDLName(Form("%s.jdl",taskname));
  
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");
  
  return plugin;
}

void LoadMacros()
{
  // Aliroot macros
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODOutputHandler.C");
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCompat.C");

  gROOT->LoadMacro("addTask/AddTaskJetAna.C");
}
