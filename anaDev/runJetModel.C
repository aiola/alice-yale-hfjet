// runJetModel.C

class AliAnalysisGrid;

//______________________________________________________________________________
void runJetModel(
	 const char   *datatype     = "AOD",                       // set the analysis type, AOD, ESD or sESD
         const char   *runtype      = "local",                     // local or grid
         const char   *gridmode     = "test",                      // set the grid run mode (can be "full", "test", "offline", "submit" or "terminate")
	 const char   *localfiles   = "",                          // set the local list file
	 UInt_t        fileshift    = 0,                           // list file shift
	 UInt_t        numfiles     = 50,                          // number of files analyzed locally
	 UInt_t        numevents    = 1234567890,                  // number of events to be analyzed
	 const char   *runperiod    = "LHC12a15e_fix",             // set the run period
         const char   *taskname     = "JetModel"                   // sets name of grid generated macros
         )
{
  //gSystem->SetFPEMask(TSystem::kInvalid | TSystem::kDivByZero | TSystem::kOverflow | TSystem::kUnderflow);
  gSystem->Setenv("ETRAIN_ROOT", "../../emcaltrain");

  //Int_t run = 170040;
  Int_t run = 169838;
  Int_t nrunnumbers = 10;
  Int_t runnumbers[] = {1,2,3,4,5,6,7,8,9,10};

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

  cout << runtype << " " << datatype << " analysis chosen" << endl;

  if (rType == kGrid && dType == kSesd) {
    cout << "Skimmed ESD analysis not available on the grid!" << endl;
    return;
  }

  TString localFiles(localfiles);
  if (rType == kLocal || rType == kPdsf || rType == kCarver) {
    if (!strcmp(localfiles, "")) {
      switch (dType) {
      case (kAod):
	localFiles = "files_aod129.txt";
	break;
      case (kEsd):
	localFiles = "files_fixprod_169838.txt";
	break;
      default:    // skimmed ESD analysis
	localFiles = "files_sesd.txt";
      }
    }
    cout << "setting local analysis for " << numfiles << " files from list " << localFiles << ", max events = " << numevents << endl;
  }

  LoadLibs();

  UInt_t prePhysSel = AliVEvent::kAny;
  //UInt_t prePhysSel = AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  //UInt_t prePhysSel = AliVEvent::kEMCEGA;
  //UInt_t prePhysSel = AliVEvent::kEMCEJE;

  UInt_t physSel = AliVEvent::kAny;
  //UInt_t physSel = AliEmcalPhysicsSelection::kEmcalOk;
  //UInt_t physSel = AliVEvent::kINT7;
  //UInt_t physSel = AliVEvent::kAnyINT | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  //UInt_t physSel = AliVEvent::kEMCEGA;
  //UInt_t physSel = AliVEvent::kEMCEJE;
  //UInt_t physSel = AliEmcalPhysicsSelection::kEmcalHT;

  // analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(taskname);

  if (dType == kAod) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
    AliAODInputHandler* aodH = AddAODHandler();
  }
  else {   //ESD or skimmed ESD
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
    AliESDInputHandler* esdH = AddESDHandler();

    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCHandler.C");
    AliMCEventHandler* mcH = AddMCHandler();
  }

  if (0) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODOutputHandler.C");
    AliAODHandler* aodoutHandler = AddAODOutputHandler();
  }

  AliPhysicsSelectionTask *physSelTask = 0;

  // Physics selection task
  if (0) {
    if (0) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
      AliPhysicsSelectionTask *physSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kFALSE,
									  prePhysSel,
									  5, 5, 10, kTRUE, -1, -1, -1, -1);
    }
    else {
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(kTRUE);
    }
   
    if (!physSelTask) {
      cout << "no physSelTask"; 
      return; 
    }
  }

  // Centrality task
  if (dType == kEsd) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *centralityTask = AddTaskCentrality(kTRUE);
    centralityTask->SelectCollisionCandidates(physSel);
  }

  // Setup task
  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
  AliEmcalSetupTask *setupTask = AddTaskEmcalSetup();

  // Compatibility task, only needed for skimmed ESD
  if (dType == kSesd) {
    gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/AddTaskEmcalCompat.C");
    AliEmcalCompatTask *comptask = AddTaskEmcalCompat();
    comptask->SelectCollisionCandidates(physSel);

    gROOT->LoadMacro("$ETRAIN_ROOT/saiola/AliPicoTrackFixer.cxx+g");
    gROOT->LoadMacro("$ETRAIN_ROOT/saiola/AddTask1PicoTrackFixer.C");
    AliPicoTrackFixer *fixertask = AddTask1PicoTrackFixer();
    fixertask->SelectCollisionCandidates(physSel);
  }

  if (1) {
    gROOT->LoadMacro("$ETRAIN_ROOT/saiola/AddTaskJetModel.C");
    AddTaskJetModel(kFALSE, datatype, runtype);
  }

  if (0) {
    gROOT->LoadMacro("$ETRAIN_ROOT/saiola/AliAnalysisTaskTwoPartCorr.cxx+g");
    gROOT->LoadMacro("$ETRAIN_ROOT/saiola/AddTaskTwoPartCorr.C");
    AddTaskTwoPartCorr();
  }
	
  mgr->SetUseProgressBar(1, 25);
	
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();

  if (rType == kGrid) {

    TString taskName(Form("%s_0",taskname));

    TArrayI runNumbers(nrunnumbers, runnumbers+1);

    TString dataPattern;
    if (dType == kAod) {
      dataPattern = "*/AliAOD.root";
    }
    else {   // ESD
      dataPattern = "*/AliESDs.root";
    }
    TString subDir = Form("/%d",run);

    AliAnalysisGrid *plugin = CreateAlienHandler(taskname, gridmode, runNumbers, runperiod, dataPattern.Data(), subDir.Data()); 
    mgr->SetGridHandler(plugin);

    // start analysis
    cout << "Starting Analysis...";
    mgr->SetDebugLevel(0);
    mgr->StartAnalysis("grid");
  }
  else {  // local, carver, pdsf

    TChain* chain = 0;
    if (dType == kAod) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/CreateAODChain.C");
      chain = CreateAODChain(localFiles.Data(), numfiles, fileshift, kFALSE);
    }
    else if (dType == kEsd) { 
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/CreateESDChain.C");
      chain = CreateESDChain(localFiles.Data(), numfiles, fileshift, kFALSE);
    }
    else { // skimmed ESD 
      gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/CreateESDChain.C");
      chain = CreateESDChain(localFiles.Data(), numfiles, fileshift, kFALSE);
    }
    
    // start analysis
    cout << "Starting Analysis...";
    mgr->SetDebugLevel(2);
    //mgr->AddClassDebug("AliEmcalPhysicsSelection",100);
    //mgr->AddClassDebug("AliJetEmbeddingFromPYTHIATask",AliLog::kDebug+3);
    //mgr->AddClassDebug("AliJetResponseMaker",100);
    //mgr->AddClassDebug("AliAnalysisTaskSAQA",100);
    //mgr->AddClassDebug("AliAnalysisTaskDeltaPt",100);
    //mgr->AddClassDebug("AliEmcalJetTask",100);
    //mgr->AddClassDebug("AliJetConstituentTagCopier",AliLog::kDebug+1);
    //mgr->AddClassDebug("AliEMCALTenderSupply",AliLog::kDebug+1);
    mgr->StartAnalysis("local", chain, numevents);
  }
}

//______________________________________________________________________________
void LoadLibs()
{
  // load ROOT libraries
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");

  gSystem->Load("libMinuit");
  gSystem->Load("libGui");
  gSystem->Load("libXMLParser");
  gSystem->Load("libMinuit2");
  gSystem->Load("libProof");

  // load AliRoot libraries
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libSTEER");
  gSystem->Load("libCORRFW");
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libPWGCaloTrackCorrBase");
  gSystem->Load("libPWGGACaloTrackCorrelations");
  gSystem->Load("libPWGGACaloTasks");
  gSystem->Load("libPWGLFforward2");
  gSystem->Load("libEMCALraw");
  gSystem->Load("libEMCALbase");
  gSystem->Load("libEMCALrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libTENDER");   
  gSystem->Load("libTENDERSupplies"); 
 
  // load fastjet libraries
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical");
  gSystem->Load("$FASTJET/lib/libfastjetplugins");

  // Aliroot jet libs
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGEMCAL");
  gSystem->Load("libPWGGAEMCALTasks");
  gSystem->Load("libPWGJE");
  gSystem->Load("libPWGJEEMCALJetTasks");
   
  // include path
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGJE/EMCALJetTasks/");
  gSystem->AddIncludePath("-I$FASTJET/include -I$FASTJET/include/fastjet");
  gSystem->AddIncludePath("-I$PWD/.");
}

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, const TArrayI &runNumbers, 
				    const char *runPeriod = "LHC12a15e", const char *dataPattern = "*/*AOD.root", const char* subdir="")
{
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->SetRunMode(gridmode);
  
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-05");
  plugin->SetAliROOTVersion("v5-04-65-AN");
  plugin->AddExternalPackage("fastjet::v2.4.2");

  // Declare input data to be processed.
  
  // Create automatically XML collections using alien 'find' command.
  // Define production directory LFN
  TString datadir;
  if (newprod) {
    datadir = newprod;
  }
  else {
    datadir = "/alice/sim/";
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
  }

  datadir += subdir;
  plugin->SetGridDataDir(datadir);

  // Set data search pattern
  plugin->SetDataPattern(dataPattern);
  plugin->SetRunPrefix("");   // sim
  // ...then add run numbers to be considered
  for (Int_t i = 0; i < runNumbers.GetSize(); i++) {
    plugin->AddRunNumber(runNumbers[i]);
  }
  plugin->SetNrunsPerMaster(1);
  plugin->SetOutputToRunNo();
  // comment out the next line when using the "terminate" option, unless
  // you want separate merged files for each run
  plugin->SetMergeViaJDL();
  
  // Define alien work directory where all files will be copied. Relative to alien $HOME.
  plugin->SetGridWorkingDir(taskname);
  
  // Declare alien output directory. Relative to working directory.
  plugin->SetGridOutputDir("out"); // In this case will be $HOME/taskname/out
  
  plugin->AddIncludePath("-I$PWD/.");
  //plugin->AddIncludePath("-I$PWD/PWGEMCAL/EMCAL/.");
  plugin->AddIncludePath("-I$PWD/PWGJEEMCALJetTasks/EMCALJetTasks/.");
  plugin->AddIncludePath("-I$ALICE_ROOT/include");
  //plugin->AddIncludePath("-I$ALICE_ROOT/PWGJE/EMCALJetTasks/");
  plugin->AddIncludePath("-I$ALICE_ROOT/JETAN/fastjet -I$ALICE_ROOT/JETAN/fastjet/fastjet");

  // Declare the analysis source files names separated by blancs. To be compiled runtime
  // using ACLiC on the worker nodes.
  //plugin->SetAnalysisSource("");
  
  // Declare all libraries (other than the default ones for the framework. These will be
  // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
  plugin->SetAdditionalRootLibs("libGui.so libXMLParser.so libMinuit2.so libProof.so");
  plugin->SetAdditionalLibs("libCDB.so libRAWDatabase.so libSTEER.so libEMCALUtils.so libEMCALraw.so libEMCALbase.so libEMCALrec.so libPHOSUtils.so libPWGCaloTrackCorrBase.so libPWGGACaloTrackCorrelations.so libPWGGACaloTasks.so libTRDbase.so libVZERObase.so libVZEROrec.so libTENDER.so libTENDERSupplies.so libJETAN.so libCGAL.so libfastjet.so libsiscone.so libSISConePlugin.so libFASTJETAN.so libPWGTools.so libPWGEMCAL.so libPWGGAEMCALTasks.so libPWGJEEMCALJetTasks.so");

  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  // To only save certain files, use SetDefaultOutputs(kFALSE), and then
  // SetOutputFiles("list.root other.filename") to choose which files to save
  plugin->SetDefaultOutputs();
  //plugin->SetOutputFiles("list.root");
  
  // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
  plugin->SetAnalysisMacro(Form("%s.C",taskname));
  
  // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
  plugin->SetSplitMaxInputFileNumber(8);
  
  // Optionally modify the executable name (default analysis.sh)
  plugin->SetExecutable(Form("%s.sh",taskname));
  
  // set number of test files to use in "test" mode
  plugin->SetNtestFiles(2);
  
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
