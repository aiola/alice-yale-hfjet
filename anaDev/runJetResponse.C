// runJetResponse.C

class AliAnalysisGrid;

//______________________________________________________________________________
void runJetResponse( 
	 const char   *datatype     = "AOD",                                       // set the analysis type, AOD, ESD or sESD
         const char   *runtype      = "local",                                     // local or grid
         const char   *gridmode     = "test",                                      // set the grid run mode (can be "full", "test", "offline", "submit" or "terminate")
	 const char   *localfiles   = "fileLists/files_LHC10f7a_fix_AOD136a.txt",  // set the local list file
	 UInt_t        numfiles     = 50,                                          // number of files analyzed locally
	 UInt_t        numevents    = 500000,                                       // number of events to be analyzed
	 const char   *runperiod    = "LHC10f7a",                                  // set the run period
         const char   *taskname     = "JetResponse",                               // sets name of grid generated macros
         Bool_t        doEmcal      = kFALSE
         )
{
  //gSystem->SetFPEMask(TSystem::kInvalid | TSystem::kDivByZero | TSystem::kOverflow | TSystem::kUnderflow);
  
  Int_t nrunnumbers = 1;
  Int_t runnumbers[] = {117112};

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

  cout << runtype << " " << datatype << " analysis chosen" << endl;

  TString localFiles(localfiles);
  if (rType == kLocal) {
    if (!strcmp(localfiles, "")) {
      Printf("You need to specify the local list files.");
      return;
    }
    cout << "setting local analysis for " << numfiles << " files from list " << localFiles << endl;
  }

  LoadLibs();

  // analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager(taskname);

  if (dType == kAod) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddAODHandler.C");
    AliAODInputHandler* aodH = AddAODHandler();
  }
  else {   //ESD
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
    AliESDInputHandler* esdH = AddESDHandler();

    if (0) { // MC
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCHandler.C");
      AliMCEventHandler* mcH = AddMCHandler();

      // Centrality selection task
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
      AliCentralitySelectionTask *centTask = AddTaskCentrality(kTRUE);

      // Physics selection task
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask *pSelTask = AddTaskPhysicsSelection(kTRUE, kTRUE, kTRUE);
    }

    if (1) { // Data
      // Centrality selection task
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
      AliCentralitySelectionTask *centTask = AddTaskCentrality(kTRUE);

      UInt_t prePhysSel = AliVEvent::kMB;

      // Physics selection task
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
      AliPhysicsSelectionTask *physSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kFALSE,
                                                                          prePhysSel,
                                                                          5, 5, 10, kTRUE, -1, -1, -1, -1);
    }
  }

  // Setup task
  if (doEmcal) {
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
    AliEmcalSetupTask *setupTask = AddTaskEmcalSetup();
  }

  // PID response
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse* pPIDtask = (AliAnalysisTaskPIDResponse*)AddTaskPIDResponse(kTRUE);

  // Analysis tasks
  if (1) {
    gROOT->LoadMacro("addTask/AddTaskJetResp.C");
    AddTaskJetResp(datatype, "local");
  }
	
  if (!mgr->InitAnalysis()) 
    return;
  mgr->PrintStatus();

  if (rType == kGrid) {

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
  else {  // local

    TChain* chain = 0;
    if (dType == kAod) {
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateAODChain.C");
      chain = CreateAODChain(localFiles.Data(), numfiles, 0, kFALSE, "AliAOD.VertexingHF.root");
    }
    else {  // ESD or skimmed ESD
      gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
      chain = CreateESDChain(localFiles.Data(), numfiles, 0, kFALSE);
    }

    // start analysis
    cout << "Starting Analysis...";
    mgr->SetUseProgressBar(1, 25);
    //mgr->SetDebugLevel(2);
    
    //mgr->AddClassDebug("AliJetContainer",100);
    //mgr->AddClassDebug("AliAnalysisTaskSAQA",AliLog::kDebug+1);
    //mgr->AddClassDebug("AliEmcalJetTask",AliLog::kDebug+1);
    //mgr->AddClassDebug("AliJetResponseMaker",100);
    //mgr->AddClassDebug("AliEmcalContainer",100);
    //mgr->AddClassDebug("AliAnalysisTaskSEDmesonsFilterCJ", AliLog::kDebug+100);
    //mgr->AddClassDebug("AliEmcalMCTrackSelector", AliLog::kDebug+5);

    TFile *pOutFile = new TFile("train.root","RECREATE");
    pOutFile->cd();
    mgr->Write();
    pOutFile->Close();
    delete pOutFile;

    mgr->StartAnalysis("local", chain, numevents);
  }
}

//______________________________________________________________________________
void LoadLibs()
{
  gSystem->Load("libCGAL");
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical");
  gSystem->Load("$FASTJET/lib/libfastjetplugins");
  gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");
  
  // Aliroot jet libs
  gSystem->Load("libPWGJE");
  gSystem->Load("libPWGJEEMCALJetTasks");
   
  // include path
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/PWGJE/EMCALJetTasks/");
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

