#include "Riostream.h"
void LoadLibraries();
void AddAnalysisTasks(); 
void QAmerge(const char *, Int_t);

Int_t iCollisionType = 0; // 0=pp, 1=PbPb
// Trigger mask.

UInt_t kTriggerInt = AliVEvent::kAnyINT;
UInt_t kTriggerMuonAll = AliVEvent::kMUL7 | AliVEvent::kMUSH7 | AliVEvent::kMUU7 | AliVEvent::kMUS7
                       | AliVEvent::kMUSPB | AliVEvent::kMUSHPB | AliVEvent::kMuonLikePB | AliVEvent::kMuonUnlikePB;
UInt_t kTriggerMuonBarell = AliVEvent::kMUU7;
UInt_t kTriggerEMC   = AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMCEJE | AliVEvent::kEMCEGA;
UInt_t kTriggerHM   = AliVEvent::kHighMult;
// Main trigger mask used:
UInt_t kTriggerMask = kTriggerInt;

Int_t runNumbers[5] = {158626};

Bool_t doCDBconnect   = 1;
Bool_t doEventStat    = 1;
Bool_t doCentrality   = 0;
Bool_t doQAsym        = 0;
Bool_t doVZERO        = 1;   // there is a 2nd file
Bool_t doVZEROPbPb    = 1; 
Bool_t doVertex       = 1;
Bool_t doSPD          = 1;   // needs RP   
Bool_t doTPC          = 1;
Bool_t doHLT          = 1;
Bool_t doSDD          = 1;   // needs RP
Bool_t doSSDdEdx      = 1;

Bool_t doTRD          = 1;
Bool_t doITS          = 1;
Bool_t doITSsaTracks  = 1; 
Bool_t doITSalign     = 1;  
Bool_t doCALO         = 1;

Bool_t doMUONTrig     = 1;
Bool_t doImpParRes    = 1;
Bool_t doMUON         = 1;
Bool_t doTOF          = 1;
Bool_t doHMPID        = 1;
Bool_t doT0           = 1;
Bool_t doZDC          = 1;
Bool_t doPIDResponse  = 1;
Bool_t doPIDqa        = 1; //new
Bool_t doFMD          = 1; // new
Bool_t doPHOS         = 1; // new
Bool_t doPHOSTrig     = 1; // new
Bool_t doEMCAL        = 0;
Bool_t doFBFqa        = 1; // new - not ported yet to revision

Bool_t doMUONEff      = 0;   // NEEDS geometry
Bool_t doV0           = 0;   // NEEDS MCtruth 

               // Debug level
Int_t       debug_level        = 1;        // Debugging
Int_t       run_number = 0;

void QAtrainsim(Int_t run = 0, 
             const char *xmlfile   = "wn.xml",
             Int_t  stage          = 0, /*0 = QA train, 1...n - merging stage*/
             const char *cdb     = "raw://")
{
  run_number = run;
  TString cdbString(cdb);
  if (cdbString.Contains("raw://"))
 {
  TGrid::Connect("alien://");
  if (!gGrid || !gGrid->IsConnected()) {
    ::Error("QAtrain", "No grid connection");
    return;
  }  
  
}  
  // Set temporary merging directory to current one
  gSystem->Setenv("TMPDIR", gSystem->pwd());
  // Set temporary compilation directory to current one
  gSystem->SetBuildDir(gSystem->pwd(), kTRUE);
  // Load libraries
  LoadLibraries();
  printf("Include path: %s\n", gSystem->GetIncludePath());
  // Create manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("PilotAnalysis_sim", "Production train");
  mgr->SetRunFromPath(run_number);
  // Input handler
  AliESDInputHandlerRP *esdHandler = new AliESDInputHandlerRP();
  esdHandler->SetReadFriends(kTRUE);
  esdHandler->SetActiveBranches("ESDfriend");
  mgr->SetInputEventHandler(esdHandler);
  mgr->SetDebugLevel(debug_level);
  
  // Monte Carlo handler
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  mgr->SetMCtruthEventHandler(mcHandler);
  mcHandler->SetPreReadMode(1);
  mcHandler->SetReadTR(kTRUE); 

  // AnalysisTasks
//  mgr->Lock();
  mgr->SetFileInfoLog("fileinfo.log"); 
  AddAnalysisTasks(cdb);
//  mgr->UnLock();
//  mcHandler = (AliMCEventHandler*)mgr->GetMCtruthEventHandler(); 
//  mcHandler->SetReadTR(kTRUE); 
//  mcHandler->SetPreReadMode(1);
  if (stage>0) {
    QAmerge(xmlfile, stage);
    return;
  }   
  // Input chain
  TChain *chain = new TChain("esdTree");
  chain->Add("AliESDs.root");
  TStopwatch timer;
  timer.Start();
  if (mgr->InitAnalysis()) {                                                                                                              
    mgr->PrintStatus(); 
    mgr->SetSkipTerminate(kTRUE);
    mgr->SetNSysInfo(1);
    mgr->StartAnalysis("local", chain);
  }
  timer.Print();
}

void LoadLibraries()
{
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TRD -I$ALICE_ROOT/PWGPP -I$ALICE_ROOT/PWGPP/TRD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libESDfilter.so");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTENDER");
  gSystem->Load("libPWGPP.so");
  gSystem->Load("libAliHLTTrigger.so");

  if (doEMCAL || doPHOS || doCALO) {
     gSystem->Load("libEMCALUtils");
     gSystem->Load("libPHOSUtils");
     gSystem->Load("libPWGCaloTrackCorrBase");
     gSystem->Load("libPWGGACaloTrackCorrelations");
     gSystem->Load("libPWGGACaloTasks");
     gSystem->Load("libPWGGAPHOSTasks");
     gSystem->Load("libPWGTools");
     gSystem->Load("libPWGEMCAL");
     gSystem->Load("libPWGGAEMCALTasks");
  }  
  if(doMUON || doMUONTrig) {
     gSystem->Load("libPWGmuon");
     gSystem->Load("libPWGPPMUONlite");
     gSystem->Load("libPWGmuondep");
  }
  if (doFMD) {
     gSystem->Load("libPWGLFforward2");
  }      
}
 

void AddAnalysisTasks(const char *cdb_location)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetCommonFileName("QAresults.root");
  // Statistics task
  mgr->AddStatisticsTask(kTriggerMask);
  //
  // CDB connection
  //
  if (doCDBconnect) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    AliTaskCDBconnect *taskCDB = AddTaskCDBconnect(cdb_location, run_number);
    if (!taskCDB) return;
  }    
  
  //
  // Event Statistics (Jan Fiete)
  //
  if (doEventStat) {
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
      AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kTRUE /*MC*/);
  }
  //
  // PIDResponse(JENS)
  //
  if (doPIDResponse) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"); 
    AliAnalysisTaskPIDResponse *PIDResponse = AddTaskPIDResponse(kTRUE);
    PIDResponse->SelectCollisionCandidates(kTriggerMask);
  }  

   
  //
  // Centrality (A. Toia)
  //
  if (doCentrality) {
//     if (!iCollisionType) {
//        printf("Disabling centrality task for p-p\n");
//        doCentrality = kFALSE;
//     } else {           
        gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
        AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
        taskCentrality->SetMCInput();        
//     }   
  }   
  
  // Vertexing (A. Dainese)
  // 
  if (doVertex) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskVertexESD.C");
    // Specific setting for MC
    AliAnalysisTaskVertexESD* taskvertexesd =  AddTaskVertexESD(kTRUE, kTriggerMask);
    taskvertexesd->SelectCollisionCandidates(kTriggerMask);
  }  

  // TPC QA (E. Sicking)
  //
  if (doQAsym) {
  // offline trigger in AddTask
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskQAsym.C");
    AliAnalysisTaskSE * taskqasim = AddTaskQAsym(0, kTriggerMask, kTriggerHM, kTriggerEMC, kTriggerMuonBarell);
  }  
  //
  // VZERO QA  (C. Cheshkov)
  //
  if (doVZERO) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskVZEROQA.C");
    AliAnalysisTaskSE * taskv0qa = AddTaskVZEROQA(0);
//  taskv0qa->SelectCollisionCandidates();
  }
  if (doVZEROPbPb) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/VZERO/AddTaskVZEROPbPb.C");
    AliAnaVZEROPbPb* taskV0PbPb = (AliAnaVZEROPbPb*)AddTaskVZEROPbPb(run_number);
//    taskV0PbPb->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU_B1-,CPBI2_B1-,CPBI1WU-,CPBI1-,CVHNWU-,CVHN-,CVLNWU-,CVLN-");
//    taskV0PbPb->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU,CPBI2,CPBI1WU-,CPBI1-,CVHNWU,CVHN,CVLNWU,CVLN");
//    taskV0PbPb->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU-,CPBI2-,CPBI2WU_B1-,CPBI2_B1-,CPBI1WU-,CPBI1-,CVHNWU-,CVHN-,CVHN_R2-,CVHNWU_R2-,CVLNWU-,CVLN-,CVLN_B2-,CVLNWU_B2-");
//    taskV0PbPb->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU-,CPBI2-,CPBI2WU_B1-,CPBI2_B1-,CPBI1WU-,CPBI1-,CVHNWU-,CVHN-,CVHN_R2-,CVHNWU_R2-,CVLNWU-,CVLN-,CVLN_R1-,CVLN_B2-,CVLNWU_R1-,CVLNWU_B2-");
//    taskV0PbPb->SetClassesNames("CTRUE-,C0HWU-,CPBI2WU-,CPBI2-,CPBI2WU_B1-,CPBI2_B1-,CPBI1WU-,CPBI1-,CVHNWU-,CVHN-,CVHN_R2-,CVHNWU_R2-,CVLNWU-,CVLN-,CVLN_R1-,CVLN_B2-,CVLNWU_R1-,CVLNWU_B2-,CSEMI_R1-,CSEMIWU_R1-,CCENT_R2-,CCENTWU_R2-");
  }
  //
  // TPC (Jacek Otwinowski & Michael Knichel)
  //
  //
  // Optionally MC information can be used by setting the 1st argument to true
  // Optionally friends information can be switched off by setting the 2st argument 
  // to false
  // Optionally highMult axis can be used by setting the 3st argument to true (for PbPb)
  if (doTPC) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/TPC/macros/AddTaskPerformanceTPCdEdxQA.C");
    AliPerformanceTask *tpcQA = 0;
    if (iCollisionType) {
       // High multiplicity Pb-Pb
       tpcQA = AddTaskPerformanceTPCdEdxQA(kTRUE, kTRUE, kTRUE);
    } else {
      // Low multiplicity (pp)
       tpcQA = AddTaskPerformanceTPCdEdxQA(kTRUE, kTRUE, kFALSE);
    }
    tpcQA->SelectCollisionCandidates(kTriggerMask);
    AliPerformanceRes::SetMergeEntriesCut(5000000); 
  }  

  // HLT (Alberica Toia)
  if (doHLT) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/TPC/macros/AddTaskPerformanceTPCdEdxQA.C");
    AliPerformanceTask *hltQA = AddTaskPerformanceTPCdEdxQA(kTRUE, kTRUE, kFALSE,0,kTRUE);
    hltQA->SelectCollisionCandidates(kTriggerMask);
  }  
  //
  // SPD (A. Mastroserio)
  //
  if (doSPD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskSPDQA.C");
    AliAnalysisTaskSPD* taskspdqa = (AliAnalysisTaskSPD*)AddTaskSPDQA();
    // Request from Annalisa
    if (iCollisionType) taskspdqa->SetHeavyIonMode();
    taskspdqa->SelectCollisionCandidates(kTriggerMask);
    taskspdqa->SetOCDBInfo(run_number, "raw://");
  }  
  //
  // SDD (F. Prino)
  //
  if (doSDD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddSDDPoints.C");
    AliAnalysisTaskSE* tasksdd = AddSDDPoints();
    tasksdd->SelectCollisionCandidates(kTriggerMask);
  }
  //
  // SSD dEdx (Marek Chojnacki)
  //
  if (doSSDdEdx) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskdEdxSSDQA.C");
    AliAnalysisTaskSE* taskssddedx = AddTaskdEdxSSDQA();
    taskssddedx->SelectCollisionCandidates(kTriggerMask);
  }

  //
  // ITS
  //
  if (doITS) {
  // hardcoded non-zero trigger mask
      gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskPerformanceITS.C");
      AliAnalysisTaskITSTrackingCheck *itsQA = 0;
      AliAnalysisTaskITSTrackingCheck *itsQACent0010 = 0;
      AliAnalysisTaskITSTrackingCheck *itsQACent3050 = 0;
      AliAnalysisTaskITSTrackingCheck *itsQACent6080 = 0;
      if(iCollisionType==0) {
        itsQA = AddTaskPerformanceITS(kTRUE);
      } else {
        itsQA = AddTaskPerformanceITS(kTRUE);
        itsQACent0010 = AddTaskPerformanceITS(kTRUE,kFALSE,kFALSE,3500,10000);
        itsQACent3050 = AddTaskPerformanceITS(kTRUE,kFALSE,kFALSE,590,1570);
        itsQACent6080 = AddTaskPerformanceITS(kTRUE,kFALSE,kFALSE,70,310);
      }
  }
  //
  // ITS saTracks, align (F.Prino)
  //
  if (doITSsaTracks) {
  // offline trigger in AddTask
     gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskITSsaTracks.C");
     AliAnalysisTaskITSsaTracks *itssaTracks = AddTaskITSsaTracks(kTRUE,kFALSE);
     itssaTracks->SelectCollisionCandidates(kTriggerMask);
  }   
  if (doITSalign) {
  // no offline trigger selection
     gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskITSAlign.C");
     AliAnalysisTaskITSAlignQA *itsAlign = AddTaskITSAlign(0,2011);
  }   
  //
  // TRD (Alex Bercuci, M. Fasel) 
  //
  if(doTRD) {
  // no offline trigger selection
      gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTrainPerformanceTRD.C");
      // steer individual TRD tasks
      Bool_t 
      doCheckESD(kTRUE),  // AliTRDcheckESD
      doCheckDET(kTRUE),  // AliTRDcheckDET
      doEffic(kTRUE),     // AliTRDefficiency
      doResolution(kTRUE),// AliTRDresolution
      doCheckPID(kTRUE),  // AliTRDcheckPID
      doV0Monitor(kFALSE);// AliTRDv0Monitor
      AddTrainPerformanceTRD(Translate(doCheckESD, doCheckDET, doEffic, doResolution, doCheckPID, doV0Monitor));
  }

  //
  // ZDC (Chiara Oppedisano) 
  //
  if(doZDC) {
  // hardcoded kMB trigger mask
     gROOT->LoadMacro("$ALICE_ROOT/PWGPP/ZDC/AddTaskZDCQA.C");
     AliAnalysisTaskSE *taskZDC = AddTaskZDCQA();
     taskZDC->SelectCollisionCandidates(kTriggerMask);
  }   
  //
  // Calorimetry (Gustavo Conesa)
  //

  if(doCALO) {
        
      gROOT->LoadMacro("$ALICE_ROOT/PWGGA/CaloTrackCorrelations/macros/QA/AddTaskCalorimeterQA.C");
      AliAnalysisTaskCaloTrackCorrelation *taskCaloQA = AddTaskCalorimeterQA("default");
      taskCaloQA->SetDebugLevel(0);
      // offline mask set in AddTask to kMB
//      taskCaloQA->SelectCollisionCandidates(kTriggerMask);
      // Add a new calo task with EMC1 trigger only
//      taskCaloQA = AddTaskCalorimeterQA("trigEMC");
//      taskCaloQA->SelectCollisionCandidates(kTriggerEMC);
     
//      gROOT->LoadMacro("$ALICE_ROOT/PWGGA/CaloTrackCorrelations/macros/QA/AddTaskCalorimeterQA.C");
//      AliAnalysisTaskCaloTrackCorrelation*taskCaloQA=  AddTaskCalorimeterQA("ESD",  2012,  kFALSE,  kTRUE);
//      taskCaloQA->SetDebugLevel(0);
      // offline mask set in AddTask to kMB
//      taskCaloQA->SelectCollisionCandidates(kTriggerMask);
      // Add a new calo task with EMC1 trigger only
//      taskCaloQA=  AddTaskCalorimeterQA("ESD",  2012,  kFALSE,  kTRUE,  "",  "EMC7");
//      taskCaloQA->SetDebugLevel(0);
//      taskCaloQA->SelectCollisionCandidates(kTriggerEMC);

  
//      gROOT->LoadMacro("$ALICE_ROOT/PWGGA/CaloTrackCorrelations/macros/QA/AddTaskCalorimeterQA.C");
//       AliAnalysisTaskCaloTrackCorrelation *taskCaloQA = AddTaskCalorimeterQA("ESD", 2011, kFALSE, kTRUE);
//      taskCaloQA->SetDebugLevel(0);
      // offline mask set in AddTask to kMB
//      taskCaloQA->SelectCollisionCandidates(kTriggerMask);
      // Add a new calo task with EMC1 trigger only
//      taskCaloQA = AddTaskCalorimeterQA("ESD", 2011, kFALSE, kTRUE, "", "EMC7");
//      taskCaloQA->SetDebugLevel(0);
//      taskCaloQA->SelectCollisionCandidates(kTriggerEMC);
  }

  //
  // Muon Trigger
  //
  
  if(doMUONTrig) {
  // no offline trigger selection
      gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskMTRchamberEfficiency.C");
      AliAnalysisTaskTrigChEff *taskMuonTrig = AddTaskMTRchamberEfficiency();
  }

  //
  // Muon Efficiency (not used)
  //

  if(doMUONEff) {
      gROOT->LoadMacro("$ALICE_ROOT/PWG3/muondep/AddTaskMUONTrackingEfficiency.C");
      AliAnalysisTaskMuonTrackingEff *taskMuonTrackEff = AddTaskMUONTrackingEfficiency(kTRUE);
  }
  
  //
  // V0-Decay Reconstruction (Ana Marin) (not used)
  // 

  if (doV0) {
      gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskV0QA.C");
      AliAnalysisTaskV0QA *taskv0QA = AddTaskV0QA(kTRUE);
  }
  //
  // Impact parameter resolution (xianbao.yuan@pd.infn.it, andrea.dainese@pd.infn.it)
  //
  if (doImpParRes) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskImpParRes.C");
    AliAnalysisTaskSE* taskimpparres=0;
    // Specific setting for MC
    if(iCollisionType==0) {
       taskimpparres= AddTaskImpParRes(kTRUE);
    } else {
       taskimpparres= AddTaskImpParRes(kTRUE,-1,kTRUE,kFALSE);
    }
    taskimpparres->SelectCollisionCandidates(kTriggerMask);
  }  
  //
  // MUON QA (Philippe Pillot)
  //
  if (doMUON) {
  // trigger analysis internal
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/PilotTrain/AddTaskMuonQA.C");
    AliAnalysisTaskSE* taskmuonqa= AddTaskMuonQA(kFALSE);
  }  
  //
  // TOF (Francesca Bellini)
  //
  if (doTOF) {
  //    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/TOF/AddTaskTOFQA.C");
 //   AliAnalysisTaskTOFqa *tofQA = AddTaskTOFQA(kFALSE);
  //  tofQA->SelectCollisionCandidates(kTriggerMask);
   gROOT->LoadMacro("$ALICE_ROOT/PWGPP/TOF/AddTaskTOFqaID.C");
   AliAnalysisTaskTOFqaID *tofQA = AddTaskTOFqaID(kFALSE, AliVEvent::kAnyINT, 0, kFALSE, "", kTRUE, 0);
//   tofQA->SelectCollisionCandidates(kTriggerMask);
 } 
  //
  // PIDqa(JENS)
  //
  if (doPIDqa) {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa *PIDQA = AddTaskPIDqa();
    PIDQA->SelectCollisionCandidates(kTriggerMask);
  }  
 
  //
  // HMPID QA (Giacomo Volpe)
  //
  if (doHMPID) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/HMPID/AddTaskHmpidQA.C");
    AliAnalysisTaskSE* taskhmpidqa= AddTaskHmpidQA(kTRUE);
      // offline mask set in AddTask to kMB
    taskhmpidqa->SelectCollisionCandidates(kTriggerMask);
  }      
  // T0 QA (Alla Mayevskaya)
  if (doT0) {
  // no offline trigger selection
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/T0/AddTaskT0QA.C");
    AliT0AnalysisTaskQA* taskt0qa= AddTaskT0QA();
    taskt0qa->SelectCollisionCandidates(kTriggerMask);
  }      
  // FMD QA (Christian Holm Christiansen)
  if (doFMD) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/AddTaskForwardQA.C");
    // Parameters: usemc, usecentrality
    AliAnalysisTaskSE *forwardQA = (AliAnalysisTaskSE *)AddTaskForwardQA(kTRUE, (Bool_t)doCentrality);
    // No offline trigger config. needed (see #84077)
  }
   //     
  // PHOS QA (Boris Polishchuk)
  //
  if (doPHOS) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/CaloCellQA/macros/AddTaskCaloCellsQA.C");
    AliAnalysisTaskCaloCellsQA *taskPHOSCellQA1 = AddTaskCaloCellsQA(4, 1, NULL,"PHOSCellsQA_AnyInt"); 
    taskPHOSCellQA1->SelectCollisionCandidates(kTriggerMask);
    taskPHOSCellQA1->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.3,1.0);
    AliAnalysisTaskCaloCellsQA *taskPHOSCellQA2 = AddTaskCaloCellsQA(4, 1, NULL,"PHOSCellsQA_PHI7"); 
    taskPHOSCellQA2->SelectCollisionCandidates(AliVEvent::kPHI7);
    taskPHOSCellQA2->GetCaloCellsQA()->SetClusterEnergyCuts(0.3,0.3,1.0);
    // Pi0 QA fo PbPb
    if (iCollisionType == 1) {
      gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_PbPbQA/macros/AddTaskPHOSPbPb.C");
      AliAnalysisTaskPHOSPbPbQA* phosPbPb = AddTaskPHOSPbPbQA(0);
    }
  } 
   if (doPHOSTrig) {
     gROOT->LoadMacro("$ALICE_ROOT/PWGGA/PHOSTasks/PHOS_TriggerQA/macros/AddTaskPHOSTriggerQA.C");
     AliAnalysisTaskPHOSTriggerQA *taskPHOSTrig = AddTaskPHOSTriggerQA(NULL);
  }   
  //
  // EMCAL QA (Gustavo Conesa)
  //
  if (doEMCAL) {
     gROOT->LoadMacro("$ALICE_ROOT/PWGGA/EMCALTasks/macros/AddTaskEMCALTriggerQA.C");
     AliAnalysisTaskEMCALTriggerQA *emctrig = AddTaskEMCALTriggerQA();
  }   
  //     
  // FLOW and BF QA (C.Perez && A.Rodriguez)
  //
  if (doFBFqa) {
    gROOT->LoadMacro("$ALICE_ROOT/PWGPP/macros/AddTaskFBFqa.C");
    AliAnalysisTaskSE *qaFBFMB = (AliAnalysisTaskSE*) AddTaskFBFqa("qaFBFmb",kFALSE);
    qaFBFMB->SelectCollisionCandidates(AliVEvent::kMB);
    AliAnalysisTaskSE *qaFBFSC = (AliAnalysisTaskSE*) AddTaskFBFqa("qaFBFsc",kFALSE);
    qaFBFSC->SelectCollisionCandidates(AliVEvent::kSemiCentral);
    AliAnalysisTaskSE *qaFBFCE = (AliAnalysisTaskSE*) AddTaskFBFqa("qaFBFce",kFALSE);
    qaFBFCE->SelectCollisionCandidates(AliVEvent::kCentral);
  }
}

void QAmerge(const char *dir, Int_t stage)
{
// Merging method
  TStopwatch timer;
  timer.Start();
  TString outputDir = dir;
  TString outputFiles = "QAresults.root,EventStat_temp.root";
  TString mergeExcludes = "";
  TObjArray *list = outputFiles.Tokenize(",");
  TIter *iter = new TIter(list);
  TObjString *str;
  TString outputFile;
  Bool_t merged = kTRUE;
  while((str=(TObjString*)iter->Next())) {
    outputFile = str->GetString();
    // Skip already merged outputs
    if (!gSystem->AccessPathName(outputFile)) {
       printf("Output file <%s> found. Not merging again.",outputFile.Data());
       continue;
    }
    if (mergeExcludes.Contains(outputFile.Data())) continue;
    merged = AliAnalysisAlien::MergeOutput(outputFile, outputDir, 10, stage);
    if (!merged) {
       printf("ERROR: Cannot merge %s\n", outputFile.Data());
       continue;
    }
  }
  TString infolog = "fileinfo.log";
  AliAnalysisAlien::MergeInfo(infolog, dir); 
  // read the analysis manager from file
  if (!outputDir.Contains("Stage")) {
    ofstream out;
    out.open("outputs_valid", ios::out);
    out.close();    
    return;
  }
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  mgr->SetRunFromPath(mgr->GetRunFromAlienPath(dir));
  mgr->SetSkipTerminate(kFALSE);
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  AliLog::SetGlobalLogLevel(AliLog::kError);
  TTree *tree = NULL;
  gROOT->cd();
  mgr->StartAnalysis("gridterminate", tree);
  ofstream out;
  out.open("outputs_valid", ios::out);
  out.close();
  timer.Print();
}
