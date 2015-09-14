// runJetSim.C

//______________________________________________________________________________
void runJetSim(UInt_t numevents = 10000)
{
  gSystem->Load("libpythia6_4_21");
  
  // analysis manager
  AliAnalysisManager* mgr = new AliAnalysisManager("JetSim");

  AliVEventHandler* esdH = new AliESDInputHandler;
  mgr->SetInputEventHandler(esdH);

  AliMCEventHandler* handler = new AliMCEventHandler;
  handler->SetReadTR(kFALSE);
  mgr->SetMCtruthEventHandler(handler);
  
  // Analysis tasks
  gROOT->LoadMacro("AddTaskJetSim.C");
  AddTaskJetSim();

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  TChain* chain = 0;
  chain = new TChain("TE");
  chain->Add("/Users/sa639/Documents/Work/ALICE/alice-yale-hfjet/sim/prodPP2010pass4SalvatoreSplitPtHard/data/kPyCharmppMNRwmi/galice.root");
  //chain->Add("/Users/sa639/Documents/Work/ALICE/alice-yale-hfjet/sim/prodPP2010pass4SalvatoreSplitPtHard/data/kPyMb/galice.root");
  //chain->Add("/Users/sa639/Documents/Work/ALICE/alice-yale-hfjet/sim/prodPP2010pass4SalvatoreSplitPtHard/data/kPyMb_trigger/galice.root");
  //chain->Add("/Users/sa639/Documents/Work/ALICE/alice-yale-hfjet/sim/prodPP2010pass4SalvatoreSplitPtHard/data/kPyCharm/galice.root");
  //chain->Add("/Users/sa639/Documents/Work/ALICE/Data/LHC12a15a/1/001/root_archive.zip#galice.root");
  //chain->Add("/Users/sa639/Documents/Work/ALICE/alice-yale-hfjet/sim/prodPP2010pass4SalvatoreSplitPtHard/galice.root");
  chain->Print();
  chain->GetListOfFiles()->Print();
  
  // start analysis
  cout << "Starting Analysis...";
  mgr->SetUseProgressBar(1, 25);
  //mgr->SetDebugLevel(1);

  //mgr->AddClassDebug("AliJetContainer",100);
  //mgr->AddClassDebug("AliAnalysisTaskSAQA",AliLog::kDebug+1);
  //mgr->AddClassDebug("AliEmcalJetTask",AliLog::kDebug+1);
  //mgr->AddClassDebug("AliJetResponseMaker",100);
  //mgr->AddClassDebug("AliEmcalContainer",100);
  //mgr->AddClassDebug("AliAnalysisTaskSEDmesonsFilterCJ", AliLog::kDebug+100);
  //mgr->AddClassDebug("AliMCHFParticleSelector", AliLog::kDebug);

  TFile *pOutFile = new TFile("train.root","RECREATE");
  pOutFile->cd();
  mgr->Write();
  pOutFile->Close();
  delete pOutFile;

  mgr->StartAnalysis("local", chain, numevents);
}
