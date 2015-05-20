// AddTaskEx01.C

//______________________________________________________________________________
AliAnalysisTaskEx01* AddTaskEx01()
{
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
 
  gROOT->LoadMacro("AliAnalysisTaskEx01.cxx+g");
  AliAnalysisTaskSE* task = new AliAnalysisTaskEx01(taskname);
  task->SelectCollisionCandidates(AliVEvent::kMB);
  mgr->AddTask(task);
    
  TString outfilename("AnalysisResults.root");
  
    // create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("coutput1", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
        
  // connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutput1);
        
  return task;
}
