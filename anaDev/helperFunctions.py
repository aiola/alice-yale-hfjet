import ROOT
from enum import Enum

class AnaMode(Enum):
    ESD = 1
    AOD = 2

def LoadMacros():
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetQA.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetSpectraQA.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskSEDmesonsFilterCJ.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskDmesonJetCorr.C")

def AddESDHandler():
    mgr = ROOT.AliAnalysisManager.GetAnalysisManager()
    
    if mgr == None:
        ROOT.Error("AddESDHandler", "No analysis manager to connect to.")
        return None

    inputHandler = mgr.GetInputEventHandler()

    handler = ROOT.AliESDInputHandler()
  
    if inputHandler == None:
        mgr.SetInputEventHandler(handler);
    else:
      ROOT.Error("AddESDHandler", "inputHandler is NOT null. ESD handler was NOT added !!!")
      
    return handler
    
def AddAODHandler():
    mgr = ROOT.AliAnalysisManager.GetAnalysisManager()
    
    if mgr == None:
        ROOT.Error("AddAODHandler", "No analysis manager to connect to.")
        return None

    inputHandler = mgr.GetInputEventHandler()

    handler = ROOT.AliAODInputHandler()
  
    if inputHandler == None:
        mgr.SetInputEventHandler(handler);
    else:
      ROOT.Error("AddAODHandler", "inputHandler is NOT null. AOD handler was NOT added !!!")
      
    return handler

def AddTaskPIDResponse(isMC=False, autoMCesd=True, tuneOnData=True, recoPass=2, cachePID=False, detResponse="",
                       useTPCEtaCorrection=True,             #Please use default value! Otherwise splines can be off
                       useTPCMultiplicityCorrection=True,    #Please use default value! Otherwise splines can be off
                       recoDataPass=-1):
    
    # Macro to connect a centrality selection task to an existing analysis manager.
    mgr = ROOT.AliAnalysisManager.GetAnalysisManager()
    
    if mgr == None:
        ROOT.Error("AddTaskPIDResponse", "No analysis manager to connect to.")
        return None

    inputHandler = mgr.GetInputEventHandler()
    
    print "========================================================================================"
    print "PIDResponse: Initialising AliAnalysisTaskPIDResponse"

    pidTask = ROOT.AliAnalysisTaskPIDResponse("PIDResponseTask")
    pidTask.SetIsMC(isMC)
    if isMC:
        if tuneOnData:
            print "             Using MC with tune on data."
            print "             !!! ATTENTION ATTENTION ATTENTION !!!"
            print("             You MUST make sure the reco pass set (", recoPass, ") corresponds to the one this MC was produced for!")
            pidTask.SetTuneOnData(kTRUE,recoPass)
            # tuning on MC is by default active on TPC and TOF, to enable it only on one of them use:
            # pidTask->SetTuneOnDataMask(AliPIDResponse::kDetTPC);
            # pidTask->SetTuneOnDataMask(AliPIDResponse::kDetTOF);
        else:
            print "             !!! ATTENTION ATTENTION ATTENTION !!!"
            print "             You are using MC without the tune on data option."
            print "             NOTE that this is not supported any longer!."
            print "             !!! ATTENTION ATTENTION ATTENTION !!!"

    pidTask.SetCachePID(cachePID)
    pidTask.SetSpecialDetectorResponse(detResponse)
    pidTask.SetUseTPCEtaCorrection(useTPCEtaCorrection)
    pidTask.SetUseTPCMultiplicityCorrection(useTPCMultiplicityCorrection)
    pidTask.SetUserDataRecoPass(recoDataPass)
    mgr.AddTask(pidTask)

    mgr.ConnectInput(pidTask, 0, mgr.GetCommonInputContainer())
      
    print "========================================================================================"

    return pidTask
