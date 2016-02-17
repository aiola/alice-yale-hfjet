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
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJetQA.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWGJE/FlavourJetTasks/macros/AddTaskDmesonJets.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEMCALTender.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskClusterizerFast.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusterMaker.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalClusTrackMatcher.C")
    ROOT.gROOT.LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskHadCorr.C")

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

def PrepareEMCAL(kPhysSel=ROOT.AliVEvent.kMB) :
    # Tender
    bDistBC         = False #switch for recalculation cluster position from bad channel
    bRecalibClus    = False
    bRecalcClusPos  = False
    bNonLinearCorr  = False
    bRemExoticCell  = False
    bRemExoticClus  = False
    bFidRegion      = False
    bCalibEnergy    = True
    bCalibTime      = True
    bRemBC          = True
    iNonLinFunct    = ROOT.AliEMCALRecoUtils.kNoCorrection
    bReclusterize   = False
    fSeedThresh     = 0.1      # 100 MeV
    fCellThresh     = 0.05     # 50 MeV
    iClusterizer    = ROOT.AliEMCALRecParam.kClusterizerv2
    bTrackMatch     = False
    bUpdateCellOnly = True
    fEMCtimeMin     = -50e-6
    fEMCtimeMax     =  50e-6
    fEMCtimeCut     =  1e6

    pTenderTask = ROOT.AddTaskEMCALTender(bDistBC, bRecalibClus, bRecalcClusPos, bNonLinearCorr, bRemExoticCell, bRemExoticClus,
                                          bFidRegion, bCalibEnergy, bCalibTime, bRemBC, iNonLinFunct, bReclusterize, fSeedThresh,
                                          fCellThresh, iClusterizer, bTrackMatch, bUpdateCellOnly, fEMCtimeMin, fEMCtimeMax, fEMCtimeCut)
    pTenderTask.SelectCollisionCandidates(kPhysSel)

    # Clusterizer
    pClusterizerTask = ROOT.AddTaskClusterizerFast("ClusterizerFast", "", "", iClusterizer, 0.05, 0.1,
                                                   fEMCtimeMin, fEMCtimeMax, fEMCtimeCut, False, False)
    pClusterizerTask.SelectCollisionCandidates(kPhysSel)

    bRemExoticClus  = True
    iNonLinFunct    = ROOT.AliEMCALRecoUtils.kBeamTestCorrected

    # Cluster maker
    pClusterMakerTask = ROOT.AddTaskEmcalClusterMaker(iNonLinFunct, bRemExoticClus, "usedefault", "", 0., False)
    pClusterMakerTask.GetClusterContainer(0).SetClusPtCut(0.)
    pClusterMakerTask.GetClusterContainer(0).SetClusECut(0.)
    pClusterMakerTask.SelectCollisionCandidates(kPhysSel)

    # Cluster-track matcher task
    pMatcherTask = ROOT.AddTaskEmcalClusTrackMatcher("usedefault", "usedefault", 0.1, False, True, True, False)
    pMatcherTask.SelectCollisionCandidates(kPhysSel)
    pMatcherTask.GetParticleContainer(0).SetParticlePtCut(0.15)
    pMatcherTask.GetClusterContainer(0).SetClusNonLinCorrEnergyCut(0.15)
    pMatcherTask.GetClusterContainer(0).SetClusECut(0.)
    pMatcherTask.GetClusterContainer(0).SetClusPtCut(0.)

    # Hadronic correction task
    pHadCorrTask = ROOT.AddTaskHadCorr("usedefault", "usedefault", "", 2.0, 0.15, 0.030, 0.015, 0, True, False)
    pHadCorrTask.SelectCollisionCandidates(kPhysSel)
    pHadCorrTask.GetClusterContainer(0).SetClusNonLinCorrEnergyCut(0.15)
    pHadCorrTask.GetClusterContainer(0).SetClusECut(0)
    pHadCorrTask.GetClusterContainer(0).SetClusPtCut(0.)
