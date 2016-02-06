// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++")

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TH1.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include <TGeoGlobalMagField.h>
#include "STEER/STEER/AliRunLoader.h"
#include "STEER/STEER/AliRun.h"
#include "STEER/STEER/AliConfig.h"
#include "STEER/STEER/AliSimulation.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "PYTHIA6/AliGenPythia.h"
#include "PYTHIA6/AliGenPythiaPlus.h"
#include "TDPMjet/AliGenDPMjet.h"
#include "EVGEN/AliGenParam.h"
#include "EVGEN/AliGenMUONlib.h"
#include "EVGEN/AliGenPHOSlib.h"
#include "EVGEN/AliGenSTRANGElib.h"
#include "EVGEN/AliGenCocktail.h"
#include "EVGEN/AliGenBox.h"
#include "THijing/AliGenHijing.h"
#include "TAmpt/AliGenAmpt.h"
#include "TUHKMgen/AliGenUHKM.h"
#include "EVGEN/AliGenSlowNucleons.h"
#include "EVGEN/AliSlowNucleonModel.h"
#include "EVGEN/AliSlowNucleonModelExp.h"
#include "STEER/STEERBase/AliMagF.h"
#include "STRUCT/AliBODY.h"
#include "STRUCT/AliMAG.h"
#include "STRUCT/AliABSOv3.h"
#include "STRUCT/AliDIPOv3.h"
#include "STRUCT/AliHALLv3.h"
#include "STRUCT/AliFRAMEv2.h"
#include "STRUCT/AliSHILv3.h"
#include "STRUCT/AliPIPEv3.h"
#include "ITS/AliITSv11.h"
#include "TPC/AliTPCv2.h"
#include "TOF/AliTOFv6T0.h"
#include "HMPID/AliHMPIDv3.h"
#include "ZDC/AliZDCv4.h"
#include "TRD/AliTRDv1.h"
#include "TRD/AliTRDgeometry.h"
#include "FMD/AliFMDv1.h"
#include "MUON/AliMUONv1.h"
#include "PHOS/AliPHOSv1.h"
#include "PHOS/AliPHOSSimParam.h"
#include "PMD/AliPMDv1.h"
#include "T0/AliT0v1.h"
#include "EMCAL/AliEMCALv2.h"
#include "ACORDE/AliACORDEv1.h"
#include "VZERO/AliVZEROv7.h"
#endif

//--- Functions ---
class AliGenPythia;

//========================//
// Set Random Number seed //
//========================//
TDatime dt;
static UInt_t seed = dt.Get();

// Comment line
static TString comment;

void Config()
{
//Random Number seed
  if (gSystem->Getenv("CONFIG_SEED")) {
    seed = atoi(gSystem->Getenv("CONFIG_SEED"));
  }

  gRandom->SetSeed(seed);
  cerr << "Seed for random number generation= " << seed << endl;

//Libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6_4_25");
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libgeant321");
#endif

  new TGeant3TGeo("C++ Interface to Geant3");

//=======================================================================
//Create the output file

  AliRunLoader* rl=0x0;

  cout << "Config.C: Creating Run Loader ..." << endl;
  rl = AliRunLoader::Open("galice.root", AliConfig::GetDefaultEventFolderName(), "recreate");
  if (rl == 0x0)
  {
    gAlice->Fatal("Config.C","Can not instatiate the Run Loader");
    return;
  }
  rl->SetCompressionLevel(2);
  rl->SetNumberOfEventsPerFile(1000);
  gAlice->SetRunLoader(rl);
  // gAlice->SetGeometryFromFile("geometry.root");
  // gAlice->SetGeometryFromCDB();

  // Set the trigger configuration: proton-proton
  AliSimulation::Instance()->SetTriggerConfig("p-p"); 

  //
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV


  gMC->SetProcess("DCAY",1);
  gMC->SetProcess("PAIR",1);
  gMC->SetProcess("COMP",1);
  gMC->SetProcess("PHOT",1);
  gMC->SetProcess("PFIS",0);
  gMC->SetProcess("DRAY",0);
  gMC->SetProcess("ANNI",1);
  gMC->SetProcess("BREM",1);
  gMC->SetProcess("MUNU",1);
  gMC->SetProcess("CKOV",1);
  gMC->SetProcess("HADR",1);
  gMC->SetProcess("LOSS",2);
  gMC->SetProcess("MULS",1);
  gMC->SetProcess("RAYL",1);

  Float_t cut = 1.e-3;        // 1MeV cut by default
  Float_t tofmax = 1.e10;

  gMC->SetCut("CUTGAM", cut);
  gMC->SetCut("CUTELE", cut);
  gMC->SetCut("CUTNEU", cut);
  gMC->SetCut("CUTHAD", cut);
  gMC->SetCut("CUTMUO", cut);
  gMC->SetCut("BCUTE",  cut);
  gMC->SetCut("BCUTM",  cut);
  gMC->SetCut("DCUTE",  cut);
  gMC->SetCut("DCUTM",  cut);
  gMC->SetCut("PPCUTM", cut);
  gMC->SetCut("TOFMAX", tofmax);

  const char *phwgproc = gSystem->Getenv("PWHGPROC");
  
  //======================//
  // Set External decayer //
  //======================//
  TVirtualMCDecayer* decayer = new AliDecayerPythia();
  if (!strcmp(phwgproc,"charm") || !strcmp(phwgproc,"beauty")) {
    decayer->SetForceDecay(kSemiElectronic);
  } else {
    decayer->SetForceDecay(kAll);
  }
  decayer->Init();
  gMC->SetExternalDecayer(decayer);

  //=========================//
  // Generator Configuration //
  //=========================//
  AliGenerator* gener = 0x0;
  comment = comment.Append(" pp: POWHEG + PYTHIA ");
  
  AliGenPythia *pythia = new AliGenPythia(-1);
  pythia->SetTune(350);
  
  comment = comment.Append(phwgproc);
  
  if (!strcmp(phwgproc,"charm")) {
    pythia->SetProcess(kPyCharmPWHG);
    cout << "Generating charm..." << endl;
  } else if (!strcmp(phwgproc,"beauty")) {
    pythia->SetProcess(kPyBeautyPWHG);
    cout << "Generating beauty..." << endl;
  } else if (!strcmp(phwgproc,"dijet")) {
    pythia->SetProcess(kPyJetsPWHG);
    cout << "Generating jets..." << endl;
  } else {
    cout << "Process " << phwgproc << " does not exist, please check!\n";
  }
  
  char *lhefile = Form("pwgevents.lhe");
  
  pythia->SetReadLHEF(lhefile);
  //pythia->SetStrucFunc(kCTEQ6m);
  //  pythia->SetEnergyCMS(7000.);
  pythia->SetProjectile("p", 1, 1);
  pythia->SetTarget(    "p", 1, 1);
  pythia->SetMomentumRange(0, 999999.);
  pythia->SetThetaRange(0., 180.);
  pythia->SetYRange(-12.,12.);
  pythia->SetPtRange(0,1000.);
//   pythia->SetEventListRange(0,4);
  pythia->SetTrackingFlag(0);
  
  gener = pythia;
  
  gener->SetVertexSmear(kPerEvent);
  gener->Init();

  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG, AliMagF::kBeamTypepp));

  printf("\n \n Comment: %s \n \n", comment.Data());

  rl->CdGAFile();

  Int_t iABSO  = 0;
  Int_t iACORDE= 0;
  Int_t iDIPO  = 0;
  Int_t iEMCAL = 0;
  Int_t iFMD   = 0;
  Int_t iFRAME = 0;
  Int_t iHALL  = 0;
  Int_t iITS   = 0;
  Int_t iMAG   = 1;
  Int_t iMUON  = 0;
  Int_t iPHOS  = 0;
  Int_t iPIPE  = 0;
  Int_t iPMD   = 0;
  Int_t iHMPID = 0;
  Int_t iSHIL  = 0;
  Int_t iT0    = 0;
  Int_t iTOF   = 0;
  Int_t iTPC   = 0;
  Int_t iTRD   = 0;
  Int_t iVZERO = 0;
  Int_t iZDC   = 0;


    //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");


  if (iMAG)
  {
        //=================== MAG parameters ============================
        // --- Start with Magnet since detector layouts may be depending ---
        // --- on the selected Magnet dimensions ---
    AliMAG *MAG = new AliMAG("MAG", "Magnet");
  }


  if (iABSO)
  {
        //=================== ABSO parameters ============================
    AliABSO *ABSO = new AliABSOv3("ABSO", "Muon Absorber");
  }

  if (iDIPO)
  {
        //=================== DIPO parameters ============================

    AliDIPO *DIPO = new AliDIPOv3("DIPO", "Dipole version 3");
  }

  if (iHALL)
  {
        //=================== HALL parameters ============================

    AliHALL *HALL = new AliHALLv3("HALL", "Alice Hall");
  }


  if (iFRAME)
  {
        //=================== FRAME parameters ============================

    AliFRAMEv2 *FRAME = new AliFRAMEv2("FRAME", "Space Frame");
    FRAME->SetHoles(1);
  }

  if (iSHIL)
  {
        //=================== SHIL parameters ============================

    AliSHIL *SHIL = new AliSHILv3("SHIL", "Shielding Version 3");
  }


  if (iPIPE)
  {
        //=================== PIPE parameters ============================

    AliPIPE *PIPE = new AliPIPEv3("PIPE", "Beam Pipe");
  }

  if (iITS)
  {
        //=================== ITS parameters ============================

    AliITS *ITS  = new AliITSv11("ITS","ITS v11");
  }

  if (iTPC)
  {
      //============================ TPC parameters =====================

    AliTPC *TPC = new AliTPCv2("TPC", "Default");
  }


  if (iTOF) {
        //=================== TOF parameters ============================

    AliTOF *TOF = new AliTOFv6T0("TOF", "normal TOF");
  }


  if (iHMPID)
  {
        //=================== HMPID parameters ===========================

    AliHMPID *HMPID = new AliHMPIDv3("HMPID", "normal HMPID");
  }


  if (iZDC)
  {
        //=================== ZDC parameters ============================

    AliZDC *ZDC = new AliZDCv4("ZDC", "normal ZDC");
  }

  if (iTRD)
  {
        //=================== TRD parameters ============================

    AliTRD *TRD = new AliTRDv1("TRD", "TRD slow simulator");
    AliTRDgeometry *geoTRD = TRD->GetGeometry();
	// Partial geometry: modules at 0,1,7,8,9,16,17
_	// starting at 3h in positive direction
    geoTRD->SetSMstatus(2,0);
    geoTRD->SetSMstatus(3,0);
    geoTRD->SetSMstatus(4,0);
    geoTRD->SetSMstatus(5,0);
    geoTRD->SetSMstatus(6,0);
    geoTRD->SetSMstatus(11,0);
    geoTRD->SetSMstatus(12,0);
    geoTRD->SetSMstatus(13,0);
    geoTRD->SetSMstatus(14,0);
    geoTRD->SetSMstatus(15,0);
    geoTRD->SetSMstatus(16,0);
  }

  if (iFMD)
  {
        //=================== FMD parameters ============================

    AliFMD *FMD = new AliFMDv1("FMD", "normal FMD");
  }

  if (iMUON)
  {
        //=================== MUON parameters ===========================
        // New MUONv1 version (geometry defined via builders)

    AliMUON *MUON = new AliMUONv1("MUON", "default");
  }

  if (iPHOS)
  {
        //=================== PHOS parameters ===========================

    AliPHOS *PHOS = new AliPHOSv1("PHOS", "noCPV_Modules123");
  }


  if (iPMD)
  {
        //=================== PMD parameters ============================

    AliPMD *PMD = new AliPMDv1("PMD", "normal PMD");
  }

  if (iT0)
  {
        //=================== T0 parameters ============================
    AliT0 *T0 = new AliT0v1("T0", "T0 Detector");
  }

  if (iEMCAL)
  {
        //=================== EMCAL parameters ============================

    AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_COMPLETEV1");
  }

  if (iACORDE)
  {
        //=================== ACORDE parameters ============================

    AliACORDE *ACORDE = new AliACORDEv1("ACORDE", "normal ACORDE");
  }

  if (iVZERO)
  {
        //=================== ACORDE parameters ============================

    AliVZERO *VZERO = new AliVZEROv7("VZERO", "normal VZERO");
  }
}

