//
// Configuration for PWG3 barrel open charm and beauty 2010
// PYTHIA pp 7000 TeV Perugia0
// 20% ccbar pair per event
//     at least one in |y|<1.5
//     D mesons decay hadronically
// 20% bbbar pair per event
//     at least one in |y|<1.5
//     D mesons decay hadronically
// 20% ccbar pair per event
//     decays not forced
//     at least one electron from charm in |y|<1.2
// 20% bbbar pair per event
//     decays not forced
//     at least one electron from charm or beauty in |y|<1.2
//  10% J/psi(|y|<1.0)->e+e-
//  10% B(|y|<2.0)->J/psi(|y|<2.0)->e+e-
//
// One can use the configuration macro in compiled mode by
// root [0] gSystem->Load("libgeant321");
// root [0] gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include\
//                   -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");
// root [0] .x grun.C(1,"Config.C++")

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TRandom.h>
#include <TDatime.h>
#include <TSystem.h>
#include <TVirtualMC.h>
#include <TGeant3TGeo.h>
#include "STEER/STEER/AliRunLoader.h"
#include "STEER/STEER/AliRun.h"
#include "STEER/STEER/AliConfig.h"
#include "STEER/STEER/AliSimulation.h"
#include "PYTHIA6/AliDecayerPythia.h"
#include "PYTHIA6/AliGenPythia.h"
#include "PYTHIA6/AliGenPythiaPlus.h"
#include "TDPMjet/AliGenDPMjet.h"
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


enum PDC06Proc_t
{
  kPythia6, kPythia6D6T, kPythia6ATLAS, kPythia6ATLAS_Flat, kPythiaPerugia0, kPhojet, kPythiaPerugia0charmbarrel, kPythiaPerugia0beautybarrel, kPythiaPerugia0chadr, kPythiaPerugia0bchadr,kPythiaPerugia2011chadrPtHard, kPythiaPerugia2011bchadrPtHard, kPythiaPerugia0cele, kPythiaPerugia0bele, kPythiaPerugia0Jpsi2e, kPythiaPerugia0BtoJpsi2e, kPythiaPerugia0HFbarrel, kPythiaPerugia2011HFbarrelPtHard, kRunMax,
};

const char * pprRunName[] = {
  "kPythia6", "kPythia6D6T", "kPythia6ATLAS", "kPythia6ATLAS_Flat", "kPythiaPerugia0", "kPhojet",  "kPythiaPerugia0charmbarrel", "kPythiaPerugia0beautybarrel", "kPythiaPerugia0chadr", "kPythiaPerugia0bchadr", "kPythiaPerugia2011chadrPtHard", "kPythiaPerugia2011bchadrPtHard","kPythiaPerugia0cele", "kPythiaPerugia0bele", "kPythiaPerugia0Jpsi2e", "kPythiaPerugia0BtoJpsi2e", "kPythiaPerugia0HFbarrel","kPythiaPerugia2011HFbarrelPtHard"
};

enum Mag_t
{
  kNoField, k5kG, kFieldMax
};

const char * pprField[] = {
  "kNoField", "k5kG"
};

enum PprTrigConf_t
{
    kDefaultPPTrig, kDefaultPbPbTrig
};

const char * pprTrigConfName[] = {
    "p-p","Pb-Pb"
};

//--- Functions ---
class AliGenPythia;

AliGenerator *MbPythiaTunePerugia2011bchadrPtHard();
AliGenerator *MbPythiaTunePerugia2011bchadrPtHard();

void ProcessEnvironmentVars();

// Generator, field, beam energy
static PDC06Proc_t   proc     = kPythiaPerugia2011HFbarrelPtHard;
static Mag_t         mag      = k5kG;
static Float_t       energy   = 7000; // energy in CMS
static Double_t      JpsiPol  = 0; // Jpsi polarisation
static Bool_t        JpsiHarderPt = kFALSE; // Jpsi harder pt spectrum (8.8 TeV)
static Int_t         runNumber = 0; 
static PprTrigConf_t strig = kDefaultPPTrig; // default pp trigger configuration
//========================//
// Set Random Number seed //
//========================//
TDatime dt;
static UInt_t seed    = dt.Get();

// Comment line
static TString comment;

const Int_t npthardbins=8;
Double_t pthardmin=-1,pthardmax=-1;
Double_t pthardbins[npthardbins+1]={5,11,21,36,57,84,117,152,191};
Int_t binpthard=-99;
void Config()
{


  // Get settings from environment variables
  ProcessEnvironmentVars();

  gRandom->SetSeed(seed);
  cerr<<"Seed for random number generation= "<<seed<<endl;

  // Libraries required by geant321
#if defined(__CINT__)
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  if (proc == kPythia6 || proc == kPhojet) {
    Printf("Loading PYTHIA 6.2");
    gSystem->Load("libpythia6");        // Pythia 6.2
  } else {
    Printf("Loading PYTHIA 6.4.25");
    gSystem->Load("libpythia6_4_25");   // Pythia 6.4
  }
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libgeant321");

#endif

  new TGeant3TGeo("C++ Interface to Geant3");

  //=======================================================================
  //  Create the output file


  AliRunLoader* rl=0x0;

  cout<<"Config.C: Creating Run Loader ..."<<endl;
  rl = AliRunLoader::Open("galice.root",
			  AliConfig::GetDefaultEventFolderName(),
			  "recreate");
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

  // Set the trigger configuration
  AliSimulation::Instance()->SetTriggerConfig(pprTrigConfName[strig]);
  cout<<"Trigger configuration is set to  "<<pprTrigConfName[strig]<<endl;

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

    // RANDOM SELECTION OF ONE OF THE SIX GENERATION TYPES
    //
    Int_t typeHF = -1;
    Float_t randHF = gRandom->Rndm();
    if(randHF<0.5) {
      typeHF=0;
    } else {
      typeHF=1;
    }
    
    // random selection of pt hard bin
    if(binpthard==-99){
      Float_t randPtHard=gRandom->Rndm();
      binpthard=(Int_t)(randPtHard*(npthardbins+1));// the +1 comes from the fact that we want a sample w/o pthard bins
    }
    if(binpthard){
      pthardmin=pthardbins[binpthard-1];// the -1 comes as a consequence of the +1 above
      pthardmax=pthardbins[binpthard];
    }

  //======================//
  // Set External decayer //
  //======================//
   
  TVirtualMCDecayer* decayer = new AliDecayerPythia();
  if(typeHF==0 || typeHF==1) {// is this really needed??
    decayer->SetForceDecay(kHadronicDWithout4BodiesWithV0);
  } else {
    decayer->SetForceDecay(kAll);
  }

  decayer->Init();
  gMC->SetExternalDecayer(decayer);

   

  //=========================//
  // Generator Configuration //
  //=========================//
  AliGenerator* gener = 0x0;

gener = MbPythiaTunePerugia2011chadrPtHard();
/*

  if (proc == kPythiaPerugia2011chadrPtHard ||
      (proc == kPythiaPerugia2011HFbarrelPtHard && typeHF==0)) {
    gener = MbPythiaTunePerugia2011chadrPtHard();
  } else if (proc == kPythiaPerugia2011bchadrPtHard ||
	     (proc == kPythiaPerugia2011HFbarrelPtHard && typeHF==1)) {
    gener = MbPythiaTunePerugia2011bchadrPtHard();
  }
*/
  //
  //
  // Size of the interaction diamond
  // Longitudinal
  Float_t sigmaz  = 5.4 / TMath::Sqrt(2.); // [cm]
  if (energy == 900)
    //sigmaz  = 10.5 / TMath::Sqrt(2.); // [cm]
    //sigmaz = 3.7;
  if (energy == 7000)
    sigmaz  = 6.3 / TMath::Sqrt(2.); // [cm]

  //
  // Transverse
  Float_t betast  = 10;                 // beta* [m]
  if (runNumber >= 117048) betast = 2.;
  if (runNumber >  122375) betast = 3.5; // starting with fill 1179

  Float_t eps     = 5.e-6;            // emittance [m]
  Float_t gamma   = energy / 2.0 / 0.938272;  // relativistic gamma [1]
  Float_t sigmaxy = TMath::Sqrt(eps * betast / gamma) / TMath::Sqrt(2.) * 100.;  // [cm]
  printf("\n \n Diamond size x-y: %10.3e z: %10.3e\n \n", sigmaxy, sigmaz);

  gener->SetOrigin(0., 0., 0.); // Taken from OCDB
  gener->SetSigma(sigmaxy, sigmaxy, sigmaz);      // Sigma in (X,Y,Z) (cm) on IP position
  gener->SetVertexSmear(kPerEvent);
  gener->Init();

  printf("\n \n Comment: %s \n \n", comment.Data());

    //	
   // FIELD
   //


 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG,
     	   	AliMagF::kBeamTypepp));

  rl->CdGAFile();

  
  /*
  Int_t iABSO  = 1;
  Int_t iACORDE= 0;
  Int_t iDIPO  = 1;
  Int_t iEMCAL = 1;
  Int_t iFMD   = 1;
  Int_t iFRAME = 1;
  Int_t iHALL  = 1;
  Int_t iITS   = 1;
  Int_t iMAG   = 1;
  Int_t iMUON  = 1;
  Int_t iPHOS  = 1;
  Int_t iPIPE  = 1;
  Int_t iPMD   = 1;
  Int_t iHMPID = 1;
  Int_t iSHIL  = 1;
  Int_t iT0    = 1;
  Int_t iTOF   = 1;
  Int_t iTPC   = 1;
  Int_t iTRD   = 1;
  Int_t iVZERO = 1;
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
	//	AliITS *ITS  = new AliITSv11Hybrid("ITS","ITS v11Hybrid");
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
	// starting at 3h in positive direction
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
	MUON->SetTriggerResponseV1(2);
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

        AliEMCAL *EMCAL = new AliEMCALv2("EMCAL", "EMCAL_FIRSTYEARV1");
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
  */
}

//
//           PYTHIA GENERATORS
//


AliGenerator* MbPythiaTunePerugia2011chadrPtHard()
{
      comment = comment.Append(" pp: cocktail with single generator, Pythia (Perugia2011) chadr (1 ccbar per event, 1 c-quark in |y|<1.5, chadrons decay to hadrons), pt hard bin");
      AliGenCocktail *cocktail = new AliGenCocktail();
      cocktail->SetProjectile("p", 1, 1);
      cocktail->SetTarget    ("p", 1, 1);
      cocktail->SetEnergyCMS(energy);
//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1);
      pythia->SetTitle("PYTHIA-HF-chadr");
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);
      pythia->SetEnergyCMS(energy);

      //    Tune
      //    350     Perugia 2011
      //    320     Perugia 0
      pythia->SetTune(320);
      pythia->UseNewMultipleInteractionsScenario();

      //pythia->SetProcess(kPyMbDefault);
      //pythia->SetProcess(kPyMb);
      pythia->SetProcess(kPyCharmppMNRwmi);
      //pythia->SetProcess(kPyCharm);
      //pythia->SetProcess(kPyJets);
      pythia->SetHeavyQuarkYRange(-1.5,1.5);
      pythia->SetTriggerParticle(4, 5., -1., 1000);      
    
      
//
//    decays
     pythia->SetForceDecay(kHadronicDWithout4BodiesWithV0);

      // set pt hard and name
      TString nameHF="pythiaPer2011_ccbarHFhadr";
      if(pthardmax>0){
	pythia->SetPtHard(pthardmin,pthardmax);
	nameHF.Append(Form("_%.0f",pthardmin));
      }
      else {
	nameHF.Append("_NoPtHard");
      }
      cocktail->AddGenerator(pythia,nameHF.Data(),1);



      return cocktail;
}
/*
AliGenerator* MbPythiaTunePerugia2011bchadrPtHard()
{
      comment = comment.Append(" pp: cocktail with single generator, Pythia (Perugia2011) cchadr (1 bbbar per event, 1 b-quark in |y|<1.5, chadrons decay to hadrons), pt hard bin");
      AliGenCocktail *cocktail = new AliGenCocktail();
      cocktail->SetProjectile("p", 1, 1);
      cocktail->SetTarget("p", 1, 1);
      cocktail->SetEnergyCMS(energy);


//
//    Pythia
      AliGenPythia* pythia = new AliGenPythia(-1);
      pythia->SetTitle("PYTHIA-HF-bchadr");
      pythia->SetMomentumRange(0, 999999.);
      pythia->SetThetaRange(0., 180.);
      pythia->SetYRange(-12.,12.);
      pythia->SetPtRange(0,1000.);
      pythia->SetEnergyCMS(energy);

//    Tune
//    350     Perugia 2011
      pythia->SetTune(350);
      pythia->UseNewMultipleInteractionsScenario();
      
      pythia->SetProcess(kPyMb);
      //pythia->SetProcess(kPyBeautyppMNRwmi);
      //pythia->SetHeavyQuarkYRange(-1.5,1.5);
      //pythia->SetTriggerParticle(4, 5., -1., 1000);      

//
//    decays
      pythia->SetForceDecay(kHadronicDWithout4BodiesWithV0);

      // set pt hard and name
      TString nameHF="pythiaPer2011_bbbarHFhadr";
      if(pthardmax>0){
	pythia->SetPtHard(pthardmin,pthardmax);
	nameHF.Append(Form("_%.0f",pthardmin));
      }
      else {
	nameHF.Append("_NoPtHard");
      }
      cocktail->AddGenerator(pythia,nameHF.Data(),1);



      return cocktail;

}


*/

void ProcessEnvironmentVars()
{
    // Run type
    if (gSystem->Getenv("CONFIG_RUN_TYPE")) {
      for (Int_t iRun = 0; iRun < kRunMax; iRun++) {
	if (strcmp(gSystem->Getenv("CONFIG_RUN_TYPE"), pprRunName[iRun])==0) {
	  proc = (PDC06Proc_t)iRun;
	  cout<<"Run type set to "<<pprRunName[iRun]<<endl;
	}
      }
    }

    // Field
    if (gSystem->Getenv("CONFIG_FIELD")) {
      for (Int_t iField = 0; iField < kFieldMax; iField++) {
	if (strcmp(gSystem->Getenv("CONFIG_FIELD"), pprField[iField])==0) {
	  mag = (Mag_t)iField;
	  cout<<"Field set to "<<pprField[iField]<<endl;
	}
      }
    }

    // Energy
    if (gSystem->Getenv("CONFIG_ENERGY")) {
      energy = atoi(gSystem->Getenv("CONFIG_ENERGY"));
      cout<<"Energy set to "<<energy<<" GeV"<<endl;
    }

    // Random Number seed
    if (gSystem->Getenv("CONFIG_SEED")) {
      seed = atoi(gSystem->Getenv("CONFIG_SEED"));
    }
    
    if(gSystem->Getenv("CONFIG_PTHARDBIN")){
      binpthard = atoi(gSystem->Getenv("CONFIG_PTHARDBIN"));
    }
    
    // Run number
    if (gSystem->Getenv("DC_RUN")) {
      runNumber = atoi(gSystem->Getenv("DC_RUN"));
    }
}
