void simSimpleGrid(Int_t nev=1, char *configfile="ConfigSimpleGrid.C") 
{
 gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/CDB -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/AOD -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_PHYSICS/../src/OADB -I$ALICE_ROOT/PWGHF -I$ALICE_PHYSICS/../src/PWGHF/base -I$ALICE_PHYSICS/../src/PWGHF/correlationHF -I$ALICE_PHYSICS/PWGHF/correlationHF/macros -I$ALICE_PHYSICS/../src/PWG/FLOW/Base -I$ALICE_PHYSICS../src/PWG/FLOW/Tasks -I$ALICE_ROOT/PWGJE -I$ALICE_ROOT/JETAN -I$ALICE_PHYSICS/../src/CORRFW -I$ALICE_ROOT/PYTHIA6 -g");
  gSystem->Load("liblhapdf");
  gSystem->Load("libpythia6_4_25.so");
  gSystem->Load("libEGPythia6.so");
  gSystem->Load("libAliPythia6.so");

  
  cout<<"==================  simulating events:"<<nev<<endl;
  AliSimulation simulator(configfile);
  simulator.SetRunNumber(156797); // set anchor run: good run from Tomas Aronsson analysis note
  //simulator.SetMakeSDigits("TRD TOF PHOS EMCAL T0 VZERO");
  //simulator.SetMakeDigitsFromHits("ITS TPC");
  simulator.SetRunGeneration(kTRUE);
  //  simulator.SetMakeDigits("TRD TOF PHOS EMCAL T0 VZERO");
//
// OCDB
//
  simulator.SetDefaultStorage("alien://Folder=/alice/data/2011/OCDB");
  // simulator.SetDefaultStorage("local://./OCDB");
  
//   simulator.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  
//   simulator.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
//
// TRD
//
//   simulator.SetSpecificStorage("TRD/Calib/DCS","alien://Folder=/alice/data/2011/OCDB");
//   simulator.SetSpecificStorage("TRD/Calib/ChamberStatus","alien://Folder=/alice/data/2011/OCDB");
//   simulator.SetSpecificStorage("TRD/Calib/ChamberVdrift","alien://Folder=/alice/data/2011/OCDB");
//   simulator.SetSpecificStorage("TRD/Calib/ChamberT0","alien://Folder=/alice/data/2011/OCDB");
//   simulator.SetSpecificStorage("TRD/Calib/ChamberExB","alien://Folder=/alice/data/2011/OCDB");
//
// Vertex and magfield
//
  //simulator.SetSpecificStorage("GRP/GRP/Data","alien://Folder=/alice/data/2011/OCDB");
  simulator.UseVertexFromCDB();
  simulator.UseMagFieldFromGRP();
//
// PHOS simulation settings
//
  //AliPHOSSimParam *simParam = AliPHOSSimParam::GetInstance();
  //simParam->SetCellNonLineairyA(0.001);
  //simParam->SetCellNonLineairyB(0.2);
  //simParam->SetCellNonLineairyC(1.02);

  simulator.SetRunQA(":");

  simulator.Run(nev);
}
