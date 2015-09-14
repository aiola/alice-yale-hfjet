void sim(Int_t nev=500000)
{
  gSystem->Load("libpythia6_4_25");
  AliSimulation simulator;

  //simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  //simulator.SetMakeDigitsFromHits("ITS TPC");
  simulator.SetRunGeneratorOnly(kTRUE);

  simulator.SetMakeSDigits("");
  simulator.SetMakeDigits("");
  simulator.SetMakeDigitsFromHits("");
  //simulator.SetRunGeneration(kTRUE);   // generation of primary particles
  //simulator.SetRunSimulation(kFALSE);
  /*
  //
  // Raw OCDB
  simulator.SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB");

  //
  // ITS  (1 Total)
  //     Alignment from Ideal OCDB 
  simulator.SetSpecificStorage("ITS/Align/Data",  "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal");
  
  //
  // MUON (1 object)
  simulator.SetSpecificStorage("MUON/Align/Data","alien://folder=/alice/simulation/2008/v4-15-Release/Ideal"); 

  //
  // TPC (7 total) 
  simulator.SetSpecificStorage("TPC/Calib/TimeGain",       "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
  simulator.SetSpecificStorage("TPC/Calib/ClusterParam",   "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
  simulator.SetSpecificStorage("TPC/Calib/AltroConfig",    "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
  simulator.SetSpecificStorage("TPC/Calib/Correction",     "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
  simulator.SetSpecificStorage("TPC/Align/Data",           "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
  simulator.SetSpecificStorage("TPC/Calib/TimeDrift",      "alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
  simulator.SetSpecificStorage("TPC/Calib/RecoParam", 	   "alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
  
  //
  // Vertex and Mag.field from OCDB
  simulator.UseVertexFromCDB();
  simulator.UseMagFieldFromGRP();
  */
  simulator.SetRunQA(":");
  
  //
  // The rest
  
  TStopwatch timer;
  timer.Start();
  Printf("starting...");
  simulator.Run(nev);
  WriteXsection();
  timer.Stop();
  timer.Print();
}

WriteXsection()
{
 TPythia6 *pythia = TPythia6::Instance();
 pythia->Pystat(1);
 Double_t xsection = pythia->GetPARI(1);
 UInt_t    ntrials  = pythia->GetMSTI(5);

 TFile *file = new TFile("pyxsec.root","recreate");
 TTree   *tree   = new TTree("Xsection","Pythia cross section");
 TBranch *branch = tree->Branch("xsection", &xsection, "X/D");
 TBranch *branch = tree->Branch("ntrials" , &ntrials , "X/i");
 tree->Fill();



 tree->Write();
 file->Close();

 cout << "Pythia cross section: " << xsection
      << ", number of trials: " << ntrials << endl;
}

