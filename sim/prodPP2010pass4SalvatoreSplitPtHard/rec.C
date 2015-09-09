void rec()
{

  AliReconstruction reco;

  reco.SetRunReconstruction("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");

  //
  // switch off cleanESD, write ESDfriends and Alignment data
  
  reco.SetCleanESD(kFALSE);
  reco.SetWriteESDfriend();
  reco.SetWriteAlignmentData();
  reco.SetFractionFriends(.1);

  //
  // ITS Efficiency and tracking errors

  reco.SetRunPlaneEff(kTRUE);
  reco.SetUseTrackingErrorsForAlignment("ITS");

  //
  // Raw OCDB
  reco.SetDefaultStorage("alien://folder=/alice/data/2010/OCDB");

  //
  // ITS (2 objects)
  reco.SetSpecificStorage("ITS/Align/Data",     	"alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
  reco.SetSpecificStorage("ITS/Calib/SPDSparseDead", 	"alien://folder=/alice/simulation/2008/v4-15-Release/Residual"); 

  //
  // MUON objects (1 object)
  reco.SetSpecificStorage("MUON/Align/Data",		"alien://folder=/alice/simulation/2008/v4-15-Release/Residual"); 
  
  //
  // TPC (7 objects)
  reco.SetSpecificStorage("TPC/Align/Data", 		"alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
  reco.SetSpecificStorage("TPC/Calib/ClusterParam", 	"alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
  reco.SetSpecificStorage("TPC/Calib/RecoParam", 	"alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
  reco.SetSpecificStorage("TPC/Calib/TimeGain", 	"alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
  reco.SetSpecificStorage("TPC/Calib/AltroConfig", 	"alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
  reco.SetSpecificStorage("TPC/Calib/TimeDrift", 	"alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");
  reco.SetSpecificStorage("TPC/Calib/Correction", 	"alien://Folder=/alice/simulation/2008/v4-15-Release/Residual/");

  reco.SetRunQA(":");

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}

