// root macro to use the class DJetCorrResponse

class DJetCorrAnalysis;

void runDJetCorrResponse(const char* options = "run plot", const char* train = "Jets_EMC_pp_MC_435", const char* inputPath = "$JETRESULTS")
{
  TGaxis::SetMaxDigits(3); 
  
  TString tracksName = "tracks";
  //TString tracksD0Name = tracksName;
  //TString tracksDStarName = tracksName;
  TString tracksD0Name = "DcandidatesAndTracksD0MCrec";
  TString tracksDStarName = "DcandidatesAndTracksDStarMCrec";

  gROOT->LoadMacro("MassFitter.cxx+g");
  gROOT->LoadMacro("DJetCorrAnalysisParams.cxx+g");
  gROOT->LoadMacro("DJetCorrBase.cxx+g");
  gROOT->LoadMacro("DJetCorrResponse.cxx+g");

  DJetCorrResponse* projDjet = new DJetCorrResponse(train);
  
  projDjet->SetOverwrite(kTRUE);
  projDjet->SetInputPath(inputPath);
  projDjet->SetOutputFileName("DJetCorr.root");

  projDjet->SetPlotFormat("pdf");
  projDjet->SetSavePlots(kTRUE);

  projDjet->AddAnalysisParams("D0", "Full", "R040", tracksD0Name);
  projDjet->AddAnalysisParams("D0", "Full", "R060", tracksD0Name);

  projDjet->AddAnalysisParams("DStar", "Charged", "R040", tracksDStarName);
  projDjet->AddAnalysisParams("DStar", "Charged", "R060", tracksDStarName);

  TString opt(options);
  TObjArray *optList = opt.Tokenize(" ");
  
  if (optList->Contains("run")) {
    projDjet->ProjectResponseMatrices();
    projDjet->SaveOutputFile();
  }

  if (optList->Contains("plot")) {
    projDjet->PlotResponseMatrices();
  }
}


