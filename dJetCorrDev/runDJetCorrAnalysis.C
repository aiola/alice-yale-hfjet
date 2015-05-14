// root macro to use the class DJetCorrAnalysis

class DJetCorrAnalysis;

void runDJetCorrAnalysis(const char* options = "plot", const char* train = "Jets_EMC_pp_300_301_302_303", const char* inputPath = "$JETRESULTS")
{
  TString tracksName = "tracks";
  //TString tracksD0Name = tracksName;
  //TString tracksDStarName = tracksName;
  TString tracksD0Name = "DcandidatesAndTracksD0rec";
  TString tracksDStarName = "DcandidatesAndTracksDStarrec";
  
  gROOT->LoadMacro("DJetCorrAnalysisParams.cxx+g");
  gROOT->LoadMacro("DJetCorrAnalysis.cxx+g");

  DJetCorrAnalysis* projDjet = new DJetCorrAnalysis(train);
  
  projDjet->SetOverwrite(kTRUE);
  projDjet->SetInputPath(inputPath);
  projDjet->SetOutputFileName("DJetCorr.root");

  TString qaListName(Form("AliAnalysisTaskSAQA_%s_TPC_histos", tracksName.Data()));
  
  projDjet->SetQAListName(qaListName);

  projDjet->SetPlotFormat("pdf");
  projDjet->SetSavePlots(kTRUE);

  projDjet->AddAnalysisParams("D0", "Full", "R040", tracksD0Name);
  projDjet->AddAnalysisParams("D0", "Full", "R060", tracksD0Name);

  projDjet->AddAnalysisParams("DStar", "Charged", "R040", tracksDStarName);
  projDjet->AddAnalysisParams("DStar", "Charged", "R060", tracksDStarName);

  TString opt(options);
  TObjArray *optList = opt.Tokenize(" ");
  
  if (optList->Contains("run")) {
    projDjet->GenerateQAHistograms();
    projDjet->GenerateDJetCorrHistograms();
    
    projDjet->SaveOutputFile();
  }

  if (optList->Contains("plot")) {
    projDjet->PlotTrackHistograms();
    projDjet->PlotDJetCorrHistograms();
  }
}


