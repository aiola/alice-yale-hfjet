// root macro to use the class DJetCorrAnalysis

class DJetCorrAnalysis;

void runDJetCorrAnalysis(const char* options = "run plot", const char* train = "Jets_EMC_pp_243_244_245_246", const char* inputPath = "$JETRESULTS")
{
  const char* tracksName = "tracks";
  
  gROOT->LoadMacro("DJetCorrAnalysisParams.cxx+g");
  gROOT->LoadMacro("DJetCorrAnalysis.cxx+g");

  DJetCorrAnalysis* projDjet = new DJetCorrAnalysis(train);
  
  projDjet->SetOverwrite(kTRUE);
  projDjet->SetInputPath(inputPath);
  projDjet->SetOutputFileName("DJetCorr.root");

  TString qaListName(Form("AliAnalysisTaskSAQA_%s_TPC_histos", tracksName));
  
  projDjet->SetQAListName(qaListName);

  projDjet->SetPlotFormat("pdf");
  projDjet->SetSavePlots(kTRUE);

  projDjet->AddAnalysisParams("D0", "Charged", "R040", tracksName);
  projDjet->AddAnalysisParams("DStar", "Charged", "R040", tracksName);
  projDjet->AddAnalysisParams("D0", "Charged", "R060", tracksName);
  projDjet->AddAnalysisParams("DStar", "Charged", "R060", tracksName);

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


