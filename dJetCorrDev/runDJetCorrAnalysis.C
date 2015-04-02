// root macro to use the class DJetCorrAnalysis

class DJetCorrAnalysis;

void runDJetCorrAnalysis(const char* options = "plot", const char* train = "Jets_EMC_pp_235_236_237_238", const char* inputPath = "$JETRESULTS")
{
  gROOT->LoadMacro("DJetCorrAnalysisParams.cxx+g");
  gROOT->LoadMacro("DJetCorrAnalysis.cxx+g");

  DJetCorrAnalysis* projDjet = new DJetCorrAnalysis(train);
  
  projDjet->SetOverwrite(kTRUE);
  projDjet->SetInputPath(inputPath);
  projDjet->SetOutputFileName("DJetCorr.root");

  projDjet->SetPlotFormat("pdf");
  projDjet->SetSavePlots(kTRUE);

  projDjet->AddAnalysisParams("D0", "Charged", "R040");
  projDjet->AddAnalysisParams("DStar", "Charged", "R040");
  projDjet->AddAnalysisParams("D0", "Charged", "R060");
  projDjet->AddAnalysisParams("DStar", "Charged", "R060");

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


