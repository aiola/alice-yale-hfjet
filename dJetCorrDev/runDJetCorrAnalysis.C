// root macro to use the class DJetCorrAnalysis

class DJetCorrAnalysis;

void runDJetCorrAnalysis(const char* options = "run plot", const char* train = "Jets_EMC_pp_225_226_228", const char* inputPath = "$JETRESULTS")
{
  gROOT->LoadMacro("DJetCorrAnalysis.cxx+g");

  DJetCorrAnalysis* projDjet = new DJetCorrAnalysis(train);
  
  projDjet->SetOverwrite(kTRUE);
  projDjet->SetInputPath(inputPath);
  projDjet->SetOutputFileName("<train>/DJetCorr.root");

  projDjet->SetPlotFormat("pdf");
  projDjet->SetSavePlots(kFALSE);

  projDjet->AddAnalysisParams("D0", "Charged", "R040");
  projDjet->AddAnalysisParams("DStar", "Charged", "R040");

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


