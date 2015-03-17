// root macro to use the class DJetCorrAnalysis

class DJetCorrAnalysis;

void Run(DJetCorrAnalysis* projDjet);
void Plot(DJetCorrAnalysis* projDjet);

void runDJetCorrAnalysis(const char* options = "run", const char* train = "Jets_EMC_pp_225_226", const char* inputPath = "$JETRESULTS")
{
  gROOT->LoadMacro("DJetCorrAnalysis.cxx+g");

  DJetCorrAnalysis* projDjet = new DJetCorrAnalysis(train);
  
  projDjet->SetOverwrite(kTRUE);
  projDjet->SetInputPath(inputPath);
  projDjet->SetOutputFileName("<train>/DJetCorr.root");

  TString opt(options);
  TObjArray *optList = opt.Tokenize(" ");

  if (optList->Contains("run")) {
    Run(projDjet);
  }

  if (optList->Contains("plot")) {
    Plot(projDjet);
  }
}

void Run(DJetCorrAnalysis* projDjet)
{
  projDjet->GenerateQAHistograms();
  projDjet->GenerateDJetCorrHistograms("D0", "Charged", "R040");
  projDjet->GenerateDJetCorrHistograms("DStar", "Charged", "R040");

  projDjet->SaveOutputFile();
}

void Plot(DJetCorrAnalysis* projDjet)
{
  projDjet->SetPlotFormat("pdf");
  projDjet->SetSavePlots(kFALSE);
  
  //projDjet->PlotTrackHistograms();
  projDjet->PlotDJetCorrHistograms("D0", "Charged", "R040");
  projDjet->PlotDJetCorrHistograms("DStar", "Charged", "R040");
}
