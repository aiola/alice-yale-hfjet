// root macro to use the class DJetCorrAnalysis

class DJetCorrAnalysis;

DJetCorrAnalysis* runDJetCorrAnalysis(const char* options = "refit plot", const char* train = "Jets_EMC_pp_MC_618_619_620_621",
                                      Bool_t loadLibs = kTRUE, Bool_t isMC = kTRUE, Bool_t isBkgSub = kTRUE,
                                      const char* inputPath = "$JETRESULTS")
{
  TGaxis::SetMaxDigits(3); 

  TString tracksName = "tracks";
  //TString tracksD0Name = tracksName;
  //TString tracksDStarName = tracksName;
  TString tracksD0Name = "DSBcandidatesAndTracksD0";
  TString tracksDStarName = "DSBcandidatesAndTracksDStar";

  if (!isMC) {
    tracksD0Name += "rec";
    tracksDStarName += "rec";
  }

  if (loadLibs) {
    gROOT->LoadMacro("HistoStyler.cxx+g");
    gROOT->LoadMacro("MassFitter.cxx+g");
    gROOT->LoadMacro("DJetCorrAnalysisParams.cxx+g");
    gROOT->LoadMacro("DJetCorrBase.cxx+g");
    gROOT->LoadMacro("DJetCorrAnalysis.cxx+g");
  }
  
  DJetCorrAnalysis* projDjet = new DJetCorrAnalysis(train);
  
  projDjet->SetOverwrite(kTRUE);
  projDjet->SetInputPath(inputPath);

  TString qaListName(Form("AliAnalysisTaskSAQA_%s_TPC_histos", tracksName.Data()));
  
  projDjet->SetQAListName(qaListName);

  projDjet->SetPlotFormat("pdf");
  projDjet->SetSavePlots(kTRUE);

  DJetCorrAnalysisParams* param = 0;

  //projDjet->AddAnalysisParams("D0", "Full", "R040", tracksD0Name, isMC, isBkgSub);
  param = projDjet->AddAnalysisParams("D0", "Full", "R060", tracksD0Name, isMC, isBkgSub);
  param->SetBackgroundOnly(kTRUE);
  //projDjet->AddAnalysisParams("DStar", "Charged", "R040", tracksDStarName, isMC, isBkgSub);
  //param = projDjet->AddAnalysisParams("DStar", "Charged", "R060", tracksDStarName, isMC, isBkgSub);
  //param->SetInvMassRebinFactor(2);

  TString opt(options);
  TObjArray *optList = opt.Tokenize(" ");

  if (optList->Contains("truth-only")) {
    projDjet->ProjectTruthSpectrum();
  }
  else {
    if (optList->Contains("run")) {
      projDjet->GenerateQAHistograms();
      projDjet->GenerateDJetCorrHistograms();
      projDjet->ProjectTruthSpectrum();
    }

    if (optList->Contains("refit")) {
      projDjet->PlotTrackHistograms();
      projDjet->PlotDJetCorrHistograms(kTRUE);
    }
    else if (optList->Contains("plot")) {
      projDjet->PlotTrackHistograms();
      projDjet->PlotDJetCorrHistograms(kFALSE);
    }
  }

  if (!opt.IsNull()) projDjet->SaveOutputFile();

  return projDjet;
}


