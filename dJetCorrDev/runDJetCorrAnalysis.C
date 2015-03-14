// root macro to use the class DJetCorrAnalysis

class DJetCorrAnalysis;

void runDJetCorrAnalysis()
{
  gROOT->LoadMacro("DJetCorrAnalysis.cxx+g");

  DJetCorrAnalysis* projDjet = new DJetCorrAnalysis("Jets_EMC_pp_223");
  projDjet->SetOverwrite(kTRUE);

  projDjet->SetUseTest(kFALSE);
  projDjet->SetInputPath("$JETRESULTS");
  projDjet->SetOutputFileName("<train>/DJetCorr.root");

  projDjet->GenerateQAHistograms();
  projDjet->GenerateDJetCorrHistograms("D0", "Charged", "R040");
  projDjet->GenerateDJetCorrHistograms("DStar", "Charged", "R040");

  projDjet->SaveOutputFile();
}
