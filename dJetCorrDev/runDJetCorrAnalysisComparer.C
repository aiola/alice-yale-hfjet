// root macro to use the class DJetCorrAnalysisComparer

class DJetCorrAnalysisComparer;

#include "runDJetCorrAnalysis.C"
#include "runDJetCorrResponse.C"

DJetCorrAnalysisComparer* runDJetCorrAnalysisComparer(const char* train1 = "Jets_EMC_pp_MC_645_646_647_648",
                                                      const char* train2 = "Jets_EMC_pp_MC_631_632_633_634",
                                                      const char* train3 = 0,
                                                      Bool_t loadLibs = kTRUE, const char* inputPath = "$JETRESULTS")
{
  if (loadLibs) {
    gROOT->LoadMacro("HistoStyler.cxx+g");
    gROOT->LoadMacro("MassFitter.cxx+g");
    gROOT->LoadMacro("DJetCorrAnalysisParams.cxx+g");
    gROOT->LoadMacro("DJetCorrBase.cxx+g");
    gROOT->LoadMacro("DJetCorrAnalysis.cxx+g");
    gROOT->LoadMacro("DJetCorrResponse.cxx+g");
    gROOT->LoadMacro("DJetCorrAnalysisComparer.cxx+g");
  }
  
  DJetCorrBase* ana1 = runDJetCorrAnalysis("", train1, kFALSE, kTRUE, kFALSE, inputPath);
  ana1->SetTitle("Signal-only");
  DJetCorrBase* ana2 = runDJetCorrAnalysis("", train2, kFALSE, kTRUE, kFALSE, inputPath);
  ana2->SetTitle("Invariant mass fit");

  DJetCorrBase* ana3 = 0;

  if (train3) {
    ana3 = runDJetCorrAnalysis("", train3, kFALSE, kTRUE, kFALSE, inputPath);
    ana3->SetTitle("");
  }
  
  DJetCorrAnalysisComparer* comparer = new DJetCorrAnalysisComparer(DJetCorrAnalysisComparer::kCompareTruth |
                                                                    DJetCorrAnalysisComparer::kCompareMeasured,
                                                                    ana1, ana2, ana3, 0, 0, 0);
  comparer->SetSavePlots(kTRUE);
  comparer->SetOverwrite(kTRUE);
  comparer->SetMakeRatios(kTRUE);
  comparer->SetNormalizationType(DJetCorrAnalysisComparer::kEvents);
  comparer->Start();
  comparer->SaveOutputFile();

  return comparer;
}
