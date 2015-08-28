// root macro to use the class DJetCorrAnalysisComparer

class DJetCorrAnalysisComparer;

#include "runDJetCorrAnalysis.C"
#include "runDJetCorrResponse.C"

DJetCorrAnalysisComparer* runDJetCorrAnalysisComparer(const char* train1 = "Jets_EMC_pp_MC_509_510_511_512", const char* train2 = "Jets_EMC_pp_MC_508",
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
  DJetCorrBase* ana2 = runDJetCorrAnalysis("", train2, kFALSE, kTRUE, kTRUE, inputPath);
  
  DJetCorrAnalysisComparer* comparer = new DJetCorrAnalysisComparer(DJetCorrAnalysisComparer::kCompareTruth, ana1, ana2, 0, 0);
  comparer->SetSavePlots(kTRUE);
  comparer->SetOverwrite(kTRUE);
  comparer->SetMakeRatios(kFALSE);
  comparer->Start();
  comparer->SaveOutputFile();

  return comparer;
}
