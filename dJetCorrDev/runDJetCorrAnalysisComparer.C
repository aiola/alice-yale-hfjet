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
  
  DJetCorrBase* ana  = runDJetCorrAnalysis("", train1, kTRUE, kFALSE, inputPath);
  DJetCorrBase* resp = runDJetCorrAnalysis("", train2, kTRUE, kFALSE, inputPath);
  
  DJetCorrAnalysisComparer* comparer = new DJetCorrAnalysisComparer(DJetCorrAnalysisComparer::kCompareTruth, ana, resp, 0, 0);
  comparer->SetMakeRatios(kFALSE);
  comparer->Start();
}
