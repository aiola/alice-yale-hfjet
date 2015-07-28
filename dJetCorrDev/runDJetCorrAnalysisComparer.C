// root macro to use the class DJetCorrAnalysisComparer

class DJetCorrAnalysisComparer;

#include "runDJetCorrAnalysis.C"

DJetCorrAnalysisComparer* runDJetCorrAnalysisComparer(const char* train1 = "Jets_EMC_pp_MC_502_503_504_505", const char* train2 = "Jets_EMC_pp_407_408_409_410",
                                                      Bool_t loadLibs = kTRUE, const char* inputPath = "$JETRESULTS")
{
  if (loadLibs) {
    gROOT->LoadMacro("HistoStyler.cxx+g");
    gROOT->LoadMacro("MassFitter.cxx+g");
    gROOT->LoadMacro("DJetCorrAnalysisParams.cxx+g");
    gROOT->LoadMacro("DJetCorrBase.cxx+g");
    gROOT->LoadMacro("DJetCorrAnalysis.cxx+g");
    gROOT->LoadMacro("DJetCorrAnalysisComparer.cxx+g");
  }
  
  DJetCorrAnalysis* ana1 = runDJetCorrAnalysis("", train1, kTRUE, kFALSE, inputPath);
  DJetCorrAnalysis* ana2 = runDJetCorrAnalysis("", train2, kFALSE, kFALSE, inputPath);
  
  DJetCorrAnalysisComparer* comparer = new DJetCorrAnalysisComparer(ana1, ana2);
  comparer->Start();
}
