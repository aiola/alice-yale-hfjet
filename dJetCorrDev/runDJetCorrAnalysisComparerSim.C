// root macro to use the class DJetCorrAnalysisComparer

class DJetCorrAnalysisComparer;

#include "runDJetCorrAnalysis.C"
#include "runDJetCorrResponse.C"

DJetCorrAnalysisComparer* runDJetCorrAnalysisComparer(const char* train1 = "kPyMbDefault",
                                                      const char* train2 = "kPyJets",
                                                      const char* train3 = "kPyCharmppMNRwmi",
                                                      Bool_t loadLibs = kTRUE, const char* inputPath = "/Users/sa639/Documents/Work/ALICE/alice-yale-hfjet/sim/prodPP2010pass4SalvatoreSplitPtHard/results")
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
  ana1->SetTitle("kPyMbDefault");
  DJetCorrBase* ana2 = runDJetCorrAnalysis("", train2, kFALSE, kTRUE, kFALSE, inputPath);
  ana2->SetTitle("kPyJets");

  DJetCorrBase* ana3 = 0;

  if (train3) {
    ana3 = runDJetCorrAnalysis("", train3, kFALSE, kTRUE, kFALSE, inputPath);
    ana3->SetTitle("kPyCharmppMNRwmi");
  }
  
  DJetCorrAnalysisComparer* comparer = new DJetCorrAnalysisComparer(DJetCorrAnalysisComparer::kCompareTruth,
                                                                    ana1, ana2, ana3, 0, 0, 0);
  comparer->SetSavePlots(kTRUE);
  comparer->SetOverwrite(kTRUE);
  comparer->SetMakeRatios(kTRUE);
  comparer->SetNormalizationType(DJetCorrAnalysisComparer::kEvents);
  comparer->Start();
  comparer->SaveOutputFile();

  return comparer;
}
