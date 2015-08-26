// root macro to use the class DJetCorrUnfold

class DJetCorrUnfold;

#include "runDJetCorrAnalysis.C"
#include "runDJetCorrResponse.C"

DJetCorrUnfold* runDJetCorrUnfold(const char* trainData = "Jets_EMC_pp_MC_509_510_511_512", const char* trainResp = "Jets_EMC_pp_MC_508",
				  Bool_t loadLibs = kTRUE, const char* inputPath = "$JETRESULTS")
{
  if (loadLibs) {    
    gROOT->LoadMacro("HistoStyler.cxx+g");
    gROOT->LoadMacro("MassFitter.cxx+g");
    gROOT->LoadMacro("DJetCorrAnalysisParams.cxx+g");
    gROOT->LoadMacro("DJetCorrBase.cxx+g");
    gROOT->LoadMacro("DJetCorrAnalysis.cxx+g");
    gROOT->LoadMacro("DJetCorrResponse.cxx+g");
    gROOT->LoadMacro("DJetCorrUnfold.cxx+g");
  }
  
  DJetCorrAnalysis* ana = runDJetCorrAnalysis("", trainData, kTRUE, kFALSE, inputPath);
  DJetCorrResponse* resp = runDJetCorrResponse("", trainResp, kFALSE, inputPath);
  
  DJetCorrUnfold* unfold = new DJetCorrUnfold(ana, resp);
  unfold->SetDataParamIndex(1);
  unfold->SetRespParamIndex(1);
  unfold->Start();

  return unfold;
}
