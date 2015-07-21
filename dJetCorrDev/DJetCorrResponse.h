// This class projects the THnSparse results of the AliJetResponseMaker into TH1/TH2/TH3
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

class TString;
class TFile;
class TList;
class THnSparse;
class TCanvas;
class TLegend;
class TVirtualPad;
class TPaveText;
class DJetCorrAnalysisParams;
class TObject;
class TList;
class TDirectoryFile;
class TH1;

#include "DJetCorrBase.h"

class DJetCorrResponse : public DJetCorrBase {
 public:
  DJetCorrResponse();
  DJetCorrResponse(const char* train, const char* path = "$JETRESULTS");

  Bool_t ProjectResponseMatrices(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseMatrices(const char* paramName);
  Bool_t ProjectResponseMatrices(Int_t i);
  Bool_t ProjectResponseMatrices();

  Bool_t PlotResponseMatrices(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseMatrices(const char* paramName);
  Bool_t PlotResponseMatrices(Int_t i);
  Bool_t PlotResponseMatrices();
  
 protected:
  Bool_t ClearInputData();
  Bool_t LoadTHnSparse();

  Bool_t ProjectResponseMatricesVsJetPt(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseMatricesVsZ(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseMatrix(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseDPtMatrix(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseJetPtMatrix(DJetCorrAnalysisParams* params, Double_t minZ, Double_t maxZ);
  Bool_t ProjectResponseZMatrix(DJetCorrAnalysisParams* params, Double_t minJetPt, Double_t maxJetPt);

  Bool_t PlotResponseMatricesVsJetPt(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseMatricesVsZ(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseMatrix(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseJetPtMatrix(DJetCorrAnalysisParams* params, Double_t minZ, Double_t maxZ);
  Bool_t PlotResponseZMatrix(DJetCorrAnalysisParams* params, Double_t minJetPt, Double_t maxJetPt);
  
  THnSparse      *fHistMatching        ;//! THnSparse that contains the matched jets
  THnSparse      *fHistJets1           ;//! THnSparse that contains all reconstructed jets
  THnSparse      *fHistJets2           ;//! THnSparse that contains all generated jets
  
 private:
   
  DJetCorrResponse(const DJetCorrResponse &source);
  DJetCorrResponse& operator=(const DJetCorrResponse& source); 

  ClassDef(DJetCorrResponse, 1);
};
