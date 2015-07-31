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

  Bool_t ProjectResponseDPtMatrix(DJetCorrAnalysisParams* params);

  Bool_t ProjectResponseMatrices4DVsDPt(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseZMatricesVsDPt(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseZMatricesVsJetPt(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseJetPtMatricesVsDPt(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseJetPtMatricesVsZ(DJetCorrAnalysisParams* params);
  
  Bool_t ProjectResponseMatrix4D(DJetCorrAnalysisParams* params, Int_t dptBin=-1);
  Bool_t ProjectResponseJetPtMatrix(DJetCorrAnalysisParams* params, Int_t zBin, Int_t dptBin=-1);
  Bool_t ProjectResponseZMatrix(DJetCorrAnalysisParams* params, Int_t jetPtBin, Int_t dptBin=-1);

  Bool_t PlotResponseDPtMatrix(DJetCorrAnalysisParams* params);
    
  Bool_t PlotResponseMatrices4DVsDPt(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseZMatricesVsDPt(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseZMatricesVsJetPt(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseJetPtMatricesVsDPt(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseJetPtMatricesVsZ(DJetCorrAnalysisParams* params);

  Bool_t PlotResponseDPtMatrix(TCanvas*& canvasResp, TCanvas*& canvasEff, DJetCorrAnalysisParams* params);
  Bool_t PlotResponseMatrix4D(TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t dptBin=-1);
  Bool_t PlotResponseJetPtMatrix(TCanvas*& canvasResp, TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t zBin, Int_t dptBin=-1);
  Bool_t PlotResponseZMatrix(TCanvas*& canvasResp, TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t jetPtBin, Int_t dptBin=-1);

  TString         fEfficiencyMode      ;//  Efficiency mode
  
  THnSparse      *fHistMatching        ;//! THnSparse that contains the matched jets
  THnSparse      *fHistJets1           ;//! THnSparse that contains all reconstructed jets
  THnSparse      *fHistJets2           ;//! THnSparse that contains all generated jets
  
 private:
  DJetCorrResponse(const DJetCorrResponse &source);
  DJetCorrResponse& operator=(const DJetCorrResponse& source); 

  ClassDef(DJetCorrResponse, 1);
};
