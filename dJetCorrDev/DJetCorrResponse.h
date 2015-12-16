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

  void SetEfficiencyMode(const char* m) { fEfficiencyMode = m; }

  Bool_t ProjectResponseMatrices(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseMatrices(const char* paramName);
  Bool_t ProjectResponseMatrices(Int_t i);
  Bool_t ProjectResponseMatrices();

  Bool_t PlotResponseMatrices(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseMatrices(const char* paramName);
  Bool_t PlotResponseMatrices(Int_t i);
  Bool_t PlotResponseMatrices();

  Bool_t Regenerate();

  TString GetTruthName(Int_t p);
  TString GetMeasuredName(Int_t p);
  TString GetResponseName(Int_t p);
  TString GetMissesName(Int_t p);
  TString GetKinMissesName(Int_t p);
  TString GetEfficiencyName(Int_t p);
  TString GetKinEfficiencyName(Int_t p);

  THnSparse* GetResponse(Int_t p, Bool_t copy);
  TH2* GetMisses(Int_t p, Bool_t copy);
  TH2* GetKinMisses(Int_t p, Bool_t copy);
  TH2* GetEfficiency(Int_t p, Bool_t copy);
  TH2* GetKinEfficiency(Int_t p, Bool_t copy);
  
 protected:
  Bool_t ClearInputData();
  Bool_t LoadTHnSparse();

  Bool_t ProjectResponseJetPtZVsDPt(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseZVsDPt(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseZVsJetPt(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseJetPtVsDPt(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseJetPtVsZ(DJetCorrAnalysisParams* params);
  Bool_t ProjectResponseDPt(DJetCorrAnalysisParams* params);

  Bool_t ProjectResponseJetPtZ(DJetCorrAnalysisParams* params, Int_t dptBin=-1);
  Bool_t ProjectResponseJetPt(DJetCorrAnalysisParams* params, Int_t zBin, Int_t dptBin=-1);
  Bool_t ProjectResponseZ(DJetCorrAnalysisParams* params, Int_t jetPtBin, Int_t dptBin=-1);

  Bool_t PlotResponseJetPtZVsDPt(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseZVsDPt(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseZVsJetPt(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseJetPtVsDPt(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseJetPtVsZ(DJetCorrAnalysisParams* params);
  Bool_t PlotResponseDPt(DJetCorrAnalysisParams* params);

  Bool_t PlotEfficiencyJetPtZ(TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t dptBin=-1);
  Bool_t PlotEfficiencyJetPt(TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t zBin, Int_t dptBin=-1);
  Bool_t PlotEfficiencyZ(TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t jetPtBin, Int_t dptBin=-1, Int_t jetPtBinStep=1, Int_t dptBinStep=1);
  Bool_t PlotEfficiencyDPt(TCanvas*& canvasEff, DJetCorrAnalysisParams* params);

  Bool_t PlotResponseMatrixJetPt(TCanvas*& canvasResp, DJetCorrAnalysisParams* params, Int_t zBin, Int_t dptBin=-1);
  Bool_t PlotResponseMatrixZ(TCanvas*& canvasResp, DJetCorrAnalysisParams* params, Int_t jetPtBin, Int_t dptBin=-1);
  Bool_t PlotResponseMatrixDPt(TCanvas*& canvasResp, DJetCorrAnalysisParams* params);

  Bool_t PlotResponseMatrixDPt(DJetCorrAnalysisParams* params);

  TString         fEfficiencyMode      ;//  Efficiency mode
  
  THnSparse      *fHistMatching        ;//! THnSparse that contains the matched jets
  THnSparse      *fHistJets1           ;//! THnSparse that contains all reconstructed jets
  THnSparse      *fHistJets2           ;//! THnSparse that contains all generated jets
  
 private:
  DJetCorrResponse(const DJetCorrResponse &source);
  DJetCorrResponse& operator=(const DJetCorrResponse& source); 

  ClassDef(DJetCorrResponse, 1);
};
