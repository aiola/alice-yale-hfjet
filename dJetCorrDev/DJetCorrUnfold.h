// Class to unfold results from DJetCorrAnalysis using response matrix from DJetCorrResponse
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

class DJetCorrAnalysis;
class DJetCorrResponse;
class RooUnfoldResponse;
class TH2;
class THnSparse;

#include <TMath.h>

#include "DJetCorrBase.h"

class DJetCorrUnfold : public DJetCorrBase {

public:
  enum EAxisType_t {
    kJetPtAxis   = 0,  // x axis
    kDmesonZAxis = 1   // y axis
  };

  class UnfoldingParams : public TNamed {

  public:

    enum EUnfoldingMethod_t {
      kNoUnfolding = 0,
      kBayesian    = 1,
      kSVD         = 2,
      kChi2        = 3
    };

    UnfoldingParams();
    UnfoldingParams(EUnfoldingMethod_t m, Int_t min=0, Int_t max=0, Int_t step=0);
    UnfoldingParams(const UnfoldingParams &source);
    UnfoldingParams& operator=(const UnfoldingParams& source);

    void   Set(EUnfoldingMethod_t m, Int_t min=0, Int_t max=0, Int_t step=0);

    void   SetRegParam(Int_t min, Int_t max, Int_t step);
    void   SetUnfoldingMethod(EUnfoldingMethod_t m)      { fUnfoldingMethod = m; }

    Int_t  RegParam(Int_t i)                       const { return i < fNRegParams ? fMinRegParam + i * fRegParamStep : -1; }
    Int_t  GetMinRegParam()                        const { return fMinRegParam; }
    Int_t  GetMaxRegParam()                        const { return fMinRegParam + (fNRegParams-1)*fRegParamStep; }
    Int_t  GetNRegParams()                         const { return fNRegParams ; }
    Int_t  GetRegParamStep()                       const { return fRegParamStep ; }
    EUnfoldingMethod_t GetUnfoldingMethod()        const { return fUnfoldingMethod ; }

    void   ResetCurrentRegIndex(Int_t i = 0)             { fCurrentRegIndex = i; }
    Int_t  NextRegParam();

    TString fgkUnfoldingMethodNames[4];

  protected:
    EUnfoldingMethod_t  fUnfoldingMethod         ; //
    Int_t               fMinRegParam             ; //
    Int_t               fRegParamStep            ; //
    Int_t               fNRegParams              ; //

    Int_t               fCurrentRegIndex         ; //!

  private:
    ClassDef(UnfoldingParams, 1);
  };

  DJetCorrUnfold();
  DJetCorrUnfold(DJetCorrAnalysis* ana, DJetCorrResponse* resp);

  void   SetDataParamIndex(Int_t i)    { fDataParamIndex    = i; }
  void   SetRespParamIndex(Int_t i)    { fRespParamIndex    = i; }
  void   SetUseEfficiency(Bool_t b)    { fUseEfficiency     = b; }
  void   SetUseKinEfficiency(Bool_t b) { fUseKinEfficiency  = b; }

  void   SetRegParam(Int_t min, Int_t max, Int_t step)                                 { fUnfoldingParams.SetRegParam(min, max, step)         ; }
  void   SetRegParam1D(EAxisType_t axis, Int_t min, Int_t max, Int_t step)             { fUnfoldingParams1D[axis].SetRegParam(min, max, step) ; }
  void   SetUnfoldingMethod1D(EAxisType_t axis, UnfoldingParams::EUnfoldingMethod_t m) { fUnfoldingParams1D[axis].SetUnfoldingMethod(m)       ; }

  Bool_t GenerateRooUnfoldResponse();
  Bool_t Unfold();

  Bool_t GenerateRooUnfoldResponse1D(EAxisType_t axis);
  Bool_t Unfold1D(EAxisType_t axis);

  void Unfold1DBayes(TH1* measured, RooUnfoldResponse* resp, UnfoldingParams params);

  Bool_t Fold(TH1* hist, RooUnfoldResponse* resp, TString name, TString title);

  Bool_t Start();
  Bool_t MakeProjectionsBeforeUnfolding();
  Bool_t MakeProjectionsAfterUnfolding();
  Bool_t MakeProjections(TString hlabel);
  Bool_t MakeProjectionsTHnSparse(TString hname);

  Bool_t PrepareData();
  Bool_t PrepareResponse();
  Bool_t SaveOutputFile();
  Bool_t MakePlots();
  void PlotResponse();

  void PlotUnfolded();
  void PlotRefolded();
  void PlotFolded();

  void PlotUnfolded(EAxisType_t axis);
  void PlotRefolded(EAxisType_t axis);
  void PlotFolded(EAxisType_t axis);

  void PlotMeasured2D();
  void PlotFolded2D();
  void PlotTruth2D();

  void AddEfficiency(RooUnfoldResponse* resp, TH2* misses);
  void AddEfficiency1D(RooUnfoldResponse* resp, TH1* misses);
  
  TH2*       GetResponseTruth()              { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseTruthName())); }
  TH2*       GetResponseMeasured()           { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseMeasuredName())); }
  THnSparse* GetResponseMatrix()             { return fOutputList == 0 ? 0 : static_cast<THnSparse*>(fOutputList->FindObject(GetResponseMatrixName())); }
  TH2*       GetResponseMisses()             { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseMissesName())); }
  TH2*       GetResponseKinMisses()          { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseKinMissesName())); }

  TH2*       GetUnfolded(Int_t regParam)     { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetUnfoldedName(regParam))); }
  TH2*       GetRefolded(Int_t regParam)     { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetRefoldedName(regParam))); }
  TH2*       GetFolded()                     { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetFoldedName())); }

  TH1*       GetTruthProj(EAxisType_t axis, Int_t bin)                    { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetTruthProjName(axis, bin))); }
  TH1*       GetMeasuredProj(EAxisType_t axis, Int_t bin)                 { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetMeasuredProjName(axis, bin))); }
  TH1*       GetUnfoldedProj(Int_t regParam, EAxisType_t axis, Int_t bin) { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetUnfoldedProjName(regParam, axis, bin))); }
  TH1*       GetRefoldedProj(Int_t regParam, EAxisType_t axis, Int_t bin) { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetRefoldedProjName(regParam, axis, bin))); }
  TH1*       GetFoldedProj(EAxisType_t axis, Int_t bin)                   { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetFoldedProjName(axis, bin))); }

  TH1*       GetResponseTruthProj(EAxisType_t axis, Int_t bin)            { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetResponseTruthProjName(axis, bin))); }
  TH1*       GetResponseMeasuredProj(EAxisType_t axis, Int_t bin)         { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetResponseMeasuredProjName(axis, bin))); }
  TH2*       GetResponseMatrixProj(EAxisType_t axis, Int_t bin)           { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseMatrixProjName(axis, bin))); }
  TH1*       GetResponseMissesProj(EAxisType_t axis, Int_t bin)           { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseMissesProjName(axis, bin))); }
  TH1*       GetResponseKinMissesProj(EAxisType_t axis, Int_t bin)        { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseKinMissesProjName(axis, bin))); }

  TString    GetResponseTruthName()              ;
  TString    GetResponseMeasuredName()           ;
  TString    GetResponseMatrixName()             ;
  TString    GetResponseMissesName()             ;
  TString    GetResponseKinMissesName()          ;

  TString    GetTruthName(Int_t=0)               ;
  TString    GetMeasuredName(Int_t=0)            ;
  TString    GetUnfoldedName(Int_t regParam)     ;
  TString    GetRefoldedName(Int_t regParam)     ;
  TString    GetFoldedName()                     ;

  TString    GetTruthProjName(EAxisType_t axis, Int_t bin)                    ;
  TString    GetMeasuredProjName(EAxisType_t axis, Int_t bin)                 ;
  TString    GetUnfoldedProjName(Int_t regParam, EAxisType_t axis, Int_t bin) ;
  TString    GetRefoldedProjName(Int_t regParam, EAxisType_t axis, Int_t bin) ;
  TString    GetFoldedProjName(EAxisType_t axis, Int_t bin)                   ;

  TString    GetResponseTruthProjName(EAxisType_t axis, Int_t bin)            ;
  TString    GetResponseMeasuredProjName(EAxisType_t axis, Int_t bin)         ;
  TString    GetResponseMatrixProjName(EAxisType_t axis, Int_t bin)           ;
  TString    GetResponseMissesProjName(EAxisType_t axis, Int_t bin)           ;
  TString    GetResponseKinMissesProjName(EAxisType_t axis, Int_t bin)        ;

  Int_t      GetNbinsX();
  Int_t      GetNbinsY();
  Int_t      GetNbins(Int_t axis);

  void       SavePlot(TCanvas* canvas);

 protected:
  Int_t                    fDataParamIndex          ; //
  Int_t                    fRespParamIndex          ; //
  Bool_t                   fForceRegeneration       ; //
  Bool_t                   fUseEfficiency           ; //
  Bool_t                   fUseKinEfficiency        ; //
  UnfoldingParams          fUnfoldingParams         ; //
  UnfoldingParams          fUnfoldingParams1D[2]    ; //

  static const TString     fgkAxisLabels[2]         ; //! axis labels
  DJetCorrAnalysis*        fAnalysis                ; //! analysis results
  DJetCorrResponse*        fResponse                ; //! response
  RooUnfoldResponse*       fRooUnfoldResponse       ; //! RooUnfold response matrix
  RooUnfoldResponse*       fRooUnfoldResponse1D[2]  ; //! RooUnfold response matrix
    
 private: 
  DJetCorrUnfold(const DJetCorrUnfold &source);
  DJetCorrUnfold& operator=(const DJetCorrUnfold& source); 

  ClassDef(DJetCorrUnfold, 1);
};
