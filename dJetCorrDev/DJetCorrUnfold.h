// Class to unfold results from DJetCorrAnalysis using response matrix from DJetCorrResponse
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

class DJetCorrAnalysis;
class DJetCorrResponse;
class RooUnfoldResponse;
class TH2;
class THnSparse;

class DJetCorrUnfold : public TNamed {
  
 public:
  enum EAxisType_t {
    kJetPtAxis   = 0,  // x axis
    kDmesonZAxis = 1   // y axis
  };

  DJetCorrUnfold();
  DJetCorrUnfold(DJetCorrAnalysis* ana, DJetCorrResponse* resp);

  void   SetDataParamIndex(Int_t i)    { fDataParamIndex    = i; }
  void   SetRespParamIndex(Int_t i)    { fRespParamIndex    = i; }
  void   SetUseEfficiency(Bool_t b)    { fUseEfficiency     = b; }
  void   SetUseKinEfficiency(Bool_t b) { fUseKinEfficiency  = b; }

  Bool_t Start();
  Bool_t Unfold();
  Bool_t MakeProjections();
  Bool_t MakeProjections(TString hlabel);
  Bool_t Fold(TH2* hist, TString name, TString title);
  Bool_t PrepareData();
  Bool_t PrepareResponse();
  Bool_t GenerateRooUnfoldResponse();
  void AddEfficiency(RooUnfoldResponse* resp, TH2* misses);
  
  TH2*       GetResponseTruth()              { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseTruthName())); }
  TH2*       GetResponseMeasured()           { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseMeasuredName())); }
  THnSparse* GetResponseMatrix()             { return fOutputList == 0 ? 0 : static_cast<THnSparse*>(fOutputList->FindObject(GetResponseMatrixName())); }
  TH2*       GetResponseMisses()             { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseMissesName())); }
  TH2*       GetResponseKinMisses()          { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetResponseKinMissesName())); }

  TH2*       GetTruth()                      { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetTruthName())); }
  TH2*       GetMeasured()                   { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetMeasuredName())); }
  TH2*       GetUnfolded(Int_t regParam)     { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetUnfoldedName(regParam))); }
  TH2*       GetRefolded(Int_t regParam)     { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetRefoldedName(regParam))); }
  TH2*       GetFolded()                     { return fOutputList == 0 ? 0 : static_cast<TH2*>(fOutputList->FindObject(GetFoldedName())); }

  TH1*       GetTruthProj(EAxisType_t axis, Int_t bin)                    { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetTruthProjName(axis, bin))); }
  TH1*       GetMeasuredProj(EAxisType_t axis, Int_t bin)                 { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetMeasuredProjName(axis, bin))); }
  TH1*       GetUnfoldedProj(Int_t regParam, EAxisType_t axis, Int_t bin) { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetUnfoldedProjName(regParam, axis, bin))); }
  TH1*       GetRefoldedProj(Int_t regParam, EAxisType_t axis, Int_t bin) { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetRefoldedProjName(regParam, axis, bin))); }
  TH1*       GetFoldedProj(EAxisType_t axis, Int_t bin)                   { return fOutputList == 0 ? 0 : static_cast<TH1*>(fOutputList->FindObject(GetFoldedProjName(axis, bin))); }

  TString    GetResponseTruthName()              ;
  TString    GetResponseMeasuredName()           ;
  TString    GetResponseMatrixName()             ;
  TString    GetResponseMissesName()             ;
  TString    GetResponseKinMissesName()          ;

  TString    GetTruthName()                      ;
  TString    GetMeasuredName()                   ;
  TString    GetUnfoldedName(Int_t regParam)     ;
  TString    GetRefoldedName(Int_t regParam)     ;
  TString    GetFoldedName()                     ;

  TString    GetTruthProjName(EAxisType_t axis, Int_t bin)                    ;
  TString    GetMeasuredProjName(EAxisType_t axis, Int_t bin)                 ;
  TString    GetUnfoldedProjName(Int_t regParam, EAxisType_t axis, Int_t bin) ;
  TString    GetRefoldedProjName(Int_t regParam, EAxisType_t axis, Int_t bin) ;
  TString    GetFoldedProjName(EAxisType_t axis, Int_t bin)                   ;


 protected:
  Int_t                    fDataParamIndex   ; //
  Int_t                    fRespParamIndex   ; //
  Bool_t                   fForceRegeneration; //
  Bool_t                   fUseEfficiency    ; //
  Bool_t                   fUseKinEfficiency ; //
  Int_t                    fMinRegParam      ; //
  Int_t                    fMaxRegParam      ; //
  Int_t                    fRegParamStep     ; //

  static const TString     fgkAxisLabels[2]  ; //! axis labels
  DJetCorrAnalysis*        fAnalysis         ; //! analysis results
  DJetCorrResponse*        fResponse         ; //! response
  RooUnfoldResponse*       fRooUnfoldResponse; //! RooUnfold response matrix
  THashList*               fOutputList       ; //! output list
    
 private: 
  DJetCorrUnfold(const DJetCorrUnfold &source);
  DJetCorrUnfold& operator=(const DJetCorrUnfold& source); 

  ClassDef(DJetCorrUnfold, 1);
};
