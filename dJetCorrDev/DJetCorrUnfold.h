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
  DJetCorrUnfold();
  DJetCorrUnfold(DJetCorrAnalysis* ana, DJetCorrResponse* resp);

  void   SetDataParamIndex(Int_t i)    { fDataParamIndex    = i; }
  void   SetRespParamIndex(Int_t i)    { fRespParamIndex    = i; }
  void   SetUseEfficiency(Bool_t b)    { fUseEfficiency     = b; }
  void   SetUseKinEfficiency(Bool_t b) { fUseKinEfficiency  = b; }

  Bool_t Start();
  Bool_t Unfold();
  Bool_t PrepareData();
  Bool_t PrepareResponse();
  RooUnfoldResponse* GenerateRooUnfoldResponse();
  void AddEfficiency(RooUnfoldResponse* resp, TH2* misses);
  
 protected:
  Int_t                    fDataParamIndex   ; //
  Int_t                    fRespParamIndex   ; //
  Bool_t                   fForceRegeneration; //
  Bool_t                   fUseEfficiency    ; //
  Bool_t                   fUseKinEfficiency ; //

  TH2*                     fTruth            ; //
  TH2*                     fMeasured         ; //
  TH2*                     fUnfolded         ; //
  TH2*                     fRefolded         ; //
  
  TH2*                     fResponseTruth    ; //
  TH2*                     fResponseMeasured ; //
  THnSparse*               fResponseMatrix   ; //
  TH2*                     fResponseMisses   ; //
  TH2*                     fResponseKinMisses; //

  DJetCorrAnalysis*        fAnalysis         ; //! analysis results
  DJetCorrResponse*        fResponse         ; //! response
    
 private: 
  DJetCorrUnfold(const DJetCorrUnfold &source);
  DJetCorrUnfold& operator=(const DJetCorrUnfold& source); 

  ClassDef(DJetCorrUnfold, 1);
};
