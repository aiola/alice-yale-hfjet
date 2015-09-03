// This class projects the THnSparse results of the AliAnalysisTaskDmesonJetCorrelations into TH1/TH2/TH3
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

class TString;

#include <TObject.h>
#include <MassFitter.h>

class DJetCorrAnalysisParams : public TObject
{
 public:
  enum {kAnyMatchingStatus, kNotMatched, kMatched};
  enum DJetCorrAnalysisType {KUndefinedAna=0, kInvMassAna, kResponseMatrixAna};

  DJetCorrAnalysisParams();
  DJetCorrAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius, const char* tracksName, DJetCorrAnalysisType anaType, Bool_t isMC=kFALSE, Bool_t isBkgSub=kFALSE);
  DJetCorrAnalysisParams(const DJetCorrAnalysisParams& p);
  
  void        SetInvMassRange(Double_t min, Double_t max)          { fInvMinMass = min        ; fInvMaxMass = max        ; }
  void        SetInvMassRange(Int_t pdg, Double_t nsigma);
  void        Set2ProngMassRange(Double_t min, Double_t max, Int_t i)       { f2ProngMinMass[i] = min     ; f2ProngMaxMass[i] = max     ; }
  void        Set2ProngMassRange(Int_t pdg, Double_t nsigma, Int_t i);
  void        SetDeltaInvMassRange(Double_t min, Double_t max)     { fDeltaInvMinMass = min   ; fDeltaInvMaxMass = max   ; }
  void        SetDeltaInvMassRange(Int_t pdg1, Int_t pdg2, Double_t nsigma);

  void        SetIsMC(Bool_t mc)      { fIsMC     = mc   ; }
  void        SetIsBkgSub(Bool_t bkg) { fIsBkgSub = bkg  ; }

  Bool_t      IsBkgSub()        const { return fIsBkgSub ; }
  Bool_t      IsMC()            const { return fIsMC     ; }

  void        SetJetType(const char* type) { fJetType = type; }
  
  const char* GetName()              const { return fName.Data()                                    ; }
  const char* GetTitle()             const { return GetName()                                       ; }
  const char* GetJetType()           const { return fJetType.Data()                                 ; }
  const char* GetJetRadius()         const { return fJetRadius.Data()                               ; }
  const char* GetDmesonName()        const { return fDmesonName.Data()                              ; }
  const char* GetInputListName()     const { return fInputListName.Data()                           ; }
  const char* GetTruthInputListName()const { return fTruthInputListName.Data()                      ; }
  Int_t       GetNDPtBins()          const { return fNDPtBins                                       ; }
  const Double_t* GetDPtBins()       const { return fDPtBins                                        ; }
  Double_t    GetDPtBin(Int_t i)     const { return i >= 0 && i <= fNDPtBins ? fDPtBins[i] : -1     ; }
  Double_t    GetMinDPt()            const { return fDPtBins[0]                                     ; }
  Double_t    GetMaxDPt()            const { return fDPtBins[fNDPtBins]                             ; }
  Int_t       GetNJetPtBins()        const { return fNJetPtBins                                     ; }
  const Double_t* GetJetPtBins()     const { return fJetPtBins                                      ; }
  Double_t    GetJetPtBin(Int_t i)   const { return i >= 0 && i <= fNJetPtBins ? fJetPtBins[i] : -1 ; }
  Double_t    GetMinJetPt()          const { return fJetPtBins[0]                                   ; }
  Double_t    GetMaxJetPt()          const { return fJetPtBins[fNJetPtBins]                         ; }
  Int_t       GetNzBins()            const { return fNzBins                                         ; }
  const Double_t* GetzBins()         const { return fzBins                                          ; }
  Double_t    GetzBin(Int_t i)       const { return i >= 0 && i <= fNzBins ? fzBins[i] : -1         ; }
  Double_t    GetMinZ()              const { return fzBins[0]                                       ; }
  Double_t    GetMaxZ()              const { return fzBins[fNzBins]                                 ; }
  Bool_t      IsD0()                 const { return (fDmesonName == "D0")                           ; }
  Bool_t      IsDStar()              const { return (fDmesonName == "DStar")                        ; }
  
  Bool_t      IsInInvMassRange(Double_t mass)    const { return (mass < fInvMaxMass && mass > fInvMinMass); }
  Double_t    GetInvMinMass()          const { return fInvMinMass                                   ; }
  Double_t    GetInvMaxMass()          const { return fInvMaxMass                                   ; }

  Bool_t      IsIn2ProngMassRange(Double_t mass, Int_t zBin=-1, Int_t dptBin=-1) const { return (mass < Get2ProngMaxMass(zBin, dptBin) && mass > Get2ProngMinMass(zBin, dptBin)); }
  Double_t    Get2ProngMinMass(Int_t zBin=-1, Int_t dptBin=-1) const;
  Double_t    Get2ProngMaxMass(Int_t zBin=-1, Int_t dptBin=-1) const;
      
  Bool_t      IsInDeltaInvMassRange(Double_t dmass)    const { return (dmass < fDeltaInvMaxMass && dmass > fDeltaInvMinMass); }
  Double_t    GetDeltaInvMinMass()     const { return fDeltaInvMinMass                              ; }
  Double_t    GetDeltaInvMaxMass()     const { return fDeltaInvMaxMass                              ; }

  Double_t    GetMinDEta()             const { return fMinDEta                                      ; }
  Double_t    GetMaxDEta()             const { return fMaxDEta                                      ; }

  Int_t       GetMinJetConstituents()  const { return fMinJetConstituents                           ; }

  TString     GetCutString(Int_t st, Int_t dptBin=-1, Int_t jetptBin=-1, Int_t dzBin=-1);

  void GetDPtBinRange(Double_t& minDPt, Double_t& maxDPt, Int_t dptBin=-1) const;
  void GetzBinRange(Double_t& minZ, Double_t& maxZ, Int_t zBin=-1) const;
  void GetJetPtBinRange(Double_t& minJetPt, Double_t& maxJetPt, Int_t jetPtBin=-1) const;

  MassFitter* CreateMassFitter(const char* name) const;
  
 protected:
  TString         fName                      ;//  object name
  Int_t           fNDPtBins                  ;//  number of D pt bins
  Double_t       *fDPtBins                   ;//[fNDPtBins+1] D pt bins
  Int_t           fNJetPtBins                ;//  number of jet pt bins
  Double_t       *fJetPtBins                 ;//[fNJetPtBins+1] jet pt bins
  Int_t           fNzBins                    ;//  number of z bins
  Double_t       *fzBins                     ;//[fNzBins+1] z bins
  
  TString         fJetType                   ;//  jet type "Charged" or "Full"
  TString         fJetRadius                 ;//  jet radius R020, R030, etc.
  TString         fDmesonName                ;//  "D0" or "DStar" etc.
  TString         fTracksName                ;//  usually "tracks" or "AODFilterTracks"
  TString         fInputListName             ;//  input list name
  TString         fTruthInputListName        ;//  truth input list name

  Double_t        fInvMinMass                ;//  inv min mass
  Double_t        fInvMaxMass                ;//  inv max mass
  Double_t       *f2ProngMinMass             ;//[fNzBins] 2-prong min mass (D* -> D0pi)
  Double_t       *f2ProngMaxMass             ;//[fNzBins] 2-prong max mass (D* -> D0pi)
  Double_t        fDeltaInvMinMass           ;//  delta inv min mass (D* -> D0pi)
  Double_t        fDeltaInvMaxMass           ;//  delta inv max mass (D* -> D0pi)

  Double_t        fMinDEta                   ;//  min eta of D meson
  Double_t        fMaxDEta                   ;//  max eta of D meson
  Int_t           fMinJetConstituents        ;//  min number of jet constituents
  MassFitter::EMassFitTypeSig fMassFitTypeSig;//
  MassFitter::EMassFitTypeBkg fMassFitTypeBkg;//

  Bool_t          fIsMC                      ;//
  Bool_t          fIsBkgSub                  ;//
  
 private:
 
  DJetCorrAnalysisParams& operator=(const DJetCorrAnalysisParams& source); 

  ClassDef(DJetCorrAnalysisParams, 1);
};
