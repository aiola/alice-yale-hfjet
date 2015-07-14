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
  DJetCorrAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius, const char* tracksName, DJetCorrAnalysisType anaType);
  DJetCorrAnalysisParams(const DJetCorrAnalysisParams& p);
  
  void        SetInvMassRange(Double_t min, Double_t max)          { fInvMinMass = min        ; fInvMaxMass = max        ; }
  void        SetInvMassRange(Int_t pdg, Double_t nsigma);
  void        Set2ProngMassRange(Double_t min, Double_t max)       { f2ProngMinMass = min     ; f2ProngMaxMass = max     ; }
  void        Set2ProngMassRange(Int_t pdg, Double_t nsigma);
  void        SetDeltaInvMassRange(Double_t min, Double_t max)     { fDeltaInvMinMass = min   ; fDeltaInvMaxMass = max   ; }
  void        SetDeltaInvMassRange(Int_t pdg1, Int_t pdg2, Double_t nsigma);

  void        SetJetType(const char* type) { fJetType = type; }
  
  const char* GetName()              const { return fName.Data()                                    ; }
  const char* GetTitle()             const { return GetName()                                       ; }
  const char* GetJetType()           const { return fJetType.Data()                                 ; }
  const char* GetJetRadius()         const { return fJetRadius.Data()                               ; }
  const char* GetDmesonName()        const { return fDmesonName.Data()                              ; }
  const char* GetInputListName()     const { return fInputListName.Data()                           ; }
  Int_t       GetNDPtBins()          const { return fNDPtBins                                       ; }
  Double_t    GetDPtBin(Int_t i)     const { return i >= 0 && i <= fNDPtBins ? fDPtBins[i] : -1     ; }
  Double_t    GetMinDPt()            const { return fDPtBins[0]                                     ; }
  Double_t    GetMaxDPt()            const { return fDPtBins[fNDPtBins]                             ; }
  Int_t       GetNJetPtBins()        const { return fNJetPtBins                                     ; }
  Double_t    GetJetPtBin(Int_t i)   const { return i >= 0 && i <= fNJetPtBins ? fJetPtBins[i] : -1 ; }
  Double_t    GetMinJetPt()          const { return fJetPtBins[0]                                   ; }
  Double_t    GetMaxJetPt()          const { return fJetPtBins[fNJetPtBins]                         ; }
  Int_t       GetNzBins()            const { return fNzBins                                         ; }
  Double_t    GetzBin(Int_t i)       const { return i >= 0 && i <= fNzBins ? fzBins[i] : -1         ; }
  Double_t    GetMinZ()              const { return fzBins[0]                                       ; }
  Double_t    GetMaxZ()              const { return fzBins[fNzBins]                                 ; }
  Bool_t      IsD0()                 const { return (fDmesonName == "D0")                           ; }
  Bool_t      IsDStar()              const { return (fDmesonName == "DStar")                        ; }
  
  Bool_t      IsInInvMassRange(Double_t mass)    const { return (mass < fInvMaxMass && mass > fInvMinMass); }
  Double_t    GetInvMinMass()          const { return fInvMinMass                                   ; }
  Double_t    GetInvMaxMass()          const { return fInvMaxMass                                   ; }

  Bool_t      IsIn2ProngMassRange(Double_t mass)    const { return (mass < f2ProngMaxMass && mass > f2ProngMinMass); }
  Double_t    Get2ProngMinMass()       const { return f2ProngMinMass                                ; }
  Double_t    Get2ProngMaxMass()       const { return f2ProngMaxMass                                ; }
      
  Bool_t      IsInDeltaInvMassRange(Double_t dmass)    const { return (dmass < fDeltaInvMaxMass && dmass > fDeltaInvMinMass); }
  Double_t    GetDeltaInvMinMass()     const { return fDeltaInvMinMass                              ; }
  Double_t    GetDeltaInvMaxMass()     const { return fDeltaInvMaxMass                              ; }

  Double_t    GetMinDEta()             const { return fMinDEta                                      ; }
  Double_t    GetMaxDEta()             const { return fMaxDEta                                      ; }

  Int_t       GetMinJetConstituents()  const { return fMinJetConstituents                           ; }

  TString     GetCutString(Int_t st, Int_t dptBin=-1, Int_t jetptBin=-1, Int_t dzBin=-1);

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

  Double_t        fInvMinMass                ;//  inv min mass
  Double_t        fInvMaxMass                ;//  inv max mass
  Double_t        f2ProngMinMass             ;//  2-prong min mass (D* -> D0pi)
  Double_t        f2ProngMaxMass             ;//  2-prong max mass (D* -> D0pi)  
  Double_t        fDeltaInvMinMass           ;//  delta inv min mass (D* -> D0pi)
  Double_t        fDeltaInvMaxMass           ;//  delta inv max mass (D* -> D0pi)

  Double_t        fMinDEta                   ;//  min eta of D meson
  Double_t        fMaxDEta                   ;//  max eta of D meson
  Int_t           fMinJetConstituents        ;//  min number of jet constituents
  MassFitter::EMassFitTypeSig fMassFitTypeSig;//
  MassFitter::EMassFitTypeBkg fMassFitTypeBkg;//
  
 private:
 
  DJetCorrAnalysisParams& operator=(const DJetCorrAnalysisParams& source); 

  ClassDef(DJetCorrAnalysisParams, 1);
};
