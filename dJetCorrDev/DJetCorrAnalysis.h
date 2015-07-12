// This class projects the THnSparse results of the AliAnalysisTaskDmesonJetCorrelations into TH1/TH2/TH3
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
class TObject;
class TList;
class TDirectoryFile;
class TH1;

#include <TMap.h>

#include "DJetCorrBase.h"

class DJetCorrAnalysis : public DJetCorrBase {
  
 public:
  DJetCorrAnalysis();
  DJetCorrAnalysis(const char* train, const char* path = "$JETRESULTS");

  void   SetQAListName(const char* lname)      { fQAListName       = lname  ; }
  
  Bool_t Init();

  Bool_t GenerateQAHistograms();
  Bool_t GenerateDJetCorrHistograms(DJetCorrAnalysisParams* params);
  Bool_t GenerateDJetCorrHistograms(const char* paramName);
  Bool_t GenerateDJetCorrHistograms(Int_t i);
  Bool_t GenerateDJetCorrHistograms();
   
  Bool_t PlotTrackHistograms();
  Bool_t PlotDJetCorrHistograms(DJetCorrAnalysisParams* params);
  Bool_t PlotDJetCorrHistograms(const char* paramName);
  Bool_t PlotDJetCorrHistograms(Int_t i);
  Bool_t PlotDJetCorrHistograms();
    
 protected:
  void            ClearInputData();
  Int_t           GetAxisIndex(TString title) { return DJetCorrBase::GetAxisIndex(title, fDmesons); }
  
  Bool_t          ProjectQA();
  Bool_t          ProjectCorrD(DJetCorrAnalysisParams* params);
  Bool_t          ProjectDJetCorr(TString prefix, TString suffix, DJetCorrAnalysisParams* params, EMatchingStatus st, Int_t dptBin=-1, Int_t jetptBin=-1, Int_t dzBin=-1, Int_t minJetConst=0);

  Bool_t          PlotInvMassHistogramsVsDPt(DJetCorrAnalysisParams* params, EMatchingStatus st, Int_t jetptBin=-1, Int_t dzBin=-1);
  Bool_t          PlotInvMassHistogramsVsDz(DJetCorrAnalysisParams* params, Int_t dptBin=-1, Int_t jetptBin=-1);
  Bool_t          PlotInvMassHistogramArray(Int_t n, TH1** histos, const char* name, const char* xTitle,
                                            Double_t minMass, Double_t maxMass, Double_t pdgMass, Double_t massLimits,
                                            DJetCorrAnalysisParams* params, Bool_t doFit=kFALSE);
    
  Bool_t          LoadQAList();
  Bool_t          LoadTHnSparse();
  
  TString         fQAListName                ;//  QA list name

  // Axis titles
  TString         fDPtAxisTitle              ;//  d meson pt axis title
  TString         fDEtaAxisTitle             ;//  d meson eta axis title
  TString         fDPhiAxisTitle             ;//  d meson phi axis title
  TString         fDInvMassAxisTitle         ;//  d meson inv mass axis title
  TString         fD2prongInvMassAxisTitle   ;//  d meson 2-prong inv mass axis title
  TString         fDDeltaInvMassAxisTitle    ;//  d meson delta inv mass axis title
  TString         fDSoftPionPtAxisTitle      ;//  d meson soft pion pt axis title
  TString         fDzAxisTitle               ;//  d meson z axis title
  TString         fDeltaRAxisTitle           ;//  delta R axis title
  TString         fDeltaEtaAxisTitle         ;//  delta eta axis title
  TString         fDeltaPhiAxisTitle         ;//  delta phi axis title
  TString         fJetPtAxisTitle            ;//  jet pt axis title
  TString         fJetEtaAxisTitle           ;//  jet eta axis title
  TString         fJetPhiAxisTitle           ;//  jet phi axis title
  TString         fJetLeadPtAxisTitle        ;//  jet leading pt axis title
  TString         fJetAreaAxisTitle          ;//  jet area axis title
  TString         fJetConstAxisTitle         ;//  jet constituents axis title
  TString         fDeltaRDaughterAxisTitle   ;//  distance of a D meson daughter from the jet axis title
  TString         fMatchingStatusAxisTitle   ;//  distance of a D meson daughter from the jet axis title

  THnSparse      *fDmesons                   ;//! THnSparse contains the results
  TMap            fTHnAxisMap                ;//! mapping of axis titles with indexes
  TList          *fInputQAList               ;//! list contains the QA histograms
  TObjArray      *fMassFitters               ;//! array containing the mass fitter objects

 private:
   
  DJetCorrAnalysis(const DJetCorrAnalysis &source);
  DJetCorrAnalysis& operator=(const DJetCorrAnalysis& source); 

  ClassDef(DJetCorrAnalysis, 1);
};
