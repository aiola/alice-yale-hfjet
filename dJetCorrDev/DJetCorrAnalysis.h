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

#include "DJetCorrBase.h"

class DJetCorrAnalysis : public DJetCorrBase {
  
 public:
  enum EInvMassPlotNorm { kPureCounts, kDivideByBinWidth, kNormalizeBackground };
  
  DJetCorrAnalysis();
  DJetCorrAnalysis(const char* train, const char* path = "$JETRESULTS");

  void   SetQAListName(const char* lname)        { fQAListName      = lname  ; }
  void   SetInvMassPlotNorm(EInvMassPlotNorm n)  { fInvMassPlotNorm = n      ; }
  
  Bool_t Init();
  Bool_t SaveOutputFile();
  
  Bool_t GenerateQAHistograms();
  Bool_t GenerateDJetCorrHistograms(DJetCorrAnalysisParams* params);
  Bool_t GenerateDJetCorrHistograms(const char* paramName);
  Bool_t GenerateDJetCorrHistograms(Int_t i);
  Bool_t GenerateDJetCorrHistograms();
   
  Bool_t PlotTrackHistograms();
  Bool_t PlotDJetCorrHistograms(DJetCorrAnalysisParams* params);
  Bool_t PlotDJetCorrHistograms(const char* paramName);
  Bool_t PlotDJetCorrHistograms(Int_t i);
  Bool_t PlotDJetCorrHistograms(Bool_t forceRefit=kFALSE);

  Bool_t PlotDPtSpectraVsJetPt(DJetCorrAnalysisParams* params, Bool_t eventScaling=kFALSE);
  Bool_t PlotDPtSpectraVsDz(DJetCorrAnalysisParams* params, Bool_t eventScaling=kFALSE);
  Bool_t PlotDPtSpectraVsMatchingStatus(DJetCorrAnalysisParams* params, Bool_t eventScaling=kFALSE);
  Bool_t PlotDzSpectraVsJetPt(DJetCorrAnalysisParams* params, Bool_t eventScaling=kFALSE);

  Bool_t Regenerate();

  Bool_t LoadTruthList(DJetCorrAnalysisParams* params);
  Bool_t ProjectTruthSpectrum();
  Bool_t ProjectTruthSpectrum(DJetCorrAnalysisParams* params);
  Bool_t ProjectTruthSpectrum(TString prefix, TString suffix, DJetCorrAnalysisParams* params, EMatchingStatus st, Int_t dptBin=-1, Int_t jetptBin=-1, Int_t dzBin=-1, Int_t minJetConst=0);

  void GenerateMeaduredSpectrum(DJetCorrAnalysisParams* params);

  TString GetDzTruthName(Int_t p);
  TString GetDzMeasuredName(Int_t p);

  TString GetDPtTruthName(Int_t p, const char* matching="AnyMatchingStatus");
  TString GetDPtMeasuredName(Int_t p, const char* matching="AnyMatchingStatus");

  TString GetDEtaTruthName(Int_t p, const char* matching="AnyMatchingStatus");
  TString GetDEtaMeasuredName(Int_t p, const char* matching="AnyMatchingStatus");

  TString GetJetPtTruthName(Int_t p);
  TString GetJetPtMeasuredName(Int_t p);

 protected:
  Bool_t          LoadOutputHistograms();
  Bool_t          ClearInputData();
  Int_t           GetAxisIndex(TString title) { return DJetCorrBase::GetAxisIndex(title, fDmesons); }
  
  Bool_t          ProjectQA();
  Bool_t          ProjectCorrD(DJetCorrAnalysisParams* params);
  Bool_t          ProjectDJetCorr(TString prefix, TString suffix, DJetCorrAnalysisParams* params, EMatchingStatus st, Int_t dptBin=-1, Int_t jetptBin=-1, Int_t dzBin=-1, Int_t minJetConst=0);

  Bool_t          PlotInvMassHistogramsVsDPt(DJetCorrAnalysisParams* params, EMatchingStatus st, Int_t jetptBin=-1, Int_t dzBin=-1);
  Bool_t          PlotInvMassHistogramsVsDz(DJetCorrAnalysisParams* params, Int_t dptBin=-1, Int_t jetptBin=-1);
  Bool_t          PlotInvMassHistogramArray(Int_t n, TH1** histos, const char* name, const char* xTitle,
                                            Double_t minMass, Double_t maxMass, Double_t pdgMass, Double_t* minMassSel, Double_t* maxMassSel,
                                            DJetCorrAnalysisParams* params, Bool_t doFit=kFALSE, TObjArray* extraInfo=0, TH1* histSpectrum=0);
  
  Bool_t          PlotSpectra(Int_t n, TH1** histSpectra, const char* name, Bool_t logY);
    
  Bool_t          LoadQAList();
  Bool_t          LoadTHnSparse();
  
  TString           fQAListName              ;//  QA list name
  EInvMassPlotNorm  fInvMassPlotNorm         ;//  type of normalization used in invariant mass plots
  Bool_t            fForceRefit              ;//  force refit of invariant mass plots

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
  TString         fMatchingStatusAxisTitle   ;//  matching status

  THnSparse      *fDmesons                   ;//! THnSparse contains the results
  TList          *fInputQAList               ;//! list contains the QA histograms
  TList          *fMassFitters               ;//! list containing the mass fitter objects

 private:
   
  DJetCorrAnalysis(const DJetCorrAnalysis &source);
  DJetCorrAnalysis& operator=(const DJetCorrAnalysis& source); 

  ClassDef(DJetCorrAnalysis, 1);
};
