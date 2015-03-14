// This class projects the THnSparse results of the AliAnalysisTaskDmesonJetCorrelations into TH1/TH2/TH3
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

class TString;
class TFile;
class TList;
class TMap;
class THnSparse;

class DJetCorrAnalysis {
  
 public:
  DJetCorrAnalysis();
  DJetCorrAnalysis(const char* train, const char* path = "$JETRESULTS");

  void   SetInputTrain(const char* train)      { fTrainName        = train  ; CloseInputFile(); }
  void   SetInputPath(const char* path)        { fInputPath        = path   ; CloseInputFile(); }
  void   SetInputFileName(const char* fname)   { fInputFileName    = fname  ; CloseInputFile(); }
  void   SetInputListName(const char* lname)   { fInputListName    = lname  ; }
  void   SetOutputFileName(const char* fname)  { fOutputFileName   = fname  ; }
  void   SetOutputPath(const char* path)       { fOutputPath       = path   ; }
  void   SetOverwrite(Bool_t ow = kTRUE)       { fOverwrite        = ow     ; }
  void   SetUseTest(Bool_t t = kTRUE)          { fUseTestResults   = t      ; }

  void   SetAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius);
  
  Bool_t Init();
  Bool_t GenerateDJetCorrHistograms(const char* dmeson, const char* jetType, const char* jetRadius);
  Bool_t GenerateQAHistograms();
  Bool_t SaveOutputFile();
  
 protected:

  void            ClearInputData();
  Int_t           GetAxisIndex(TString title);

  Bool_t          ProjectQA();
  Bool_t          ProjectCorrD();
  Bool_t          ProjectDJetCorr(TString prefix, TString suffix, Bool_t /*doCorrPlots*/, 
                                  Double_t minJetPt, Double_t maxJetPt,
                                  Double_t minDPt, Double_t maxDPt, Double_t minDEta, Double_t maxDEta);

  void            GenerateAxisMap(THnSparse* hn);
  Bool_t          OpenInputFile();
  Bool_t          LoadInputList();
  Bool_t          LoadQAList();
  Bool_t          LoadTHnSparse();

  void            CloseInputFile();

  Bool_t          fUseTestResults            ;//  whether to use test results
  TString         fTrainName                 ;//  train name
  TString         fInputPath                 ;//  train path
  TString         fInputFileName             ;//  input file name
  TString         fInputDirFileName          ;//  input dir file name
  TString         fInputListName             ;//  input list name
  TString         fOutputPath                ;//  output path
  TString         fOutputFileName            ;//  output file name
  Bool_t          fOverwrite                 ;//  whether the output file should be overwritten
  TString         fJetType                   ;//  jet type "Charged" or "Full"
  TString         fJetRadius                 ;//  jet radius R020, R030, etc.
  TString         fDmesonName                ;//  "D0" or "DStar" etc.
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

  Int_t           fNDPtBins                  ;//  number of D pt bins
  Double_t       *fDPtBins                   ;//[fNDPtBins] D pt bins
  Int_t           fNJetPtBins                ;//  number of jet pt bins
  Double_t       *fJetPtBins                 ;//[fNJetPtBins] jet pt bins
  
  TMap            fTHnAxisMap                ;//! mapping of axis titles with indexes
  Bool_t          fTHnSparseMapGenerated     ;//! whether or not the axis map has been generated for the THnSparse
  TFile          *fInputFile                 ;//! input file
  TDirectoryFile *fInputDirectoryFile        ;//! input directory file
  TList          *fInputList                 ;//! list contains the input hisograms
  TList          *fInputQAList               ;//! list contains the QA histograms
  THnSparse      *fDmesons                   ;//! THnSparse contains the results
  TList          *fOutputList                ;//! list contains the output histograms
  
  static const Double_t fgkEpsilon;
};
