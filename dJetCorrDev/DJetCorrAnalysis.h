// This class projects the THnSparse results of the AliAnalysisTaskDmesonJetCorrelations into TH1/TH2/TH3
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

class TString;
class TFile;
class TList;
class TMap;
class THnSparse;
class TCanvas;
class TLegend;
class TVirtualPad;
class TPaveText;
class DJetCorrAnalysisParams;
class TObject;
class TList;

class DJetCorrAnalysis : public TObject {
  
 public:
  DJetCorrAnalysis();
  DJetCorrAnalysis(const char* train, const char* path = "$JETRESULTS");

  void   SetInputTrain(const char* train)      { fTrainName        = train  ; CloseInputFile(); }
  void   SetInputPath(const char* path)        { fInputPath        = path   ; CloseInputFile(); }
  void   SetInputFileName(const char* fname)   { fInputFileName    = fname  ; CloseInputFile(); }
  void   SetOutputFileName(const char* fname)  { fOutputFileName   = fname  ; }
  void   SetQAListName(const char* lname)      { fQAListName       = lname  ; }
  void   SetOutputPath(const char* path)       { fOutputPath       = path   ; }
  void   SetOverwrite(Bool_t ow = kTRUE)       { fOverwrite        = ow     ; }
  void   SetPlotFormat(const char* f)          { fPlotFormat       = f      ; }
  void   SetSavePlots(Bool_t s)                { fSavePlots        = s      ; }

  void   AddAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius, const char* tracksName);
  void   AddAnalysisParams(DJetCorrAnalysisParams* params);
  
  Bool_t Init();

  Bool_t GenerateQAHistograms();
  Bool_t GenerateDJetCorrHistograms(DJetCorrAnalysisParams* params);
  Bool_t GenerateDJetCorrHistograms(const char* paramName);
  Bool_t GenerateDJetCorrHistograms(Int_t i);
  Bool_t GenerateDJetCorrHistograms();
  
  Bool_t SaveOutputFile();
  
  Bool_t PlotTrackHistograms();
  Bool_t PlotDJetCorrHistograms(DJetCorrAnalysisParams* params);
  Bool_t PlotDJetCorrHistograms(const char* paramName);
  Bool_t PlotDJetCorrHistograms(Int_t i);
  Bool_t PlotDJetCorrHistograms();
  
 protected:

  void            ClearInputData();
  Int_t           GetAxisIndex(TString title);

  Bool_t          ProjectQA();
  Bool_t          ProjectCorrD(DJetCorrAnalysisParams* params);
  Bool_t          ProjectDJetCorr(TString prefix, TString suffix, Bool_t doCorrPlots, Bool_t doPtPlots, Bool_t dozPlots, 
                                  Double_t minJetPt, Double_t maxJetPt,
                                  Double_t minDPt, Double_t maxDPt, Double_t minz, Double_t maxz, Double_t minDEta, Double_t maxDEta,
                                  Double_t minInvMass, Double_t maxInvMass, Double_t min2ProngMass, Double_t max2ProngMass, Double_t minDeltaInvMass, Double_t maxDeltaInvMass);
  Bool_t          GenerateRatios(const char* nname, const char* dname);
  Bool_t          PlotInvMassHistogramsVsDPt(DJetCorrAnalysisParams* params, Double_t minJetPt, Double_t maxJetPt);
  Bool_t          PlotInvMassHistogramsVsDz(DJetCorrAnalysisParams* params);
  Bool_t          PlotInvMassHistogramArray(Int_t n, TH1** histos, const char* name, const char* xTitle, Double_t minMass, Double_t maxMass, Double_t pdgMass=-1);
  
  Bool_t          PlotObservable(DJetCorrAnalysisParams* params, TString obsName, Double_t xmin, Double_t xmax,
                                 Double_t minJetPt, Double_t maxJetPt, Double_t minDPt, Double_t maxDPt, Double_t minZ, Double_t maxZ,
                                 Int_t step, Int_t rebin, Int_t norm);
    
  Bool_t          Plot1DHistos(TString cname, Int_t n, TH1** histos, Double_t xmin, Double_t xmax);
  void            GenerateAxisMap(THnSparse* hn);
  Bool_t          OpenInputFile();
  Bool_t          LoadInputList(const char* lname);
  Bool_t          LoadQAList();
  Bool_t          LoadTHnSparse();
  Bool_t          LoadOutputHistograms();
  TVirtualPad*    SetUpPad(TVirtualPad* pad,
                           const char* xTitle, Double_t minX, Double_t maxX, Bool_t logX,
                           const char* yTitle, Double_t minY, Double_t maxY, Bool_t logY);
  TCanvas*        SetUpCanvas(const char* name,
                              const char* xTitle, Double_t minX, Double_t maxX, Bool_t logX,
                              const char* yTitle, Double_t minY, Double_t maxY, Bool_t logY,
                              Double_t w = 700, Double_t h = 500, Int_t rows=1, Int_t cols=1);
  
  TLegend*        SetUpLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t textSize);
  TPaveText*      SetUpPaveText(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t textSize, const char* text = 0);
  
  void            SavePlot(TCanvas* canvas);

  void            CloseInputFile();

  TString         fTrainName                 ;//  train name
  TString         fInputPath                 ;//  train path
  TString         fInputFileName             ;//  input file name
  TString         fQAListName                ;//  QA list name
  TString         fInputDirFileName          ;//  input dir file name
  TString         fOutputPath                ;//  output path
  TString         fOutputFileName            ;//  output file name
  Bool_t          fOverwrite                 ;//  whether the output file should be overwritten

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

  TString         fPlotFormat                ;//  plot format (pdf, eps, png...)
  Bool_t          fSavePlots                 ;//  true if save plots

  TList          *fAnalysisParams            ;//  list of analysis params
  
  TMap            fTHnAxisMap                ;//! mapping of axis titles with indexes
  Bool_t          fTHnSparseMapGenerated     ;//! whether or not the axis map has been generated for the THnSparse
  TFile          *fInputFile                 ;//! input file
  TDirectoryFile *fInputDirectoryFile        ;//! input directory file
  TList          *fInputList                 ;//! list contains the input hisograms
  TList          *fInputQAList               ;//! list contains the QA histograms
  THnSparse      *fDmesons                   ;//! THnSparse contains the results
  TList          *fOutputList                ;//! list contains the output histograms
  
  static const Double_t fgkEpsilon;

 private:
   
  DJetCorrAnalysis(const DJetCorrAnalysis &source);
  DJetCorrAnalysis& operator=(const DJetCorrAnalysis& source); 

  ClassDef(DJetCorrAnalysis, 1);
};
