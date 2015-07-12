// Base class for D jet correlation analysis
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
class DJetCorrAnalysisParams;

#include <TMap.h>

class DJetCorrBase : public TObject {
  
 public:
  DJetCorrBase();
  DJetCorrBase(const char* train, const char* path = "$JETRESULTS");

  enum EMatchingStatus {kAnyMatchingStatus, kNotMatched, kMatched};

  void   SetInputTrain(const char* train)      { fTrainName        = train  ; CloseInputFile(); }
  void   SetInputPath(const char* path)        { fInputPath        = path   ; CloseInputFile(); }
  void   SetInputFileName(const char* fname)   { fInputFileName    = fname  ; CloseInputFile(); }
  void   SetOutputFileName(const char* fname)  { fOutputFileName   = fname  ; }
  void   SetOutputPath(const char* path)       { fOutputPath       = path   ; }
  void   SetOverwrite(Bool_t ow = kTRUE)       { fOverwrite        = ow     ; }
  void   SetPlotFormat(const char* f)          { fPlotFormat       = f      ; }
  void   SetSavePlots(Bool_t s)                { fSavePlots        = s      ; }

  DJetCorrAnalysisParams*   AddAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius, const char* tracksName);
  DJetCorrAnalysisParams*   AddAnalysisParams(DJetCorrAnalysisParams* params);
  
  virtual Bool_t Init();
  
  Bool_t SaveOutputFile();
  
 protected:

  virtual void    ClearInputData();
  Int_t           GetAxisIndex(TString title, THnSparse* hn);

  Bool_t          GenerateRatios(const char* nname, const char* dname);

  Bool_t          PlotObservable(DJetCorrAnalysisParams* params, TString obsName, Double_t xmin, Double_t xmax,
                                 Double_t minDPt, Double_t maxDPt, Double_t minJetPt, Double_t maxJetPt, Double_t minZ, Double_t maxZ,
                                 Int_t step, Int_t rebin, Int_t norm, Int_t plotStats=0);
    
  Bool_t          Plot1DHistos(TString cname, Int_t n, TH1** histos, Double_t xmin, Double_t xmax, Int_t plotStats=0);
  TMap*           GenerateAxisMap(THnSparse* hn);
  Bool_t          OpenInputFile();
  Bool_t          LoadInputList(const char* lname);
  virtual Bool_t  LoadTHnSparse()=0;
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
  TString         fInputDirFileName          ;//  input dir file name
  TString         fOutputPath                ;//  output path
  TString         fOutputFileName            ;//  output file name
  Bool_t          fOverwrite                 ;//  whether the output file should be overwritten

  TList          *fAnalysisParams            ;//  list of analysis params
  
  TString         fPlotFormat                ;//  plot format (pdf, eps, png...)
  Bool_t          fSavePlots                 ;//  true if save plots
  
  TList           fTHnSparseAxisMaps         ;//! axis maps for THnSparse
  TFile          *fInputFile                 ;//! input file
  TDirectoryFile *fInputDirectoryFile        ;//! input directory file
  TList          *fInputList                 ;//! list contains the input hisograms
  TList          *fOutputList                ;//! list contains the output histograms
  
  static const Double_t fgkEpsilon;

 private:
   
  DJetCorrBase(const DJetCorrBase &source);
  DJetCorrBase& operator=(const DJetCorrBase& source); 

  ClassDef(DJetCorrBase, 1);
};
