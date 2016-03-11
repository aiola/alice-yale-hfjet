// Base class for D jet correlation analysis
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#ifndef DJETCORRBASE_H
#define DJETCORRBASE_H

class TString;
class TFile;
class TList;
class TCanvas;
class TLegend;
class TVirtualPad;
class TPaveText;
class TObject;
class TList;
class TDirectoryFile;
class DJetCorrAnalysisParams;

#include <TH1.h>
#include <TH2.h>
#include <TMap.h>
#include <THnSparse.h>
#include <THashList.h>
#include <TGraph.h>
#include <TPad.h>
#include <TCanvas.h>

#include "DJetCorrAnalysisParams.h"

class DJetCorrBase : public TNamed {
  
 public:
  DJetCorrBase();
  DJetCorrBase(const char* train, const char* path = "$JETRESULTS");

  void   SetInputTrain(const char* train)      { fTrainName        = train  ; CloseInputFile(); }
  void   SetInputPath(const char* path)        { fInputPath        = path   ; CloseInputFile(); }
  void   SetInputFileName(const char* fname)   { fInputFileName    = fname  ; CloseInputFile(); }
  void   SetOutputFileName(const char* fname)  { fOutputFileName   = fname  ; }
  void   SetOutputPath(const char* path)       { fOutputPath       = path   ; }
  void   SetOverwrite(Bool_t ow = kTRUE)       { fOverwrite        = ow     ; }
  void   SetPlotFormat(const char* f)          { fPlotFormat       = f      ; }
  void   SetSavePlots(Bool_t s)                { fSavePlots        = s      ; }
  void   SetAddTrainToCanvasName(Bool_t s)     { fAddTrainToCanvasName= s   ; }

  //DJetCorrAnalysisParams* AddAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius, const char* tracksName, Bool_t isMC=kFALSE, Bool_t isBkgSub=kFALSE);
  DJetCorrAnalysisParams* AddAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius, const char* trigger);
  DJetCorrAnalysisParams* AddAnalysisParams(DJetCorrAnalysisParams* params);

  const char* GetParamName(Int_t i) const;
  
  virtual Bool_t Init();
  
  virtual Bool_t SaveOutputFile();
  TString GetOutpuFileName() const;

  TH1* GetOutputHistogram(const char* name) { return fOutputList == 0 ? 0 : dynamic_cast<TH1*>(fOutputList->FindObject(name)); }
  THnSparse* GetOutputSparseHistogram(const char* name) { return fOutputList == 0 ? 0 : dynamic_cast<THnSparse*>(fOutputList->FindObject(name)); }
  TCanvas* GetCanvas(const char* name)      { return fCanvases == 0 ? 0 : static_cast<TCanvas*>(fCanvases->FindObject(name)); }

  virtual Bool_t LoadOutputHistograms();
  void UpdateAllCanvases();

  virtual Bool_t Regenerate() { return kTRUE;}

  TH2* GetTruth(Int_t p=0, Bool_t copy=kFALSE);
  TH2* GetMeasured(Int_t p=0, Bool_t copy=kFALSE);

  TH1* GetDzTruth(Int_t p, Bool_t copy);
  TH1* GetDzMeasured(Int_t p, Bool_t copy);

  TH1* GetDPtTruth(Int_t p, Bool_t copy);
  TH1* GetDPtMeasured(Int_t p, Bool_t copy);

  TH1* GetDEtaTruth(Int_t p, Bool_t copy);
  TH1* GetDEtaMeasured(Int_t p, Bool_t copy);

  TH1* GetJetPtTruth(Int_t p, Bool_t copy);
  TH1* GetJetPtMeasured(Int_t p, Bool_t copy);

  virtual TString GetTruthName(Int_t /*p*/) { return ""; }
  virtual TString GetMeasuredName(Int_t /*p*/) { return ""; }

  virtual TString GetDzTruthName(Int_t /*p*/) { return ""; }
  virtual TString GetDzMeasuredName(Int_t /*p*/) { return ""; }

  virtual TString GetDPtTruthName(Int_t /*p*/) { return ""; }
  virtual TString GetDPtMeasuredName(Int_t /*p*/) { return ""; }

  virtual TString GetDEtaTruthName(Int_t /*p*/) { return ""; }
  virtual TString GetDEtaMeasuredName(Int_t /*p*/) { return ""; }

  virtual TString GetJetPtTruthName(Int_t /*p*/) { return ""; }
  virtual TString GetJetPtMeasuredName(Int_t /*p*/) { return ""; }

  static void FitObjectInPad(TObject* obj, TVirtualPad* pad, Option_t* opt="", Bool_t copyAxisTitle=kFALSE, Double_t extraFactor=1.8);
  static void FitGraphInPad(TGraph* graph, TVirtualPad* pad, Option_t* opt="", Bool_t copyAxisTitle=kFALSE, Double_t extraFactor=1.8);
  static void FitHistogramInPad(TH1* hist, TVirtualPad* pad, Option_t* opt="", Bool_t copyAxisTitle=kFALSE, Double_t extraFactor=1.8);
  static void GetMinMax(TGraph* graph, Double_t& miny, Double_t& maxy);
  static void GetMinMax(TH1* hist, Double_t& miny, Double_t& maxy);
  static TLegend* GetLegend(TPad* pad);

  static THnSparse* Rebin(THnSparse* orig, const char* name, const Int_t* nbins, const Double_t** bins);
  static TH1* Rebin(TH1* orig, const char* name, Int_t nbins, const Double_t* bins);
  static TH2* Rebin(TH2* orig, const char* name, Int_t nbinsx, const Double_t* binsx, Int_t nbinsy, const Double_t* binsy);
  static TH1* Rebin(TH1* orig, const char* name, Int_t nbins, Double_t min, Double_t max);
  static TH2* Rebin(TH2* orig, const char* name, Int_t nbinsx, Double_t minX, Double_t maxX, Int_t nbinsy, Double_t minY, Double_t maxY);
  static Bool_t CheckExactRebin(TAxis* orig, TAxis* dest);
  static void GetBinCenter(THnSparse* hn, Int_t* coord_ind, Double_t* coord);
  static void MakeBinomialConsistent(TH1* pass, TH1* total);
  static Double_t* GenerateFixedArray(Int_t n, Double_t min, Double_t max);
  static TH2* Normalize(TH2* orig, const char* name);
  static TMap* GenerateAxisMap(THnSparse* hn);

 protected:

  virtual Bool_t  SaveOutputFile(TObjArray& arr);

  Bool_t          GenerateRatios(const char* nname, const char* dname);

  Bool_t          PlotObservable(DJetCorrAnalysisParams* params, TString obsName, Double_t xmin, Double_t xmax,
                                 Double_t minDPt, Double_t maxDPt, Double_t minJetPt, Double_t maxJetPt, Double_t minZ, Double_t maxZ,
                                 Int_t step, Int_t rebin, Int_t norm, Int_t plotStats=0);
    
  Bool_t          Plot1DHistos(TString cname, TObjArray& histos, Double_t xmin, Double_t xmax, Int_t plotStats=0);
  Bool_t          OpenInputFile();
  TFile*          OpenOutputFile();
  virtual Bool_t  LoadTHnSparse() {return kTRUE;}

  TVirtualPad*    SetUpPad(TVirtualPad* pad,
                           const char* xTitle, Double_t minX, Double_t maxX, Bool_t logX,
                           const char* yTitle, Double_t minY, Double_t maxY, Bool_t logY,
                           Double_t lmar=0.12, Double_t rmar=0.08, Double_t bmar=0.12, Double_t tmar=0.08);
  TCanvas*        SetUpCanvas(const char* name,
                              const char* xTitle, Double_t minX, Double_t maxX, Bool_t logX,
                              const char* yTitle, Double_t minY, Double_t maxY, Bool_t logY,
                              Double_t w = 700, Double_t h = 500, Int_t rows=1, Int_t cols=1,
                              Double_t lmar=0.12, Double_t rmar=0.08, Double_t bmar=0.12, Double_t tmar=0.08);
  TCanvas*        SetUpCanvas(TH1* histo, Bool_t logX, Bool_t logY,
                              Double_t w = 700, Double_t h = 500, Int_t rows=1, Int_t cols=1,
                              Double_t lmar=0.12, Double_t rmar=0.08, Double_t bmar=0.12, Double_t tmar=0.08);
  
  TLegend*        SetUpLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t textSize);
  TPaveText*      SetUpPaveText(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t textSize, const char* text = 0);
  
  virtual void    SavePlot(TCanvas* canvas);

  void            CloseInputFile();

  TString         fTrainName                 ;//  train name
  TString         fInputPath                 ;//  train path
  TString         fInputFileName             ;//  input file name
  TString         fOutputPath                ;//  output path
  TString         fOutputFileName            ;//  output file name
  Bool_t          fOverwrite                 ;//  whether the output file should be overwritten

  TList          *fAnalysisParams            ;//  list of analysis params
  
  TString         fPlotFormat                ;//  plot format (pdf, eps, png...)
  Bool_t          fSavePlots                 ;//  true if save plots
  Bool_t          fAddTrainToCanvasName      ;//  automatically add the train string to each canvas' name

  DJetCorrAnalysisParams::DJetCorrAnalysisType fAnaType;//  analysis type
  
  TFile          *fInputFile                 ;//! input file
  THashList      *fOutputList                ;//! list contains the output histograms
  TList          *fCanvases                  ;//! list contains the canvases
  
  static const Double_t fgkEpsilon;

 private:
   
  DJetCorrBase(const DJetCorrBase &source);
  DJetCorrBase& operator=(const DJetCorrBase& source); 

  ClassDef(DJetCorrBase, 1);
};

#endif
