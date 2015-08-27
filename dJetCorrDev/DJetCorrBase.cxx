// Base class for D jet correlation analysis
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <algorithm>

#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <Riostream.h>
#include <TParameter.h>
#include <TObjString.h>
#include <TString.h>
#include <TDirectoryFile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>

#include "DJetCorrAnalysisParams.h"

#include "DJetCorrBase.h"

const Double_t DJetCorrBase::fgkEpsilon = 1E-6;

ClassImp(DJetCorrBase);

//____________________________________________________________________________________
DJetCorrBase::DJetCorrBase() :
  TNamed(),
  fTrainName(),
  fInputPath(),
  fInputFileName(),
  fInputDirFileName(),
  fOutputPath(),
  fOutputFileName(),
  fOverwrite(kFALSE),
  fAnalysisParams(0),
  fPlotFormat(),
  fSavePlots(kFALSE),
  fAddTrainToCanvasName(kFALSE),
  fAnaType(DJetCorrAnalysisParams::KUndefinedAna),
  fTHnSparseAxisMaps(),
  fInputFile(0),
  fInputDirectoryFile(0),
  fInputList(0),
  fOutputList(0),
  fCanvases(0),
  fEvents(0.)
{
  // Default ctr.
  
}

//____________________________________________________________________________________
DJetCorrBase::DJetCorrBase(const char* train, const char* path) :
  TNamed(train, train),
  fTrainName(train),
  fInputPath(path),
  fInputFileName(),
  fInputDirFileName(),
  fOutputPath("../data/"),
  fOutputFileName("DJetCorr.root"),
  fOverwrite(kFALSE),
  fAnalysisParams(new TList()),
  fPlotFormat("pdf"),
  fSavePlots(kFALSE),
  fAddTrainToCanvasName(kFALSE),
  fAnaType(DJetCorrAnalysisParams::KUndefinedAna),
  fTHnSparseAxisMaps(),
  fInputFile(0),
  fInputDirectoryFile(0),
  fInputList(0),
  fOutputList(0),
  fCanvases(0),
  fEvents(0.)
{
  // Standard ctr.

}

//____________________________________________________________________________________
Bool_t DJetCorrBase::ClearInputData()
{
  // Clear the input data.

  fTHnSparseAxisMaps.Clear();

  return kTRUE;
}

//____________________________________________________________________________________
DJetCorrAnalysisParams* DJetCorrBase::AddAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius, const char* tracksName, Bool_t isMC)
{
  // Add the analysis params for the D meson type, jet type and jet radius provided.

  if (!fAnalysisParams) fAnalysisParams = new TList();

  DJetCorrAnalysisParams* params = new DJetCorrAnalysisParams(dmeson, jetType, jetRadius, tracksName, fAnaType, isMC);

  fAnalysisParams->Add(params);

  return params;
}

//____________________________________________________________________________________
DJetCorrAnalysisParams* DJetCorrBase::AddAnalysisParams(DJetCorrAnalysisParams* params)
{
  // Add the analysis params.

  if (!fAnalysisParams) fAnalysisParams = new TList();

  DJetCorrAnalysisParams* paramsCopy = new DJetCorrAnalysisParams(*params);
  
  fAnalysisParams->Add(paramsCopy);

  return paramsCopy;
}

//____________________________________________________________________________________
TMap* DJetCorrBase::GenerateAxisMap(THnSparse* hn)
{
  if (!hn) return 0;
  
  TMap* axisMap = new TMap();
  axisMap->SetOwnerKeyValue();
  
  for (Int_t i = 0; i < hn->GetNdimensions(); i++) {
    TObjString* key = new TObjString(hn->GetAxis(i)->GetTitle());
    TParameter<Int_t>* value = new TParameter<Int_t>(hn->GetAxis(i)->GetTitle(), i);
    axisMap->Add(key, value);
  }

  fTHnSparseAxisMaps.Add(axisMap);
  
  return axisMap;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::Init()
{
  // Init class.

  Printf("Info-DJetCorrBase::Init : Now initializing.");
  
  if (fOutputList) {
    delete fOutputList;
    fOutputList = 0;
  }
  fOutputList = new TList();

  if (fCanvases) {
    delete fCanvases;
    fCanvases = 0;
  }
  fCanvases = new TList();  
    
  Printf("Info-DJetCorrBase::Init : Initialization done.");

  return kTRUE;
}

//____________________________________________________________________________________
TFile* DJetCorrBase::OpenOutputFile()
{
  TString fname(fOutputPath);
  
  if (fTrainName.IsNull()) {
    fname += "/test/";
  }
  else {
    fname += "/";
    fname += fTrainName;
    fname += "/";
  }

  fname += fOutputFileName;
  
  TFile* outputFile = TFile::Open(fname);

  if (!outputFile || outputFile->IsZombie()) {
    Printf("Error-DJetCorrBase::LoadOutputHistograms : Could not open file '%s' to read.", fname.Data()); 
    outputFile = 0;
  }

  return outputFile;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::LoadOutputHistograms()
{
  TFile* outputFile = OpenOutputFile();

  if (!outputFile) return kFALSE;

  fOutputList = static_cast<TList*>(outputFile->Get("fOutputList"));

  outputFile->Close();
  delete outputFile;
  outputFile = 0;

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrBase::SavePlot(TCanvas* canvas)
{
  if (!canvas) return;
  
  TString fname(fOutputPath);

  if (fTrainName.IsNull()) {
    fname += "/test/";
  }
  else {
    fname += "/";
    fname += fTrainName;
    fname += "/";
  }

  fname += canvas->GetTitle();
  fname += ".";
  fname += fPlotFormat;

  canvas->SaveAs(fname);
}

//____________________________________________________________________________________
TVirtualPad* DJetCorrBase::SetUpPad(TVirtualPad* pad,
                                    const char* xTitle, Double_t minX, Double_t maxX, Bool_t logX,
                                    const char* yTitle, Double_t minY, Double_t maxY, Bool_t logY,
                                    Double_t lmar, Double_t rmar, Double_t bmar, Double_t tmar)
{
  if (!pad) return 0;

  pad->SetLeftMargin(lmar);
  pad->SetRightMargin(rmar);
  pad->SetBottomMargin(bmar);
  pad->SetTopMargin(tmar);
  
  pad->SetLogx(logX);
  pad->SetLogy(logY);
  pad->cd();

  TString blankHistName(Form("%s_blankHist", pad->GetName()));
  TH1* blankHist = new TH1D(blankHistName, blankHistName, 1000, minX, maxX);

  blankHist->SetMinimum(minY);
  blankHist->SetMaximum(maxY);
  
  blankHist->GetXaxis()->SetTitle(xTitle);
  blankHist->GetYaxis()->SetTitle(yTitle);
  
  blankHist->Draw("AXIS");

  return pad;
}

//____________________________________________________________________________________
TCanvas* DJetCorrBase::SetUpCanvas(TH1* histo, Bool_t logX, Bool_t logY,
                                   Double_t w, Double_t h, Int_t rows, Int_t cols,
                                   Double_t lmar, Double_t rmar, Double_t bmar, Double_t tmar)
{
  return SetUpCanvas(histo->GetName(),
                     histo->GetXaxis()->GetTitle(), histo->GetXaxis()->GetXmin(), histo->GetXaxis()->GetXmax(), logX,
                     histo->GetYaxis()->GetTitle(), histo->GetYaxis()->GetXmin(), histo->GetYaxis()->GetXmax(), logY,
                     w, h, rows, cols, lmar, rmar, bmar, tmar);
}

//____________________________________________________________________________________
TCanvas* DJetCorrBase::SetUpCanvas(const char* name,
                                   const char* xTitle, Double_t minX, Double_t maxX, Bool_t logX,
                                   const char* yTitle, Double_t minY, Double_t maxY, Bool_t logY,
                                   Double_t w, Double_t h, Int_t rows, Int_t cols,
                                   Double_t lmar, Double_t rmar, Double_t bmar, Double_t tmar)
{
  Printf("Info-DJetCorrBase::SetUpCanvas : Setting up canvas '%s'", name);

  TString cname;

  if (fAddTrainToCanvasName) {
    cname = Form("%s_%s", GetName(), name);
  }
  else {
    cname = name;
  }

  TCanvas* canvas = static_cast<TCanvas*>(fCanvases->FindObject(cname));
  if (canvas) {
    Printf("Warning-DJetCorrBase::SetUpCanvas : Canvas '%s' already exists.", cname.Data());
    return canvas;
  }
  
  canvas = new TCanvas(cname, cname, w, h);
  
  if (rows == 1 && cols == 1) {
    SetUpPad(canvas, xTitle, minX, maxX, logX, yTitle, minY, maxY, logY);    
  }
  else {
    canvas->Divide(cols, rows);
    Int_t n = rows * cols;
    for (Int_t i = 1; i <= n; i++) {
      TVirtualPad* pad = canvas->cd(i);
      SetUpPad(pad, xTitle, minX, maxX, logX, yTitle, minY, maxY, logY,
               lmar, rmar, bmar, tmar); 
    }
  }

  fCanvases->Add(canvas);
  
  return canvas;
}

//____________________________________________________________________________________
TLegend* DJetCorrBase::SetUpLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t textSize)
{
  TLegend* leg = new TLegend(x1, y1, x2, y2);
  leg->SetTextFont(43);
  leg->SetTextSize(textSize);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  return leg;
}

//____________________________________________________________________________________
TPaveText* DJetCorrBase::SetUpPaveText(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t textSize, const char* text)
{
  TPaveText* pave = new TPaveText(x1, y1, x2, y2, "brNDC");
  pave->SetTextFont(43);
  pave->SetTextAlign(11);
  pave->SetTextSize(textSize);
  pave->SetFillStyle(0);
  pave->SetBorderSize(0);
  
  if (text) {
    TString stext(text);
    TObjArray* lines = stext.Tokenize("\n");
    TIter next(lines);
    TObject *obj = 0;
    while ((obj = next())) {
      pave->AddText(obj->GetName());
    }
    delete lines;
  }

  return pave;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::PlotObservable(DJetCorrAnalysisParams* params, TString obsName, Double_t xmin, Double_t xmax,
                                        Double_t minDPt, Double_t maxDPt, Double_t minJetPt, Double_t maxJetPt, Double_t minZ, Double_t maxZ,
                                        Int_t step, Int_t rebin, Int_t norm, Int_t plotStats)
{
  if (!fOutputList) return kFALSE;

  TString cname(Form("fig_%s_%s",obsName.Data(), params->GetName()));
  
  TString jetCuts;
  TString dCuts;
  TString zCuts;

  Int_t nj = params->GetNJetPtBins();
  Int_t nd = params->GetNDPtBins();
  Int_t nz = params->GetNzBins();

  if (maxJetPt > 0) {
    jetCuts = Form("JetPt_%03.0f_%03.0f", minJetPt, maxJetPt);
    jetCuts.ReplaceAll(".", "");
    nj = 1;
  }

  if (maxDPt > 0) {
    dCuts = Form("DPt_%02.0f_%02.0f", minDPt, maxDPt);
    dCuts.ReplaceAll(".", "");
    nd = 1;
  }

  if (maxZ > 0) {
    zCuts = Form("z_%.1f_%.1f", minZ, maxZ);
    zCuts.ReplaceAll(".", "");
    nz = 1;
  }
  
  if (!jetCuts.IsNull()) {
    cname += "_";
    cname += jetCuts;
  }

  if (!dCuts.IsNull()) {
    cname += "_";
    cname += dCuts;
  }

  if (!zCuts.IsNull()) {
    cname += "_";
    cname += zCuts;
  }

  TString hname(obsName);
  TString prefix(params->GetName());

  TH1** histos = new TH1*[nj*nd*nz];
  Int_t ih = 0;
  TString jetTitle;
  TString dTitle;
  TString zTitle;
  
  for (Int_t ij = 0; ij < nj; ij+=step) {
    if (maxJetPt < 0) {
      jetCuts = Form("JetPt_%03.0f_%03.0f", params->GetJetPtBin(ij), params->GetJetPtBin(ij+1));
      jetCuts.ReplaceAll(".", "");
      jetTitle = Form("%.1f < #it{p}_{T,jet} < %.1f GeV/#it{c}", params->GetJetPtBin(ij), params->GetJetPtBin(ij+1));
    }
    for (Int_t id = 0; id < nd; id+=step) {
      if (maxDPt < 0) {
        dCuts = Form("DPt_%02.0f_%02.0f", params->GetDPtBin(id), params->GetDPtBin(id+1));
        dCuts.ReplaceAll(".", "");
        dTitle = Form("%.1f < #it{p}_{T,D} < %.1f GeV/#it{c}", params->GetDPtBin(id), params->GetDPtBin(id+1));
      }
      for (Int_t iz = 0; iz < nz; iz+=step) {
        if (maxZ < 0) {
          zCuts = Form("z_%.1f_%.1f", params->GetzBin(iz), params->GetzBin(iz+1));
          zCuts.ReplaceAll(".", "");
          zTitle = Form("%.1f < z < %.1f", params->GetzBin(iz), params->GetzBin(iz+1));
        }
            
        TString cuts(Form("%s_%s_%s", dCuts.Data(), jetCuts.Data(), zCuts.Data()));
        TString htitle;

        htitle = jetTitle;
        if (!dTitle.IsNull()) {
          if (!htitle.IsNull()) htitle += ", ";
          htitle += dTitle;
        }
        if (!zTitle.IsNull()) {
          if (!htitle.IsNull()) htitle += ", ";
          htitle += zTitle;
        }
    
        TString objname(Form("h%s_%s_%s_Matched", prefix.Data(), hname.Data(), cuts.Data()));
        //Printf("Info-DJetCorrAnalysis::PlotObservable : Retrieving histogram '%s'", objname.Data());
        TH1* hist = dynamic_cast<TH1*>(fOutputList->FindObject(objname));
        if (!hist) {
          Printf("Error-DJetCorrAnalysis::PlotObservable : Histogram '%s' not found!", objname.Data());
          continue;
        }
        TString newName(objname);
        newName += "_copy";
        //Printf("Info-DJetCorrAnalysis::PlotObservable : Cloning histogram '%s'", objname.Data());
        histos[ih] = static_cast<TH1*>(hist->Clone(newName));
        if (norm == 1) {
          histos[ih]->Scale(1. / histos[ih]->Integral(), "width");
          histos[ih]->GetYaxis()->SetTitle("Probability density");
        }
        else if (norm == 2) {
          histos[ih]->Scale(1. / histos[ih]->Integral(), "");
          TF1 f1("f1", "TMath::TwoPi() * x", 0, histos[ih]->GetXaxis()->GetXmax()*1.5);
          for (Int_t ibin = 1; ibin <= histos[ih]->GetNbinsX(); ibin++) {
            Double_t integ = f1.Integral(histos[ih]->GetBinLowEdge(ibin), histos[ih]->GetBinLowEdge(ibin+1));
            histos[ih]->SetBinContent(ibin, histos[ih]->GetBinContent(ibin) / integ);
            histos[ih]->SetBinError(ibin, histos[ih]->GetBinError(ibin) / integ);
          }
          histos[ih]->GetYaxis()->SetTitle("Prob. density / 2#pi#Delta R");
        }
        if (rebin > 1) histos[ih]->Rebin(4);
        histos[ih]->SetTitle(htitle);
        ih++;
      }
    }
  }

  if (ih == 0) return kFALSE;
  
  return Plot1DHistos(cname, ih, histos, xmin, xmax, plotStats);
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::Plot1DHistos(TString cname, Int_t n, TH1** histos, Double_t xmin, Double_t xmax, Int_t plotStats)
{
  if (!histos[0]) return kFALSE;
  
  TString xTitle(histos[0]->GetXaxis()->GetTitle());
  TString yTitle(histos[0]->GetYaxis()->GetTitle());
  if (xmax - xmin < fgkEpsilon) {
    xmin = histos[0]->GetXaxis()->GetXmin();
    xmax = histos[0]->GetXaxis()->GetXmax();
  }
  TCanvas *canvas = SetUpCanvas(cname, xTitle, xmin, xmax, kFALSE, yTitle, 0, 1, kFALSE);
  TH1* hblank = dynamic_cast<TH1*>(canvas->GetListOfPrimitives()->At(0));
  Color_t colors[12] = {kRed+1, kBlue+1, kMagenta, kGreen+2, kOrange+1, kCyan+2, kPink+1, kTeal+1, kViolet, kAzure+1, kYellow+2, kSpring+3};

  TLegend* leg = 0;
  Int_t neff = n;
  if (plotStats) neff *= 2;
  if (neff > 8) leg = SetUpLegend(0.38, 0.45, 0.88, 0.88, 15);
  else leg = SetUpLegend(0.38, 0.69, 0.88, 0.88, 15);
  leg->SetNColumns(2);
  Double_t maxY = 0;
  for (Int_t i = 0; i < n; i++) {
    if (maxY < histos[i]->GetMaximum()) maxY = histos[i]->GetMaximum();
    histos[i]->SetLineColor(colors[i]);
    histos[i]->SetMarkerColor(colors[i]);
    histos[i]->SetMarkerStyle(kFullCircle);
    histos[i]->SetMarkerSize(0.5);
    histos[i]->Draw("same");
    TLegendEntry* legEntry = leg->AddEntry(histos[i], histos[i]->GetTitle(), "pe");
    legEntry->SetLineColor(colors[i]);
    if (plotStats) {
      TString statsStr(Form("#mu = %.3f, #sigma = %.3f", histos[i]->GetMean(), histos[i]->GetRMS()));
      leg->AddEntry((TObject*)0, statsStr, "");
    }
  }

  if (hblank && maxY > 0) {
    hblank->GetYaxis()->SetRangeUser(0, maxY*1.6);
  }

  leg->Draw();

  if (fSavePlots) {
    SavePlot(canvas);
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Int_t DJetCorrBase::GetAxisIndex(TString title, THnSparse* hn, Bool_t messageOnFail)
{
  TMap* axisMap = static_cast<TMap*>(fTHnSparseAxisMaps.FindObject(hn->GetName()));
  if (!axisMap) {
    axisMap = GenerateAxisMap(hn);
    if (!axisMap) return -1;
  }

  TParameter<Int_t>* par = static_cast<TParameter<Int_t>*>(axisMap->GetValue(title));
  if (!par) {
    if (messageOnFail) Printf("Warning-DJetCorrBase::GetAxisIndex : Could not find axis with title '%s' in histogram '%s'.", title.Data(), hn->GetName());
    return -1;
  }
  
  return par->GetVal();
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::GenerateRatios(const char* nname, const char* dname)
{
  TIter next(fOutputList);
  TObject* obj = 0;

  while ((obj = next())) {
    TH1* hist1D = dynamic_cast<TH1*>(obj);
    if (!hist1D) continue;

    TString hdname(hist1D->GetName());
    if (!hdname.Contains(dname)) continue;

    TString hnname(hdname);
    hnname.ReplaceAll(dname, nname);

    TH1* num = dynamic_cast<TH1*>(fOutputList->FindObject(hnname));

    if (!num) continue;

    TString rname(Form("Ratio_%s_%s", hnname.Data(), hdname.Data()));

    if (fOutputList->Contains(rname)) continue;

    Printf("Info-DJetCorrBase::GenerateRatios : Now calcutaling ratio '%s' over '%s'", hnname.Data(), hdname.Data());

    if (hist1D->GetSumw2N() == 0) hist1D->Sumw2();
    if (num->GetSumw2N() == 0) num->Sumw2();

    TH2* hist2D = dynamic_cast<TH2*>(hist1D);
    if (hist2D) {
      TH2* num2D = dynamic_cast<TH2*>(num);
      if (!num2D) continue;
      TH2* ratio2D = static_cast<TH2*>(num2D->Clone(rname));
      ratio2D->Divide(hist2D);
      fOutputList->Add(ratio2D);
    }
    else {
      TGraphAsymmErrors* graph = new TGraphAsymmErrors(num, hist1D);
      graph->SetName(rname);
      fOutputList->Add(graph);
    }
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::OpenInputFile()
{
  // Open the input file and set the input list.
  
  TString fname;
  
  fname = fInputPath;
  fname += "/";
  fname += fTrainName;
  fname += "/";
  fname += fInputFileName;

  if (fInputDirectoryFile == fInputFile) {
    fInputDirectoryFile = 0; // to avoid double delete
  }

  if (!fInputFile) {
    Printf("Info-DJetCorrBase::OpenInputFile : Opening file '%s'", fname.Data()); 
    fInputFile = TFile::Open(fname);
  }
  else {
    Printf("Info-DJetCorrBase::OpenInputFile : File '%s' already open", fInputFile->GetName()); 
  }

  if (!fInputFile || fInputFile->IsZombie()) {
    Printf("Error-DJetCorrBase::OpenInputFile : Could not open file '%s'", fname.Data()); 
    fInputFile = 0;
    return kFALSE;
  }

  if (!fInputDirectoryFile || fInputDirFileName != fInputDirectoryFile->GetName()) {
    if (fInputDirFileName.IsNull()) {
      fInputDirectoryFile = fInputFile;
    }
    else {
      ClearInputData();
    
      delete fInputDirectoryFile;
      Printf("Info-DJetCorrBase::OpenInputFile : Getting directory '%s' from file '%s'", fInputDirFileName.Data(), fInputFile->GetName()); 
      fInputDirectoryFile = dynamic_cast<TDirectoryFile*>(fInputFile->Get(fInputDirFileName));

      delete fInputList;
      fInputList = 0;

      if (!fInputDirectoryFile) {
        Printf("Error-DJetCorrBase::OpenInputFile : Could not get directory '%s' from file '%s'", fInputDirFileName.Data(), fInputFile->GetName()); 
        return kFALSE;
      }
    }
  }

  Printf("Info-DJetCorrBase::OpenInputFile : Success.");

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::LoadInputList(const char* inputListName)
{
  // Load the input list.

  if (!OpenInputFile()) return kFALSE;
  
  if (!fInputList || strcmp(inputListName, fInputList->GetName()) != 0) {
    ClearInputData();
    
    delete fInputList;
    Printf("Info-DJetCorrBase::OpenInputFile : Getting list '%s' from directory '%s' of file '%s'", inputListName, fInputDirectoryFile->GetName(), fInputFile->GetName()); 
    fInputList = dynamic_cast<TList*>(fInputDirectoryFile->Get(inputListName));
  }
  
  if (!fInputList) {
    Printf("Error-DJetCorrBase::OpenInputFile : Could not get list '%s' from directory '%s' of file '%s'", inputListName, fInputDirectoryFile->GetName(), fInputFile->GetName());
    return kFALSE;
  }

  Printf("Info-DJetCorrBase::LoadInputList : Success.");

  GetEvents();

  Printf("Info-DJetCorrBase::LoadInputList : Total number of events: %.0f.", fEvents);

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrBase::CloseInputFile()
{
  // Close the input file.

  ClearInputData();
  
  if (fInputFile) {
    fInputFile->Close();
    delete fInputFile;
    fInputFile = 0;
  }

  if (fInputDirectoryFile) {
    delete fInputDirectoryFile;
    fInputDirectoryFile = 0;
  }

  if (fInputList) {
    delete fInputList;
    fInputList = 0;
  }
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::SaveOutputFile()
{
  // Save the output in a file.

  TObjArray arr;
  arr.SetOwner(kFALSE);

  fOutputList->SetName("fOutputList");
    
  arr.Add(fOutputList);

  return SaveOutputFile(arr);
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::SaveOutputFile(TObjArray& arr)
{
  // Save the content arr in a file.

  TString opt("create");
  if (fOverwrite) opt = "recreate";
  
  TString fname(fOutputPath);

  if (fTrainName.IsNull()) {
    fname += "/test/";
  }
  else {
    fname += "/";
    fname += fTrainName;
    fname += "/";
  }
  
  fname += fOutputFileName;

  TString path(fname);
  path.Remove(path.Last('/'));
  if (gSystem->AccessPathName(path)) {
    Printf("Info-DJetCorrBase::SaveOutputFile : creating directory '%s'.", path.Data()); 
    gSystem->mkdir(path.Data(), kTRUE);
  }
  
  TFile* outputFile = TFile::Open(fname, opt);

  if (!outputFile || outputFile->IsZombie()) {
    Printf("Error-DJetCorrBase::SaveOutputFile : could not open file '%s' to write (if file exists, you may need to set the overwrite option).", fname.Data()); 
    outputFile = 0;
    return kFALSE;
  }

  outputFile->cd();

  Printf("Info-DJetCorrBase::SaveOutputFile : Now streaming results.");
  TIter next(&arr);
  TObject* obj = 0;
  while ((obj = next())) {
    TCollection* coll = dynamic_cast<TCollection*>(obj);
    if (coll) {
      coll->Write(coll->GetName(), TObject::kSingleKey);
    }
    else {
      obj->Write();
    }
  }

  Printf("Info-DJetCorrBase::SaveOutputFile : Closing the output file."); 
  outputFile->Close();
  delete outputFile;
  outputFile = 0;
  
  return kTRUE;
}

//____________________________________________________________________________________
Double_t DJetCorrBase::GetEvents(Bool_t recalculate)
{
  if (fEvents == 0 || recalculate) {
    fEvents = 0;

    if (fOutputList) {

      TH1* hevents = dynamic_cast<TH1*>(fOutputList->FindObject("hEvents"));

      if (!hevents && fInputList) {
        TH1* hevents_temp = static_cast<TH1*>(fInputList->FindObject("fHistEventCount"));
        if (hevents_temp) {
          hevents = static_cast<TH1*>(hevents_temp->Clone("hEvents"));
          hevents->SetTitle("hEvents");
          fOutputList->Add(hevents);
        }
      }
      
      if (hevents) {
        fEvents = hevents->GetBinContent(1);
      }
    }
  }

  return fEvents;
}

//____________________________________________________________________________________
void DJetCorrBase::FitGraphInPad(TGraph* graph, TVirtualPad* pad)
{
  TH1* blankHist = dynamic_cast<TH1*>(pad->GetListOfPrimitives()->At(0));
  if (blankHist) {
    Double_t miny = blankHist->GetMinimum();
    Double_t maxy = blankHist->GetMaximum() / 1.8;
  
    GetMinMax(graph, miny, maxy);
    
    blankHist->SetMinimum(miny);
    blankHist->SetMaximum(maxy * 1.8);
  }
  else {
    Printf("Error-DJetCorrBase::FitGraphInPad : Could not find blank histogram!");
  }
}

//____________________________________________________________________________________
void DJetCorrBase::FitHistogramInPad(TH1* hist, TVirtualPad* pad)
{
  TH1* blankHist = dynamic_cast<TH1*>(pad->GetListOfPrimitives()->At(0));
  if (blankHist) {
    Double_t miny = blankHist->GetMinimum();
    Double_t maxy = blankHist->GetMaximum() / 1.8;

    GetMinMax(hist, miny, maxy);

    blankHist->SetMinimum(miny);
    blankHist->SetMaximum(maxy * 1.8);

    hist->Draw("same");
  }
  else {
    Printf("Error-DJetCorrBase::FitHistogramInPad : Could not find blank histogram!");
  }
}

//____________________________________________________________________________________
void DJetCorrBase::GetMinMax(TGraph* graph, Double_t& miny, Double_t& maxy)
{
  Double_t* array = graph->GetY();

  for (Int_t i = 0; i < graph->GetN(); i++) {
    if (miny > graph->GetY()[i] - graph->GetEYlow()[i]) miny = graph->GetY()[i] - graph->GetEYlow()[i];
    if (maxy < graph->GetY()[i] + graph->GetEYhigh()[i]) maxy = graph->GetY()[i] + graph->GetEYhigh()[i];
  } 
}

//____________________________________________________________________________________
void DJetCorrBase::GetMinMax(TH1* hist, Double_t& miny, Double_t& maxy)
{
  Int_t minBin = hist->GetMinimumBin();
  miny = TMath::Min(hist->GetBinContent(minBin) - hist->GetBinError(minBin), miny);

  Int_t maxBin = hist->GetMaximumBin();
  maxy = TMath::Max(hist->GetBinContent(maxBin) + hist->GetBinError(maxBin), maxy);
}

//____________________________________________________________________________________
TLegend* DJetCorrBase::GetLegend(TPad* pad)
{
  TLegend* leg = 0;

  for (Int_t i = 0; i < pad->GetListOfPrimitives()->GetEntries(); i++) {
    leg = dynamic_cast<TLegend*>(pad->GetListOfPrimitives()->At(i));
    if (leg) break;
  }

  return leg;
}

//____________________________________________________________________________________
Bool_t DJetCorrBase::CheckExactRebin(TAxis* orig, TAxis* dest)
{
  for (Int_t i = 0; i <= orig->GetNbins(); i++) {
    Double_t xlow    = orig->GetBinLowEdge(i);
    Double_t xup     = orig->GetBinUpEdge(i);
    Int_t    xlowBin = dest->FindBin(xlow);
    Int_t    xupBin  = dest->FindBin(xup);
    if (TMath::Abs(xup - dest->GetBinLowEdge(xupBin)) < fgkEpsilon) {
      xupBin--;
    }
    if (xlowBin != xupBin) {
      Printf("Bin [%.2f, %.2f] -> [%.2f, %.2f] , [%.2f, %.2f]", xlow, xup, dest->GetBinLowEdge(xlowBin), dest->GetBinUpEdge(xlowBin), dest->GetBinLowEdge(xupBin), dest->GetBinUpEdge(xupBin));
      return kFALSE;
    }
  }

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrBase::GetBinCenter(THnSparse* hn, Int_t* coord_ind, Double_t* coord)
{
  for (Int_t idim = 0; idim < hn->GetNdimensions(); idim++) {
    coord[idim] = hn->GetAxis(idim)->GetBinCenter(coord_ind[idim]);
  }
}

//____________________________________________________________________________________
THnSparse* DJetCorrBase::Rebin(THnSparse* orig, const char* name, const Int_t* nbins, const Double_t** bins)
{
  THnSparse* result = new THnSparseD(name, orig->GetTitle(), orig->GetNdimensions(), nbins);
  for (Int_t idim = 0; idim < orig->GetNdimensions(); idim++) {
    result->GetAxis(idim)->SetTitle(orig->GetAxis(idim)->GetTitle());
    result->GetAxis(idim)->Set(nbins[idim], bins[idim]);
    if (!CheckExactRebin(orig->GetAxis(idim), result->GetAxis(idim))) {
      Printf("WARNING: unable to exact rebin axis %s of histogram %s", orig->GetAxis(idim)->GetTitle(), name);
    }
  }

  Int_t* coord_ind = new Int_t[orig->GetNdimensions()];
  Double_t* coord = new Double_t[orig->GetNdimensions()];
    
  for (Int_t ibin = 0; ibin < orig->GetNbins(); ibin++) {
    Double_t content = orig->GetBinContent(ibin, coord_ind);
    GetBinCenter(orig, coord_ind, coord);
    Int_t newbin = result->GetBin(coord);
    result->AddBinContent(newbin, content);
    result->AddBinError2(newbin, orig->GetBinError2(ibin));
  }

  delete[] coord_ind;
  delete[] coord;

  return result;
}

//____________________________________________________________________________________
TH1* DJetCorrBase::Rebin(TH1* orig, const char* name, Int_t nbins, const Double_t* bins)
{
  TH1* dest = new TH1D(name, orig->GetTitle(), nbins, bins);

  if (!CheckExactRebin(orig->GetXaxis(), dest->GetXaxis())) Printf("WARNING: unable to exact rebin axis %s of histogram %s", orig->GetXaxis()->GetTitle(), name);
  dest->GetXaxis()->SetTitle(orig->GetXaxis()->GetTitle());
  
  for (Int_t x = 0; x <= orig->GetNbinsX(); x++) {
    Double_t xlow    = orig->GetXaxis()->GetBinLowEdge(x);
    Double_t xup     = orig->GetXaxis()->GetBinUpEdge(x);
    Int_t    xlowBin = dest->GetXaxis()->FindBin(xlow);
    Int_t    xupBin  = dest->GetXaxis()->FindBin(xup);
    if (TMath::Abs(xup - dest->GetXaxis()->GetBinLowEdge(xupBin)) < fgkEpsilon) xupBin--;
    
    Double_t content = orig->GetBinContent(x);
    Double_t err     = orig->GetBinError(x);

    dest->SetBinContent(xlowBin, dest->GetBinContent(xlowBin)+content);
    dest->SetBinError(xlowBin, TMath::Sqrt(dest->GetBinError(xlowBin)*dest->GetBinError(xlowBin)+err*err));
  }

  return dest;
}

//____________________________________________________________________________________
TH2* DJetCorrBase::Rebin(TH2* orig, const char* name, Int_t nbinsx, const Double_t* binsx, Int_t nbinsy, const Double_t* binsy)
{
  TH2* dest = new TH2D(name, orig->GetTitle(), nbinsx, binsx, nbinsy, binsy);

  if (!CheckExactRebin(orig->GetXaxis(), dest->GetXaxis())) Printf("WARNING: unable to exact rebin axis %s of histogram %s", orig->GetXaxis()->GetTitle(), name);
  if (!CheckExactRebin(orig->GetYaxis(), dest->GetYaxis())) Printf("WARNING: unable to exact rebin axis %s of histogram %s", orig->GetYaxis()->GetTitle(), name);

  dest->GetXaxis()->SetTitle(orig->GetXaxis()->GetTitle());
  dest->GetYaxis()->SetTitle(orig->GetYaxis()->GetTitle());
  
  for (Int_t x = 0; x <= orig->GetNbinsX(); x++) {
    Double_t xlow    = orig->GetXaxis()->GetBinLowEdge(x);
    Double_t xup     = orig->GetXaxis()->GetBinUpEdge(x);
    Int_t    xlowBin = dest->GetXaxis()->FindBin(xlow);
    Int_t    xupBin  = dest->GetXaxis()->FindBin(xup);
    if (TMath::Abs(xup - dest->GetXaxis()->GetBinLowEdge(xupBin)) < fgkEpsilon) xupBin--;

    for (Int_t y = 0; y <= orig->GetNbinsY(); y++) {
      Double_t ylow    = orig->GetYaxis()->GetBinLowEdge(y);
      Double_t yup     = orig->GetYaxis()->GetBinUpEdge(y);
      Int_t    ylowBin = dest->GetYaxis()->FindBin(ylow);
      Int_t    yupBin  = dest->GetYaxis()->FindBin(yup);
      if (TMath::Abs(yup - dest->GetYaxis()->GetBinLowEdge(yupBin)) < fgkEpsilon) yupBin--;
    
      Double_t content = orig->GetBinContent(x,y);
      Double_t err     = orig->GetBinError(x,y);

      dest->SetBinContent(xlowBin, ylowBin, dest->GetBinContent(xlowBin, ylowBin)+content);
      dest->SetBinError(xlowBin, ylowBin, TMath::Sqrt(dest->GetBinError(xlowBin, ylowBin)*dest->GetBinError(xlowBin, ylowBin)+err*err));
    }
  }

  return dest;
}

//____________________________________________________________________________________
TH2* DJetCorrBase::GetTruth(Int_t p, Bool_t copy)
{
  TH2* hist = dynamic_cast<TH2*>(GetOutputHistogram(GetTruthName(p)));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH2*>(hist->Clone(hname));
  }

  return hist;
}

//____________________________________________________________________________________
TH2* DJetCorrBase::GetMeasured(Int_t p, Bool_t copy)
{
  TH2* hist = dynamic_cast<TH2*>(GetOutputHistogram(GetMeasuredName(p)));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH2*>(hist->Clone(hname));
  }

  return hist;
}
