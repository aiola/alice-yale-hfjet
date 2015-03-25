// This class projects the THnSparse results of the AliAnalysisTaskDmesonJetCorrelations into TH1/TH2/TH3
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMap.h>
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

#include "DJetCorrAnalysis.h"

const Double_t DJetCorrAnalysis::fgkEpsilon = 1E-6;

ClassImp(DJetCorrAnalysis);

//____________________________________________________________________________________
DJetCorrAnalysis::DJetCorrAnalysis() :
  TObject(),
  fTrainName(),
  fInputPath(),
  fInputFileName(),
  fQAListName(),
  fInputDirFileName(),
  fOutputPath(),
  fOutputFileName(),
  fOverwrite(kFALSE),
  fDPtAxisTitle(),
  fDEtaAxisTitle(),
  fDPhiAxisTitle(),
  fDInvMassAxisTitle(),
  fD2prongInvMassAxisTitle(),
  fDDeltaInvMassAxisTitle(),
  fDSoftPionPtAxisTitle(),
  fDzAxisTitle(),
  fDeltaRAxisTitle(),
  fDeltaEtaAxisTitle(),
  fDeltaPhiAxisTitle(),
  fJetPtAxisTitle(),
  fJetEtaAxisTitle(),
  fJetPhiAxisTitle(),
  fJetLeadPtAxisTitle(),
  fJetAreaAxisTitle(),
  fJetConstAxisTitle(),
  fPlotFormat(),
  fSavePlots(kFALSE),
  fAnalysisParams(0),
  fTHnAxisMap(),
  fTHnSparseMapGenerated(kFALSE),
  fInputFile(0),
  fInputDirectoryFile(0),
  fInputList(0),
  fInputQAList(0),
  fDmesons(0),
  fOutputList(0)
{
  // Default ctr.

  fTHnAxisMap.SetOwnerKeyValue();
}

//____________________________________________________________________________________
DJetCorrAnalysis::DJetCorrAnalysis(const char* train, const char* path) :
  TObject(),
  fTrainName(train),
  fInputPath(path),
  fInputFileName("AnalysisResults.root"),
  fQAListName("AliAnalysisTaskSAQA_AODFilterTracks_TPC_histos"),
  fInputDirFileName("SA_DmesonJetCorr"),
  fOutputPath("../data/"),
  fOutputFileName("DJetCorr.root"),
  fOverwrite(kFALSE),
  fDPtAxisTitle("#it{p}_{T,D} (GeV/#it{c})"),
  fDEtaAxisTitle("#eta_{D}"),
  fDPhiAxisTitle("#phi_{D} (rad)"),
  fDInvMassAxisTitle("#it{M}_{D} (GeV/#it{c}^{2})"),
  fD2prongInvMassAxisTitle("#it{M}_{2-prong} (GeV/#it{c}^{2})"),
  fDDeltaInvMassAxisTitle("#it{M}_{D*} - #it{M}_{D_{0}} (GeV/#it{c}^{2})"),
  fDSoftPionPtAxisTitle("#it{p}_{T,#pi} (GeV/#it{c})"),
  fDzAxisTitle("#it{z}_{D}"),
  fDeltaRAxisTitle("#Delta R_{D-jet}"),
  fDeltaEtaAxisTitle("#eta_{D} - #eta_{jet}"),
  fDeltaPhiAxisTitle("#phi_{D} - #phi_{jet} (rad)"),
  fJetPtAxisTitle("#it{p}_{T,jet} (GeV/#it{c})"),
  fJetEtaAxisTitle("#eta_{jet}"),
  fJetPhiAxisTitle("#phi_{jet} (rad)"),
  fJetLeadPtAxisTitle("#it{p}_{T,particle}^{leading} (GeV/#it{c})"),
  fJetAreaAxisTitle("#it{A}_{jet}"),
  fJetConstAxisTitle("No. of constituents"),
  fPlotFormat("pdf"),
  fSavePlots(kFALSE),
  fAnalysisParams(new TList()),
  fTHnAxisMap(),
  fTHnSparseMapGenerated(kFALSE),
  fInputFile(0),
  fInputDirectoryFile(0),
  fInputList(0),
  fInputQAList(0),
  fDmesons(0),
  fOutputList(0)
{
  // Standard ctr.

  fTHnAxisMap.SetOwnerKeyValue();
}

//____________________________________________________________________________________
void DJetCorrAnalysis::AddAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius)
{
  // Add the analysis params for the D meson type, jet type and jet radius provided.

  if (!fAnalysisParams) fAnalysisParams = new TList();

  fAnalysisParams->Add(new DJetCorrAnalysisParams(dmeson, jetType, jetRadius));
}

//____________________________________________________________________________________
void DJetCorrAnalysis::AddAnalysisParams(DJetCorrAnalysisParams* params)
{
  // Add the analysis params.

  if (!fAnalysisParams) fAnalysisParams = new TList();

  fAnalysisParams->Add(new DJetCorrAnalysisParams(*params));
}

//____________________________________________________________________________________
void DJetCorrAnalysis::GenerateAxisMap(THnSparse* hn)
{
  fTHnAxisMap.Delete();

  if (!LoadTHnSparse()) {
    return;
  }
  
  for (Int_t i = 0; i < hn->GetNdimensions(); i++) {
    TObjString* key = new TObjString(hn->GetAxis(i)->GetTitle());
    TParameter<Int_t>* value = new TParameter<Int_t>(hn->GetAxis(i)->GetTitle(), i);
    fTHnAxisMap.Add(key, value);
  }

  fTHnSparseMapGenerated = kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::Init()
{
  // Init class.

  Printf("Info-DJetCorrAnalysis::Init : Now initializing.");
  
  if (fOutputList) {
    delete fOutputList;
    fOutputList = 0;
  }
  fOutputList = new TList();

  Printf("Info-DJetCorrAnalysis::Init : Initialization done.");

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrAnalysis::ClearInputData()
{
  // Clear the input data.
  
  fDmesons = 0;
  fTHnAxisMap.Delete();
  fTHnSparseMapGenerated = kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateQAHistograms()
{
  // Generate QA histograms.

  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;
  }
  
  result = LoadQAList();
  if (!result) return kFALSE;
  
  result = ProjectQA();
  if (!result) return kFALSE;

  TH1::AddDirectory(addDirStatus);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateDJetCorrHistograms()
{
  // Generate D-jet correlation histograms.

  DJetCorrAnalysisParams* params = 0;
  TIter next(fAnalysisParams);

  while ((params = static_cast<DJetCorrAnalysisParams*>(next()))) GenerateDJetCorrHistograms(params);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateDJetCorrHistograms(Int_t i)
{
  // Generate D-jet correlation histograms.

  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(i));

  if (params) return GenerateDJetCorrHistograms(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateDJetCorrHistograms(const char* paramName)
{
  // Generate D-jet correlation histograms.

  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->FindObject(paramName));

  if (params) return GenerateDJetCorrHistograms(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateDJetCorrHistograms(DJetCorrAnalysisParams* params)
{
  // Generate D-jet correlation histograms.
  
  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;
  }
  
  result = LoadInputList(params->GetInputListName());
  if (!result) return kFALSE;

  result = ProjectCorrD(params);
  if (!result) return kFALSE;

  TH1::AddDirectory(addDirStatus);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDJetCorrHistograms()
{
  DJetCorrAnalysisParams* params = 0;
  TIter next(fAnalysisParams);

  while ((params = static_cast<DJetCorrAnalysisParams*>(next()))) PlotDJetCorrHistograms(params);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDJetCorrHistograms(Int_t i)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(i));

  if (params) return PlotDJetCorrHistograms(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDJetCorrHistograms(const char* paramName)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->FindObject(paramName));

  if (params) return PlotDJetCorrHistograms(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDJetCorrHistograms(DJetCorrAnalysisParams* params)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;

    result = LoadOutputHistograms();
    if (!result) return kFALSE;
  }
  
  // Invariant mass
  Printf("Info-DJetCorrAnalysis::PlotDJetCorrHistograms : Plotting invariant mass vs D pt for all D mesons.");
  PlotInvMassHistogramsVsDPt(params, 0, params->GetJetPtBin(params->GetNJetPtBins()));
  Printf("Info-DJetCorrAnalysis::PlotDJetCorrHistograms : Plotting invariant mass vs D pt for uncorrelated D mesons.");
  PlotInvMassHistogramsVsDPt(params, 0, 1);
  Printf("Info-DJetCorrAnalysis::PlotDJetCorrHistograms : Plotting invariant mass vs D pt for correlated D mesons.");
  PlotInvMassHistogramsVsDPt(params, params->GetJetPtBin(0), params->GetJetPtBin(params->GetNJetPtBins()));
  Printf("Info-DJetCorrAnalysis::PlotDJetCorrHistograms : Plotting invariant mass vs z for correlated D mesons.");
  PlotInvMassHistogramsVsDz(params);

  // Delta R
  PlotDeltaRVsDPt(params, params->GetJetPtBin(0), params->GetJetPtBin(params->GetNJetPtBins()));
  
  TH1::AddDirectory(addDirStatus);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::LoadOutputHistograms()
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
    Printf("Error-DJetCorrAnalysis::LoadOutputHistograms : Could not open file '%s' to read.", fname.Data()); 
    outputFile = 0;
    return kFALSE;
  }
  
  TList *keys = outputFile->GetListOfKeys();
  if (!keys) {
    Printf("Error-DJetCorrAnalysis::LoadOutputHistograms : Could not get keys from file '%s'.", fname.Data()); 
    return kFALSE;
  }

  TIter next(keys);
  TObject* key = 0;
  while ((key = next())) {
    TString objname(key->GetName());
    TObject* obj = outputFile->Get(objname);
    fOutputList->Add(obj);
  }

  outputFile->Close();
  delete outputFile;
  outputFile = 0;

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrAnalysis::SavePlot(TCanvas* canvas)
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

  fname += canvas->GetTitle();
  fname += ".";
  fname += fPlotFormat;

  canvas->SaveAs(fname);
}

//____________________________________________________________________________________
TVirtualPad* DJetCorrAnalysis::SetUpPad(TVirtualPad* pad,
                                        const char* xTitle, Double_t minX, Double_t maxX, Bool_t logX,
                                        const char* yTitle, Double_t minY, Double_t maxY, Bool_t logY)
{
  if (!pad) return 0;

  pad->SetLogx(logX);
  pad->SetLogy(logY);
  pad->cd();

  TString blankHistName(Form("%s_blankHist", pad->GetName()));
  TH1* blankHist = new TH1D(blankHistName, blankHistName, 1000, minX, maxX);
  blankHist->GetXaxis()->SetTitle(xTitle);
  blankHist->GetYaxis()->SetTitle(yTitle);
  blankHist->GetYaxis()->SetRangeUser(minY, maxY);
  blankHist->Draw();

  return pad;
}

//____________________________________________________________________________________
TCanvas* DJetCorrAnalysis::SetUpCanvas(const char* name,
                                       const char* xTitle, Double_t minX, Double_t maxX, Bool_t logX,
                                       const char* yTitle, Double_t minY, Double_t maxY, Bool_t logY,
                                       Double_t w, Double_t h, Int_t rows, Int_t cols)
{
  Printf("Info-DJetCorrAnalysis::SetUpCanvas : Setting up canvas '%s'", name);
  
  TCanvas* canvas = new TCanvas(name, name, w, h);
  
  if (rows == 1 && cols == 1) {
    SetUpPad(canvas, xTitle, minX, maxX, logX, yTitle, minY, maxY, logY);    
  }
  else {
    canvas->Divide(cols, rows);
    Int_t n = rows * cols;
    for (Int_t i = 1; i <= n; i++) {
      SetUpPad(canvas->cd(i), xTitle, minX, maxX, logX, yTitle, minY, maxY, logY); 
    }
  }
  
  return canvas;
}

//____________________________________________________________________________________
TLegend* DJetCorrAnalysis::SetUpLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t textSize)
{
  TLegend* leg = new TLegend(x1, y1, x2, y2);
  leg->SetTextFont(43);
  leg->SetTextSize(textSize);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  return leg;
}

//____________________________________________________________________________________
TPaveText* DJetCorrAnalysis::SetUpPaveText(Double_t x1, Double_t y1, Double_t x2, Double_t y2, Int_t textSize, const char* text)
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
Bool_t DJetCorrAnalysis::PlotTrackHistograms()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;

    result = LoadOutputHistograms();
    if (!result) return kFALSE;
  }

  TString labels[4] = {"Hybrid", "Global", "Constr", "ConstrNoITS"};
  Color_t colors[4] = {kBlack, kBlue+2, kRed+2, kGreen+2};
  TString titles[4] = {"Hybrid tracks", "Global tracks", "Constrained, w/ ITS refit", "Constrained, w/o ITS refit"};

  TCanvas* cTrackPt = 0;
  TCanvas* cTrackEta = 0;
  TCanvas* cTrackPhi = 0;

  TLegend* legPt = SetUpLegend(0.58, 0.68, 0.88, 0.88, 14);
  TLegend* legEta = SetUpLegend(0.19, 0.30, 0.49, 0.50, 14);
  TLegend* legPhi = SetUpLegend(0.58, 0.68, 0.88, 0.88, 14);
  
  for (Int_t i = 0; i < 4; i++) {
    TH1* hTrEta = static_cast<TH1*>(fOutputList->FindObject(Form("hTracks_%s_Eta", labels[i].Data())));
    TH1* hTrPhi = static_cast<TH1*>(fOutputList->FindObject(Form("hTracks_%s_Phi", labels[i].Data())));
    TH1* hTrPt = static_cast<TH1*>(fOutputList->FindObject(Form("hTracks_%s_Pt", labels[i].Data())));

    if (i == 0) {
      cTrackPt = SetUpCanvas("fig_TrackPt", "#it{p}_{T} (GeV/#it{c})", 0, 50, kFALSE, "counts", 1e-1, hTrPt->GetMaximum()*3, kTRUE);
      cTrackEta = SetUpCanvas("fig_TrackEta", "#eta", -1, 1, kFALSE, "counts", 1e-1, hTrEta->GetMaximum()*1.2, kFALSE);
      cTrackPhi = SetUpCanvas("fig_TrackPhi", "#phi", 0, TMath::TwoPi(), kFALSE, "counts", 1e-1, hTrPhi->GetMaximum()*1.5, kFALSE);  
    }
    
    cTrackPt->cd();
    TH1* hTrPt_copy = hTrPt->DrawCopy("same hist");
    hTrPt_copy->SetLineColor(colors[i]);
    hTrPt_copy->SetTitle(titles[i]);

    cTrackEta->cd();
    TH1* hTrEta_copy = hTrEta->DrawCopy("same hist");
    hTrEta_copy->SetLineColor(colors[i]);
    hTrEta_copy->SetTitle(titles[i]);

    cTrackPhi->cd();
    TH1* hTrPhi_copy = hTrPhi->DrawCopy("same hist");
    hTrPhi_copy->SetLineColor(colors[i]);
    hTrPhi_copy->SetTitle(titles[i]);

    legPt->AddEntry(hTrPt_copy, titles[i], "l");
    legEta->AddEntry(hTrEta_copy, titles[i], "l");
    legPhi->AddEntry(hTrPhi_copy, titles[i], "l");
  }

  cTrackPt->cd();
  legPt->Draw();
  cTrackPt->Update();

  cTrackEta->cd();
  legEta->Draw();
  cTrackEta->Update();

  cTrackPhi->cd();
  legPhi->Draw();
  cTrackPhi->Update();

  if (fSavePlots) {
    SavePlot(cTrackPt);
    SavePlot(cTrackEta);
    SavePlot(cTrackPhi);
  }

  TH1::AddDirectory(addDirStatus);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDeltaRVsDPt(DJetCorrAnalysisParams* params, Double_t minJetPt, Double_t maxJetPt)
{
  if (!fOutputList) return kFALSE;

  TString jetCuts;
  
  TString cname("fig_DeltaRVsDPt_");
  cname += params->GetName();
  
  if (maxJetPt > 1 && minJetPt > 0) {
    jetCuts = Form("%.1f < #it{p}_{T,jet} < %.1f GeV/#it{c}", minJetPt, maxJetPt);
    cname += Form("_JetPt_%03.0f_%03.0f", minJetPt, maxJetPt);
  }
  else {
    return kFALSE;
  }

  TString hname("DeltaR");
  TString xTitle("#DeltaR = #sqrt{#Delta#eta^{2} + #Delta#phi^{2}}");
  TString yTitle("Probability density");
  TCanvas *canvas = SetUpCanvas(cname, xTitle, 0, 0.5, kFALSE, yTitle, 0, 1, kFALSE);
  TH1* hblank = dynamic_cast<TH1*>(canvas->GetListOfPrimitives()->At(0));
  Color_t colors[12] = {kRed+1, kBlue+1, kMagenta, kGreen+2, kOrange+1, kCyan+2, kPink+1, kTeal+1, kViolet, kAzure+1, kYellow+2, kSpring+3};

  TString prefix(params->GetName());
  TLegend* leg = SetUpLegend(0.32, 0.49, 0.85, 0.81, 13);
  leg->SetNColumns(2);
  Double_t maxY = 0;
  for (Int_t i = 0; i < params->GetNDPtBins()-1; i+=2) {
    TString cuts(Form("JetPt_%03.0f_%03.0f_DPt_%02.0f_%02.0f_z_%.1f_%.1f", minJetPt, maxJetPt, params->GetDPtBin(i), params->GetDPtBin(i+1), params->GetzBin(0), params->GetzBin(params->GetNzBins())));
    cuts.ReplaceAll(".", "");
    TString objname(Form("h%s_%s_%s", prefix.Data(), hname.Data(), cuts.Data()));
    //Printf("Info-DJetCorrAnalysis::PlotDeltaRVsDPt : Retrieving histogram '%s'", objname.Data());
    TH1* hist = static_cast<TH1*>(fOutputList->FindObject(objname));
    if (!hist) {
      Printf("Error-DJetCorrAnalysis::PlotDeltaRVsDPt : Histogram '%s' not found!", objname.Data());
      continue;
    }
    TString newName(objname);
    newName += "_copy";
    //Printf("Info-DJetCorrAnalysis::PlotDeltaRVsDPt : Cloning histogram '%s'", objname.Data());
    TH1* hcopy = static_cast<TH1*>(hist->Clone(newName));
    hcopy->Rebin(4);
    hcopy->Scale(1. / hcopy->Integral(), "width");
    if (maxY < hcopy->GetMaximum()) maxY = hcopy->GetMaximum();
    hcopy->SetLineColor(colors[i/2]);
    hcopy->SetMarkerColor(colors[i/2]);
    hcopy->SetMarkerStyle(kFullCircle);
    hcopy->SetMarkerSize(0.5);
    hcopy->Draw("same");
    //Printf("Info-DJetCorrAnalysis::PlotDeltaRVsDPt : Setting title of histogram '%s'", histos[i]->GetName());
    TString htitle(Form("%.1f < #it{p}_{T,D} < %.1f GeV/#it{c}", params->GetDPtBin(i), params->GetDPtBin(i+1)));
    TLegendEntry* legEntry = leg->AddEntry(hcopy, htitle, "pe");
    legEntry->SetLineColor(colors[i/2]);
  }

  if (hblank && maxY > 0) {
    hblank->GetYaxis()->SetRangeUser(0, maxY*1.4);
  }

  leg->Draw();

  if (fSavePlots) {
    SavePlot(canvas);
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotInvMassHistogramsVsDPt(DJetCorrAnalysisParams* params, Double_t minJetPt, Double_t maxJetPt)
{
  if (!fOutputList) return kFALSE;
  
  Double_t D0mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(421))->Mass();
  Double_t Dstarmass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(413))->Mass();

  TString jetCuts;
  
  TString cname("fig_InvMassVsDPt");
  cname += params->GetName();
  
  if (maxJetPt > 1 && minJetPt > 0) {
    jetCuts = Form("%.1f < #it{p}_{T,jet} < %.1f GeV/#it{c}", minJetPt, maxJetPt);
    cname += Form("_JetPt_%03.0f_%03.0f", minJetPt, maxJetPt);
  }
  else if (maxJetPt > 1) {
    jetCuts = "";
  }
  else {
    jetCuts = "No correlated jet";
    cname += "_NoCorrJet";
  }
  
  TString hname;
  TString xTitle;
  Double_t minMass = 0;
  Double_t maxMass = 0;
  Double_t pdgMass = -1;
  if (params->IsD0()) {
    xTitle = "#it{m}(K#pi) (GeV/#it{c}^{2})";
    hname = "InvMass";
    minMass = D0mass - 0.1;
    maxMass = D0mass + 0.1;
    pdgMass = D0mass;
  }
  else if (params->IsDStar()) {
    xTitle = "#it{m}(K#pi#pi) - #it{m}(K#pi) (GeV/#it{c}^{2})";
    hname = "DeltaInvMass";
    pdgMass = Dstarmass - D0mass;
    minMass = pdgMass - 0.04;
    maxMass = pdgMass + 0.04;
  }
  else {
    Printf("Error-DJetCorrAnalysis::PlotInvMassHistogramsVsDPt : Meson type '%s' not recognized!", params->GetDmesonName());
    return kFALSE;
  }
  
  TString prefix(params->GetName());
  
  TH1** histos = new TH1*[params->GetNDPtBins()];
  Int_t n = 0;
  for (Int_t i = 0; i < params->GetNDPtBins(); i++) {
    TString cuts(Form("JetPt_%03.0f_%03.0f_DPt_%02.0f_%02.0f_z_%.1f_%.1f", minJetPt, maxJetPt, params->GetDPtBin(i), params->GetDPtBin(i+1), params->GetzBin(0), params->GetzBin(params->GetNzBins())));
    cuts.ReplaceAll(".", "");
    TString objname(Form("h%s_%s_%s", prefix.Data(), hname.Data(), cuts.Data()));
    //Printf("Info-DJetCorrAnalysis::PlotInvMassHistogramsVsDPt : Retrieving histogram '%s'", objname.Data());
    TH1* hist = static_cast<TH1*>(fOutputList->FindObject(objname));
    if (!hist) {
      Printf("Error-DJetCorrAnalysis::PlotInvMassHistogramsVsDPt : Histogram '%s' not found!", objname.Data());
      continue;
    }
    TString newName(objname);
    newName += "_copy";
    //Printf("Info-DJetCorrAnalysis::PlotInvMassHistogramsVsDPt : Cloning histogram '%s'", objname.Data());
    histos[n] = static_cast<TH1*>(hist->Clone(newName));
    //Printf("Info-DJetCorrAnalysis::PlotInvMassHistogramsVsDPt : Setting title of histogram '%s'", histos[i]->GetName());
    TString htitle(Form("%.1f < #it{p}_{T,D} < %.1f GeV/#it{c}\n%s", params->GetDPtBin(i), params->GetDPtBin(i+1), jetCuts.Data()));
    histos[n]->SetTitle(htitle);
    n++;
  }
  
  Bool_t result = PlotInvMassHistogramArray(n, histos, cname, xTitle, minMass, maxMass, pdgMass);

  delete[] histos;
  
  return result;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotInvMassHistogramsVsDz(DJetCorrAnalysisParams* params)
{
  if (!fOutputList) return kFALSE;
  
  Double_t D0mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(421))->Mass();
  Double_t Dstarmass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(413))->Mass();

  TString jetCuts;
  
  TString cname("fig_InvMassVsDz");
  cname += params->GetName();
  
  jetCuts = Form("%.1f < #it{p}_{T,jet} < %.1f GeV/#it{c}", params->GetJetPtBin(0), params->GetJetPtBin(params->GetNJetPtBins()));
  cname += Form("_JetPt_%03.0f_%03.0f", params->GetJetPtBin(0), params->GetJetPtBin(params->GetNJetPtBins()));
  
  TString hname;
  TString xTitle;
  Double_t minMass = 0;
  Double_t maxMass = 0;
  Double_t pdgMass = -1;
  if (params->IsD0()) {
    xTitle = "#it{m}(K#pi) (GeV/#it{c}^{2})";
    hname = "InvMass";
    pdgMass = D0mass;
    minMass = D0mass - 0.04;
    maxMass = D0mass + 0.04;
  }
  else if (params->IsDStar()) {
    xTitle = "#it{m}(K#pi#pi) - #it{m}(K#pi) (GeV/#it{c}^{2})";
    hname = "DeltaInvMass";
    pdgMass = Dstarmass - D0mass;
    minMass = pdgMass - 0.04;
    maxMass = pdgMass + 0.04;
  }
  else {
    Printf("Error-DJetCorrAnalysis::PlotInvMassHistogramsVsDz : Meson type '%s' not recognized!", params->GetDmesonName());
    return kFALSE;
  }
  
  TString prefix(params->GetName());
  
  TH1** histos = new TH1*[params->GetNzBins()];
  Int_t n = 0;
  for (Int_t i = 0; i < params->GetNzBins(); i++) {
    TString cuts(Form("JetPt_%03.0f_%03.0f_DPt_%02.0f_%02.0f_z_%.1f_%.1f", params->GetJetPtBin(0), params->GetJetPtBin(params->GetNJetPtBins()), params->GetDPtBin(0), params->GetDPtBin(params->GetNDPtBins()), params->GetzBin(i), params->GetzBin(i+1)));
    cuts.ReplaceAll(".", "");
    TString objname(Form("h%s_%s_%s", prefix.Data(), hname.Data(), cuts.Data()));
    //Printf("Info-DJetCorrAnalysis::PlotInvMassHistogramsVsDz : Retrieving histogram '%s'", objname.Data());
    TH1* hist = static_cast<TH1*>(fOutputList->FindObject(objname));
    if (!hist) {
      Printf("Error-DJetCorrAnalysis::PlotInvMassHistogramsVsDz : Histogram '%s' not found!", objname.Data());
      continue;
    }
    TString newName(objname);
    newName += "_copy";
    //Printf("Info-DJetCorrAnalysis::PlotInvMassHistogramsVsDz : Cloning histogram '%s'", objname.Data());
    histos[n] = static_cast<TH1*>(hist->Clone(newName));
    //Printf("Info-DJetCorrAnalysis::PlotInvMassHistogramsVsDz : Setting title of histogram '%s'", histos[i]->GetName());
    TString htitle(Form("%.1f < #it{z}_{D} < %.1f \n%s", params->GetzBin(i), params->GetzBin(i+1), jetCuts.Data()));
    histos[n]->SetTitle(htitle);
    n++;
  }
  
  Bool_t result = PlotInvMassHistogramArray(n, histos, cname, xTitle, minMass, maxMass, pdgMass);

  delete[] histos;

  return result;
}    

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotInvMassHistogramArray(Int_t n, TH1** histos, const char* name, const char* xTitle, Double_t minMass, Double_t maxMass, Double_t pdgMass)
{
  // Plot invariant mass histograms contained in histos.

  if (n == 0) return kFALSE;
  
  Printf("Info-DJetCorrAnalysis::PlotInvMassHistograms : Plotting invariant mass histograms '%s'", name);
  
  Int_t cols = TMath::FloorNint(TMath::Sqrt(n));
  Int_t rows = TMath::CeilNint(1. * n / cols);

  TString cname(name);

  Double_t w = cols*250;
  Double_t h = rows*250;
  
  TCanvas* canvas = SetUpCanvas(cname, xTitle, minMass, maxMass, kFALSE, "counts", 0, 1, kFALSE, h, w, cols, rows);
  for (Int_t i = 0; i < n; i++) {
    TVirtualPad* pad = canvas->cd(i+1);
    TH1* blank = dynamic_cast<TH1*>(pad->GetListOfPrimitives()->At(0));
    if (!blank) {
      Printf("Error-DJetCorrAnalysis::PlotInvMassHistograms : Could not find blank histogram!");
      continue;
    }
    blank->GetYaxis()->SetRangeUser(0, histos[i]->GetMaximum()*1.4);

    histos[i]->Sumw2();
    Printf("Info-DJetCorrAnalysis::PlotInvMassHistograms : Now plotting '%s'", histos[i]->GetName());
    histos[i]->DrawCopy("same p");

    TPaveText* pave = SetUpPaveText(0.15, 0.68, 0.90, 0.86, 13, histos[i]->GetTitle());
    pave->Draw();

    if (pdgMass > 0) {
      TLine *line = new TLine(pdgMass, 0, pdgMass, histos[i]->GetMaximum());
      line->SetLineColor(kRed);
      line->SetLineWidth(1);
      line->Draw();
    }
  }

  if (fSavePlots) SavePlot(canvas);

  return kTRUE;
}

//____________________________________________________________________________________
Int_t DJetCorrAnalysis::GetAxisIndex(TString title)
{
  if (!fTHnSparseMapGenerated) {
    GenerateAxisMap(fDmesons);
    if (!fTHnSparseMapGenerated) return -1;
  }

  TParameter<Int_t>* par = static_cast<TParameter<Int_t>*>(fTHnAxisMap.GetValue(title));
  if (!par) {
    //Printf("Warning-DJetCorrAnalysis::GetAxisIndex : Could not find axis with title '%s'", title.Data());
    return -1;
  }
  
  return par->GetVal();
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectQA()
{
  // Project QA histograms.

  TH1* evCount = dynamic_cast<TH1*>(fInputQAList->FindObject("fHistEventCount"));
  evCount->SetName("hEventCount");
  fOutputList->Add(evCount);

  TString labels[4] = {"Global", "Constr", "ConstrNoITS", "Hybrid"};
  
  for (Int_t i = 0; i <= 3; i++) {
    TH3* hTrPhiEtaPt = dynamic_cast<TH3*>(fInputQAList->FindObject(Form("fHistTrPhiEtaPt_0_%d",i)));
    if (!hTrPhiEtaPt) {
      Printf("Error-DJetCorrAnalysis::ProjectQA : Could not find track histogram for track type %d!", i);
      continue;
    }
    TH1* hTrEta = hTrPhiEtaPt->ProjectionX(Form("hTracks_%s_Eta", labels[i].Data()));
    TH1* hTrPhi = hTrPhiEtaPt->ProjectionY(Form("hTracks_%s_Phi", labels[i].Data()));
    TH1* hTrPt = hTrPhiEtaPt->ProjectionZ(Form("hTracks_%s_Pt", labels[i].Data()));
    
    fOutputList->Add(hTrEta);
    fOutputList->Add(hTrPhi);
    fOutputList->Add(hTrPt);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectCorrD(DJetCorrAnalysisParams* params)
{
  // Project histograms related to the D meson correlated to a jet.

  TString prefix(params->GetName());
  TString suffix("");

  ProjectDJetCorr(prefix, suffix, kFALSE, kTRUE, kTRUE,
                  0, params->GetJetPtBin(params->GetNJetPtBins()),
                  params->GetDPtBin(0), params->GetDPtBin(params->GetNDPtBins()), params->GetzBin(0), params->GetzBin(params->GetNzBins()), -0.9, 0.9,
                  params->GetInvMinMass(), params->GetInvMaxMass(), params->Get2ProngMinMass(), params->Get2ProngMaxMass(), params->GetDeltaInvMinMass(), params->GetDeltaInvMaxMass());

  ProjectDJetCorr(prefix, suffix, kFALSE, kTRUE, kTRUE, 
                  0, 1,
                  params->GetDPtBin(0), params->GetDPtBin(params->GetNDPtBins()), params->GetzBin(0), params->GetzBin(params->GetNzBins()), -0.9, 0.9,
                  params->GetInvMinMass(), params->GetInvMaxMass(), params->Get2ProngMinMass(), params->Get2ProngMaxMass(), params->GetDeltaInvMinMass(), params->GetDeltaInvMaxMass());

  ProjectDJetCorr(prefix, suffix, kTRUE, kTRUE, kTRUE, 
                  params->GetJetPtBin(0), params->GetJetPtBin(params->GetNJetPtBins()),
                  params->GetDPtBin(0), params->GetDPtBin(params->GetNDPtBins()), params->GetzBin(0), params->GetzBin(params->GetNzBins()), -0.9, 0.9,
                  params->GetInvMinMass(), params->GetInvMaxMass(), params->Get2ProngMinMass(), params->Get2ProngMaxMass(), params->GetDeltaInvMinMass(), params->GetDeltaInvMaxMass());

  for (Int_t i = 0; i < params->GetNDPtBins(); i++) {
    ProjectDJetCorr(prefix, suffix, kFALSE, kFALSE, kTRUE, 
                    0, params->GetJetPtBin(params->GetNJetPtBins()),
                    params->GetDPtBin(i), params->GetDPtBin(i+1), params->GetzBin(0), params->GetzBin(params->GetNzBins()), -0.9, 0.9,
                    params->GetInvMinMass(), params->GetInvMaxMass(), params->Get2ProngMinMass(), params->Get2ProngMaxMass(), params->GetDeltaInvMinMass(), params->GetDeltaInvMaxMass());
  }

  for (Int_t i = 0; i < params->GetNDPtBins(); i++) {
    ProjectDJetCorr(prefix, suffix, kFALSE, kFALSE, kTRUE, 
                    0, 1,
                    params->GetDPtBin(i), params->GetDPtBin(i+1), params->GetzBin(0), params->GetzBin(params->GetNzBins()), -0.9, 0.9,
                    params->GetInvMinMass(), params->GetInvMaxMass(), params->Get2ProngMinMass(), params->Get2ProngMaxMass(), params->GetDeltaInvMinMass(), params->GetDeltaInvMaxMass());
  }
  
  for (Int_t i = 0; i < params->GetNDPtBins(); i++) {
    ProjectDJetCorr(prefix, suffix, kTRUE, kFALSE, kTRUE, 
                    params->GetJetPtBin(0), params->GetJetPtBin(params->GetNJetPtBins()),
                    params->GetDPtBin(i), params->GetDPtBin(i+1), params->GetzBin(0), params->GetzBin(params->GetNzBins()), -0.9, 0.9,
                    params->GetInvMinMass(), params->GetInvMaxMass(), params->Get2ProngMinMass(), params->Get2ProngMaxMass(), params->GetDeltaInvMinMass(), params->GetDeltaInvMaxMass());
  }

  for (Int_t i = 0; i < params->GetNJetPtBins(); i++) {
    ProjectDJetCorr(prefix, suffix, kTRUE, kTRUE, kTRUE, 
                    params->GetJetPtBin(i), params->GetJetPtBin(i+1),
                    params->GetDPtBin(0), params->GetDPtBin(params->GetNDPtBins()), params->GetzBin(0), params->GetzBin(params->GetNzBins()), -0.9, 0.9,
                    params->GetInvMinMass(), params->GetInvMaxMass(), params->Get2ProngMinMass(), params->Get2ProngMaxMass(), params->GetDeltaInvMinMass(), params->GetDeltaInvMaxMass());
  }

  for (Int_t i = 0; i < params->GetNzBins(); i++) {
    ProjectDJetCorr(prefix, suffix, kTRUE, kTRUE, kFALSE, 
                    params->GetJetPtBin(0), params->GetJetPtBin(params->GetNJetPtBins()),
                    params->GetDPtBin(0), params->GetDPtBin(params->GetNDPtBins()), params->GetzBin(i), params->GetzBin(i+1), -0.9, 0.9,
                    params->GetInvMinMass(), params->GetInvMaxMass(), params->Get2ProngMinMass(), params->Get2ProngMaxMass(), params->GetDeltaInvMinMass(), params->GetDeltaInvMaxMass());
  }

  TString cutsN(Form("JetPt_%03.0f_%03.0f", params->GetJetPtBin(0), params->GetJetPtBin(params->GetNJetPtBins())));
  cutsN.ReplaceAll(".", "");

  TString cutsD(Form("JetPt_%03.0f_%03.0f", 0., params->GetJetPtBin(params->GetNJetPtBins())));
  cutsD.ReplaceAll(".", "");

  GenerateRatios(cutsN, cutsD);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateRatios(const char* nname, const char* dname)
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

    Printf("Info-DJetCorrAnalysis::GenerateRatios : Now calcutaling ratio '%s' over '%s'", hnname.Data(), hdname.Data());

    hist1D->Sumw2();
    num->Sumw2();

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
Bool_t DJetCorrAnalysis::ProjectDJetCorr(TString prefix, TString suffix, Bool_t doCorrPlots, Bool_t doPtPlots, Bool_t dozPlots,
                                         Double_t minJetPt, Double_t maxJetPt,
                                         Double_t minDPt, Double_t maxDPt, Double_t minz, Double_t maxz, Double_t minDEta, Double_t maxDEta,
                                         Double_t minInvMass, Double_t maxInvMass, Double_t min2ProngMass, Double_t max2ProngMass, Double_t minDeltaInvMass, Double_t maxDeltaInvMass)
{
  // Project histograms related to the D meson with specified cuts.

  if (!LoadTHnSparse()) {
    return kFALSE;
  }

  TString jetCuts(Form("JetPt_%03.0f_%03.0f", minJetPt, maxJetPt));
  jetCuts.ReplaceAll(".", "");
  
  TString dCuts(Form("DPt_%02.0f_%02.0f", minDPt, maxDPt));
  dCuts.ReplaceAll(".", "");

  TString zCuts(Form("z_%.1f_%.1f", minz, maxz));
  zCuts.ReplaceAll(".", "");
    
  TString cuts(Form("%s_%s_%s", jetCuts.Data(), dCuts.Data(), zCuts.Data()));

  if (!suffix.IsNull()) suffix.Prepend("_");

  Int_t dPtAxis = GetAxisIndex(fDPtAxisTitle);
  if (dPtAxis < 0) return kFALSE;
  
  Int_t dEtaAxis = GetAxisIndex(fDEtaAxisTitle);
  if (dEtaAxis < 0) return kFALSE;
  
  Int_t dPhiAxis = GetAxisIndex(fDPhiAxisTitle);
  if (dPhiAxis < 0) return kFALSE;
  
  Int_t jetPtAxis = GetAxisIndex(fJetPtAxisTitle);
  if (jetPtAxis < 0) return kFALSE;
    
  Int_t jetEtaAxis = GetAxisIndex(fJetEtaAxisTitle);
  if (jetEtaAxis < 0) return kFALSE;
  
  Int_t jetPhiAxis = GetAxisIndex(fJetPhiAxisTitle);
  if (jetPhiAxis < 0) return kFALSE;

  Int_t dInvMassAxis = GetAxisIndex(fDInvMassAxisTitle);
  Int_t d2ProngInvMassAxis = GetAxisIndex(fD2prongInvMassAxisTitle);
  Int_t dDeltaInvMassAxis = GetAxisIndex(fDDeltaInvMassAxisTitle);
  Int_t dSoftPionPtAxis = GetAxisIndex(fDSoftPionPtAxisTitle);
  
  Int_t deltaRAxis = GetAxisIndex(fDeltaRAxisTitle);
  Int_t deltaEtaAxis = GetAxisIndex(fDeltaEtaAxisTitle);
  Int_t deltaPhiAxis = GetAxisIndex(fDeltaPhiAxisTitle);
  Int_t dzAxis = GetAxisIndex(fDzAxisTitle);

  fDmesons->GetAxis(jetPtAxis)->SetRangeUser(minJetPt * (1+fgkEpsilon), maxJetPt * (1-fgkEpsilon));
  fDmesons->GetAxis(dPtAxis)->SetRangeUser(minDPt * (1+fgkEpsilon), maxDPt * (1-fgkEpsilon));
  fDmesons->GetAxis(dEtaAxis)->SetRangeUser(minDEta * (1+fgkEpsilon), maxDEta * (1-fgkEpsilon));

  if (0 && dInvMassAxis >= 0) {
    fDmesons->GetAxis(dInvMassAxis)->SetRangeUser(minInvMass * (1+fgkEpsilon), maxInvMass * (1-fgkEpsilon));    
  }

  if (0 && d2ProngInvMassAxis >= 0) {
    fDmesons->GetAxis(d2ProngInvMassAxis)->SetRangeUser(min2ProngMass * (1+fgkEpsilon), max2ProngMass * (1-fgkEpsilon));    
  }

  if (dDeltaInvMassAxis >= 0) {
    fDmesons->GetAxis(dDeltaInvMassAxis)->SetRangeUser(minDeltaInvMass * (1+fgkEpsilon), maxDeltaInvMass * (1-fgkEpsilon));    
  }
  
  if (dzAxis >= 0) {
    fDmesons->GetAxis(dzAxis)->SetRangeUser(minz * (1+fgkEpsilon), maxz * (1-fgkEpsilon));    
  }

  TString hname;
  
  if (doPtPlots) {
    hname = Form("h%s_MesonPt_%s_%s%s", prefix.Data(), jetCuts.Data(), zCuts.Data(), suffix.Data());
    if (!fOutputList->Contains(hname)) {
      TH1* hDpt = fDmesons->Projection(dPtAxis, "EO");
      hDpt->SetName(hname);
      fOutputList->Add(hDpt);

      if (dSoftPionPtAxis >= 0) {
        TH2* hdsoftpionptVsDpt = fDmesons->Projection(dSoftPionPtAxis, dPtAxis, "EO");
        hdsoftpionptVsDpt->SetName(Form("h%s_SoftPionVsDPt_%s_%s%s", prefix.Data(), jetCuts.Data(), zCuts.Data(), suffix.Data()));
        fOutputList->Add(hdsoftpionptVsDpt);
      }
    }
  }
  
  TH2* hDpos = fDmesons->Projection(dPhiAxis, dEtaAxis, "EO");
  hDpos->SetName(Form("h%s_MesonPhiVsEta_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
  fOutputList->Add(hDpos);

  if (dInvMassAxis >= 0) {
    TH1* hdinvmass = fDmesons->Projection(dInvMassAxis, "EO");
    //hdinvmass->Rebin(5);
    hdinvmass->SetName(Form("h%s_InvMass_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
    Printf("Info-DJetCorrAnalysis::ProjectDJetCorr : Adding histogram '%s'", hdinvmass->GetName());
    fOutputList->Add(hdinvmass);
  }

  if (d2ProngInvMassAxis >= 0) {
    TH1* hd2pronginvmass = fDmesons->Projection(d2ProngInvMassAxis, "EO");
    //hd2pronginvmass->Rebin(5);
    hd2pronginvmass->SetName(Form("h%s_D0InvMass_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
    Printf("Info-DJetCorrAnalysis::ProjectDJetCorr : Adding histogram '%s'", hd2pronginvmass->GetName());
    fOutputList->Add(hd2pronginvmass);
  }

  if (dDeltaInvMassAxis >= 0) {
    TH1* hddeltainvmass = fDmesons->Projection(dDeltaInvMassAxis, "EO");
    //hddeltainvmass->Rebin(5);
    hddeltainvmass->SetName(Form("h%s_DeltaInvMass_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
    Printf("Info-DJetCorrAnalysis::ProjectDJetCorr : Adding histogram '%s'", hddeltainvmass->GetName());
    fOutputList->Add(hddeltainvmass);

    if (dSoftPionPtAxis >= 0) {
      TH2* hdsoftpionptVsDeltaInvMass = fDmesons->Projection(dSoftPionPtAxis, dDeltaInvMassAxis, "EO");
      hdsoftpionptVsDeltaInvMass->SetName(Form("h%s_SoftPionVsDeltaInvMass_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hdsoftpionptVsDeltaInvMass);
    }
  }

  if (dSoftPionPtAxis >= 0) {
    TH1* hdsoftpionpt = fDmesons->Projection(dSoftPionPtAxis, "EO");
    hdsoftpionpt->SetName(Form("h%s_SoftPion_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
    fOutputList->Add(hdsoftpionpt);
  }

  if (doCorrPlots) {    
    if (dzAxis >= 0 && dozPlots) {
      TH1* hdz = fDmesons->Projection(dzAxis, "EO");
      hdz->SetName(Form("h%s_MesonZ_%s_%s%s", prefix.Data(), jetCuts.Data(), dCuts.Data(), suffix.Data()));
      fOutputList->Add(hdz);

      TH2* hdzVsDPt = fDmesons->Projection(dzAxis, dPtAxis, "EO");
      hdzVsDPt->SetName(Form("h%s_MesonZvsDPt_%s_%s%s", prefix.Data(), jetCuts.Data(), dCuts.Data(), suffix.Data()));
      fOutputList->Add(hdzVsDPt);

      if (jetPtAxis >= 0) {
        TH2* hdzVsJetPt = fDmesons->Projection(dzAxis, jetPtAxis, "EO");
        hdzVsJetPt->SetName(Form("h%s_MesonZvsJetPt_%s_%s%s", prefix.Data(), jetCuts.Data(), dCuts.Data(), suffix.Data()));
        fOutputList->Add(hdzVsJetPt);
      }
      
      if (deltaRAxis >= 0) {
        TH2* hdzVsDeltaR = fDmesons->Projection(dzAxis, deltaRAxis, "EO");
        hdzVsDeltaR->SetName(Form("h%s_MesonZvsDeltaR_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hdzVsDeltaR);
      }
    }
        
    if (deltaRAxis >= 0) {
      TH1* hDeltaR = fDmesons->Projection(deltaRAxis, "EO");
      hDeltaR->SetName(Form("h%s_DeltaR_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDeltaR);

      TH1* hDeltaRVsDpt = fDmesons->Projection(deltaRAxis, dPtAxis, "EO");
      hDeltaRVsDpt->SetName(Form("h%s_DeltaRvsDPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDeltaRVsDpt);
    }

    if (deltaEtaAxis >= 0) {
      TH1* hDeltaEta = fDmesons->Projection(deltaEtaAxis, "EO");
      hDeltaEta->SetName(Form("h%s_DeltaEta_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDeltaEta);
    }

    if (deltaPhiAxis >= 0) {
      TH1* hDeltaPhi = fDmesons->Projection(deltaPhiAxis, "EO");
      hDeltaPhi->SetName(Form("h%s_DeltaPhi_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDeltaPhi);
    }

    if (jetPtAxis >= 0) {
      TH1* hJetPt = fDmesons->Projection(jetPtAxis, "EO");
      hJetPt->SetName(Form("h%s_JetPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hJetPt);
      
      TH2* hDPtVsJetPt = fDmesons->Projection(dPtAxis, jetPtAxis, "EO");
      hDPtVsJetPt->SetName(Form("h%s_MesonPtvsJetPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDPtVsJetPt);
    }

    if (jetEtaAxis >= 0 && jetPhiAxis >= 0) {
      TH2* hJetpos = fDmesons->Projection(jetPhiAxis, jetEtaAxis, "EO");
      hJetpos->SetName(Form("h%s_JetPhiVsEta_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hJetpos);
    }
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::OpenInputFile()
{
  // Open the input file and set the input list.
  
  TString fname;
  
  fname = fInputPath;
  fname += "/";
  fname += fTrainName;
  fname += "/";
  fname += fInputFileName;

  if (!fInputFile) {
    Printf("Info-DJetCorrAnalysis::OpenInputFile : Opening file '%s'", fname.Data()); 
    fInputFile = TFile::Open(fname);
  }
  else {
    Printf("Info-DJetCorrAnalysis::OpenInputFile : File '%s' already open", fInputFile->GetName()); 
  }

  if (!fInputFile || fInputFile->IsZombie()) {
    Printf("Error-DJetCorrAnalysis::OpenInputFile : Could not open file '%s'", fname.Data()); 
    fInputFile = 0;
    return kFALSE;
  }

  if (!fInputDirectoryFile || fInputDirFileName != fInputDirectoryFile->GetName()) {
    ClearInputData();
    
    delete fInputDirectoryFile;
    Printf("Info-DJetCorrAnalysis::OpenInputFile : Getting directory '%s' from file '%s'", fInputDirFileName.Data(), fInputFile->GetName()); 
    fInputDirectoryFile = dynamic_cast<TDirectoryFile*>(fInputFile->Get(fInputDirFileName));

    delete fInputList;
    fInputList = 0;
  }

  if (!fInputDirectoryFile) {
    Printf("Error-DJetCorrAnalysis::OpenInputFile : Could not get directory '%s' from file '%s'", fInputDirFileName.Data(), fInputFile->GetName()); 
    return kFALSE;
  }

  Printf("Info-DJetCorrAnalysis::OpenInputFile : Success.");

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::LoadInputList(const char* inputListName)
{
  // Load the input list.

  if (!OpenInputFile()) return kFALSE;
  
  if (!fInputList || strcmp(inputListName, fInputList->GetName()) != 0) {
    ClearInputData();
    
    delete fInputList;
    Printf("Info-DJetCorrAnalysis::OpenInputFile : Getting list '%s' from directory '%s' of file '%s'", inputListName, fInputDirectoryFile->GetName(), fInputFile->GetName()); 
    fInputList = dynamic_cast<TList*>(fInputDirectoryFile->Get(inputListName));
  }
  
  if (!fInputList) {
    Printf("Error-DJetCorrAnalysis::OpenInputFile : Could not get list '%s' from directory '%s' of file '%s'", inputListName, fInputDirectoryFile->GetName(), fInputFile->GetName());
    return kFALSE;
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::LoadQAList()
{
  // Load the input QA list.

  if (!OpenInputFile()) return kFALSE;
  
  if (!fInputQAList || fQAListName != fInputQAList->GetName()) {
    delete fInputQAList;
    Printf("Info-DJetCorrAnalysis::OpenInputFile : Getting list '%s' from file '%s'", fQAListName.Data(), fInputFile->GetName()); 
    fInputQAList = dynamic_cast<TList*>(fInputFile->Get(fQAListName));
  }
  
  if (!fInputQAList) {
    Printf("Error-DJetCorrAnalysis::OpenInputFile : Could not get list '%s' from file '%s'", fQAListName.Data(), fInputFile->GetName());
    return kFALSE;
  }
    
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::LoadTHnSparse()
{
  if (!fDmesons) {
    fDmesons = static_cast<THnSparse*>(fInputList->FindObject("fDmesons"));
  }
  if (!fDmesons) {
    Printf("Error-DJetCorrAnalysis::LoadTHnSparse : could not open find THnSparse 'fDmesons'"); 
    return kFALSE;
  }

  fTHnSparseMapGenerated = kFALSE;

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrAnalysis::CloseInputFile()
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
Bool_t DJetCorrAnalysis::SaveOutputFile()
{
  // Save the content of fOutputList in a file.
  
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
    Printf("Info-DJetCorrAnalysis::SaveOutputFile : creating directory '%s'.", path.Data()); 
    gSystem->mkdir(path.Data(), kTRUE);
  }
  
  TFile* outputFile = TFile::Open(fname, opt);

  if (!outputFile || outputFile->IsZombie()) {
    Printf("Error-DJetCorrAnalysis::SaveOutputFile : could not open file '%s' to write (if file exists, you may need to set the overwrite option).", fname.Data()); 
    outputFile = 0;
    return kFALSE;
  }

  outputFile->cd();

  Printf("Info-DJetCorrAnalysis::SaveOutputFile : Now streaming results."); 
  fOutputList->Write();

  Printf("Info-DJetCorrAnalysis::SaveOutputFile : Closing the output file."); 
  outputFile->Close();
  delete outputFile;
  outputFile = 0;
  
  return kTRUE;
}
