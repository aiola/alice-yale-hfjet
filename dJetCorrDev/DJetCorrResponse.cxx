// This class projects the THnSparse results of the AliJetResponseMaker into TH1/TH2/TH3
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

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
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>

#include "DJetCorrAnalysisParams.h"

#include "DJetCorrResponse.h"

ClassImp(DJetCorrResponse);

//____________________________________________________________________________________
DJetCorrResponse::DJetCorrResponse() :
  DJetCorrBase(),
  fHistMatching(0),
  fHistJets1(0),
  fHistJets2(0)
{
  // Default ctr.

  fAnaType = DJetCorrAnalysisParams::kResponseMatrixAna;
}

//____________________________________________________________________________________
DJetCorrResponse::DJetCorrResponse(const char* train, const char* path) :
  DJetCorrBase(train, path),
  fHistMatching(0),
  fHistJets1(0),
  fHistJets2(0)
{
  // Standard ctr.

  fAnaType = DJetCorrAnalysisParams::kResponseMatrixAna;
  fInputFileName = "AnalysisResults.root";
  fInputDirFileName = "";
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::LoadTHnSparse()
{
  // Loads the THnSparse histograms.

  Printf("Info-DJetCorrResponse::LoadTHnSparse : Entering function");

  if (!fInputList) return kFALSE;
  
  if (!fHistMatching) {
    fHistMatching = static_cast<THnSparse*>(fInputList->FindObject("fHistMatching"));
    if (!fHistMatching) {
      Printf("Error-DJetCorrResponse::LoadTHnSparse : could not open find THnSparse 'fHistMatching'"); 
      return kFALSE;
    }
  }
  
  if (!fHistJets1) {
    fHistJets1 = static_cast<THnSparse*>(fInputList->FindObject("fHistJets1"));
    if (!fHistJets1) {
      Printf("Error-DJetCorrResponse::LoadTHnSparse : could not open find THnSparse 'fHistJets1'"); 
      return kFALSE;
    }
  }
  
  if (!fHistJets2) {
    fHistJets2 = static_cast<THnSparse*>(fInputList->FindObject("fHistJets2"));
    if (!fHistJets2) {
      Printf("Error-DJetCorrResponse::LoadTHnSparse : could not open find THnSparse 'fHistJets2'"); 
      return kFALSE;
    }
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ClearInputData()
{
  // Clear the input data.

  if (!DJetCorrBase::ClearInputData()) return kFALSE;

  fHistMatching = 0;
  fHistJets1 = 0;
  fHistJets2 = 0;

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseMatrices()
{
  // Generate D-jet correlation histograms.

  DJetCorrAnalysisParams* params = 0;
  TIter next(fAnalysisParams);

  while ((params = static_cast<DJetCorrAnalysisParams*>(next()))) ProjectResponseMatrices(params);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseMatrices(Int_t i)
{
  // Generate D-jet correlation histograms.

  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(i));

  if (params) return ProjectResponseMatrices(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseMatrices(const char* paramName)
{
  // Generate D-jet correlation histograms.

  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->FindObject(paramName));

  if (params) return ProjectResponseMatrices(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseMatrices(DJetCorrAnalysisParams* params)
{
  // Generate D-jet correlation histograms.

  Printf("Info-DJetCorrResponse::ProjectResponseMatrices : Entering function");

  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;
  }
  
  result = LoadInputList(params->GetInputListName());
  if (!result) return kFALSE;

  if (!LoadTHnSparse()) {
    return kFALSE;
  }

  Printf("Info-DJetCorrResponse::ProjectResponseMatrices : Start projections");
  
  ProjectResponseMatrix(params);
  ProjectResponseMatricesVsJetPt(params);
  ProjectResponseMatricesVsZ(params);
  
  TH1::AddDirectory(addDirStatus);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseMatricesVsJetPt(DJetCorrAnalysisParams* params)
{
  for (Int_t i = 0; i < params->GetNJetPtBins(); i++) {
    ProjectResponseZMatrix(params, params->GetJetPtBin(i), params->GetJetPtBin(i+1));
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseMatricesVsZ(DJetCorrAnalysisParams* params)
{
  for (Int_t i = 0; i < params->GetNzBins(); i++) {
    ProjectResponseJetPtMatrix(params, params->GetzBin(i), params->GetzBin(i+1));
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseMatrix(DJetCorrAnalysisParams* params)
{
  TString hname;
  TString htitle;

  if (!fHistMatching || !fHistJets2) return kFALSE;

  Int_t jetPt1Axis = GetAxisIndex("p_{T,1}", fHistMatching, kTRUE);
  Int_t jetPt2Axis = GetAxisIndex("p_{T,2}", fHistMatching, kTRUE);
  Int_t z1Axis = GetAxisIndex("z_{flavour,1}", fHistMatching, kTRUE);
  Int_t z2Axis = GetAxisIndex("z_{flavour,2}", fHistMatching, kTRUE);
  Int_t dPt1Axis = GetAxisIndex("p_{T,1}^{D}", fHistMatching, kTRUE);
  Int_t dPt2Axis = GetAxisIndex("p_{T,2}^{D}", fHistMatching, kTRUE);

  Int_t jetPtPartAxis = GetAxisIndex("p_{T}", fHistJets2, kTRUE);
  Int_t zPartAxis = GetAxisIndex("z_{flavour}", fHistJets2, kTRUE);
  Int_t dPtPartAxis = GetAxisIndex("p_{T}^{D}", fHistJets2, kTRUE);

  if (jetPt1Axis < 0 || jetPt2Axis < 0 || z1Axis < 0 || z2Axis < 0 || dPt1Axis < 0 || dPt2Axis < 0 ||
      jetPtPartAxis < 0 || zPartAxis < 0 || dPtPartAxis < 0) return kFALSE;

  Int_t minzBin = fHistMatching->GetAxis(z2Axis)->FindBin(params->GetMinZ());
  Int_t maxzBin = fHistMatching->GetAxis(z2Axis)->FindBin(params->GetMaxZ());

  if (fHistMatching->GetAxis(z2Axis)->GetBinUpEdge(minzBin) <= params->GetMinZ()) minzBin++;
  if (fHistMatching->GetAxis(z2Axis)->GetBinLowEdge(maxzBin) >= params->GetMaxZ()) maxzBin--;
  
  if (fHistMatching->GetAxis(z2Axis)->GetBinUpEdge(maxzBin) == 1.0) { // if the up edge == 1.0 includes the next bin
    maxzBin++;
  }

  if (fHistMatching->GetAxis(z2Axis)->GetBinLowEdge(minzBin) == 1.0) { // if the low edge == 1.0 excludes the first bin
    minzBin++;
  }

  // Apply the axis cuts
  fHistMatching->GetAxis(jetPt1Axis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(jetPt2Axis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(z1Axis)->SetRange(minzBin, maxzBin);
  fHistMatching->GetAxis(z2Axis)->SetRange(minzBin, maxzBin);
  fHistMatching->GetAxis(dPt1Axis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());
  fHistMatching->GetAxis(dPt2Axis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());

  fHistJets2->GetAxis(jetPtPartAxis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistJets2->GetAxis(zPartAxis)->SetRange(minzBin, maxzBin);
  fHistJets2->GetAxis(dPtPartAxis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());

  // Project 4D matrix
  hname = Form("%s_ResponseMatrix_JetPt_Z", params->GetName());
  htitle = Form("Response matrix");
  Int_t dims[4] = {jetPt1Axis, z1Axis, jetPt2Axis, z2Axis};
  THnSparse* resp = fHistMatching->Projection(4, dims, "");

  resp->SetName(hname);
  resp->SetTitle(htitle);

  hname = Form("%s_Efficiency_JetPt_Z", params->GetName());
  htitle = Form("Efficiency");
  TH2* eff = resp->Projection(3, 2);

  eff->SetName(hname);
  eff->SetTitle(htitle);
  eff->GetXaxis()->SetTitle("#it{p}_{jet}^{part} GeV/#it{c}");
  eff->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  eff->GetZaxis()->SetTitle("Efficiency");
  
  TH2* truth_proj = fHistJets2->Projection(zPartAxis, jetPtPartAxis, "");

  eff->Divide(truth_proj);
  delete truth_proj;
  truth_proj = 0;

  fOutputList->Add(resp);
  fOutputList->Add(eff);

  // Project 2D matrix for D meson pT
  hname = Form("%s_ResponseMatrix_DpT", params->GetName());
  htitle = Form("Response matrix for D meson #it{p}_{T}");
  TH2* respDpT = fHistMatching->Projection(dPt2Axis, dPt1Axis, "");

  respDpT->SetName(hname);
  respDpT->SetTitle(htitle);

  hname = Form("%s_Efficiency_DpT", params->GetName());
  htitle = Form("Efficiency vs D meson #it{p}_{T}");
  TH1* effDpT = respDpT->ProjectionY(hname);

  effDpT->SetTitle(htitle);
  effDpT->GetXaxis()->SetTitle("#it{p}_{T,D}^{part} GeV/#it{c}");
  effDpT->GetYaxis()->SetTitle("Efficiency");
  
  TH1* truth_projDpT = fHistJets2->Projection(dPtPartAxis, "");

  effDpT->Divide(truth_projDpT);
  delete truth_projDpT;
  truth_projDpT = 0;

  fOutputList->Add(respDpT);
  fOutputList->Add(effDpT);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseJetPtMatrix(DJetCorrAnalysisParams* params, Double_t minZ, Double_t maxZ)
{
  TString hname;
  TString htitle;

  if (!fHistMatching || !fHistJets2) return kFALSE;

  Int_t jetPt1Axis = GetAxisIndex("p_{T,1}", fHistMatching, kTRUE);
  Int_t jetPt2Axis = GetAxisIndex("p_{T,2}", fHistMatching, kTRUE);
  Int_t z1Axis = GetAxisIndex("z_{flavour,1}", fHistMatching, kTRUE);
  Int_t z2Axis = GetAxisIndex("z_{flavour,2}", fHistMatching, kTRUE);
  Int_t dPt1Axis = GetAxisIndex("p_{T,1}^{D}", fHistMatching, kTRUE);
  Int_t dPt2Axis = GetAxisIndex("p_{T,2}^{D}", fHistMatching, kTRUE);

  Int_t jetPtPartAxis = GetAxisIndex("p_{T}", fHistJets2, kTRUE);
  Int_t zPartAxis = GetAxisIndex("z_{flavour}", fHistJets2, kTRUE);
  Int_t dPtPartAxis = GetAxisIndex("p_{T}^{D}", fHistJets2, kTRUE);

  if (jetPt1Axis < 0 || jetPt2Axis < 0 || z1Axis < 0 || z2Axis < 0 || dPt1Axis < 0 || dPt2Axis < 0 ||
      jetPtPartAxis < 0 || zPartAxis < 0 || dPtPartAxis < 0) return kFALSE;

  Int_t minzBin = fHistMatching->GetAxis(z2Axis)->FindBin(minZ);
  Int_t maxzBin = fHistMatching->GetAxis(z2Axis)->FindBin(maxZ);

  if (fHistMatching->GetAxis(z2Axis)->GetBinUpEdge(minzBin) <= minZ) minzBin++;
  if (fHistMatching->GetAxis(z2Axis)->GetBinLowEdge(maxzBin) >= maxZ) maxzBin--;
  
  if (fHistMatching->GetAxis(z2Axis)->GetBinUpEdge(maxzBin) == 1.0) { // if the up edge == 1.0 includes the next bin
    maxzBin++;
  }

  if (fHistMatching->GetAxis(z2Axis)->GetBinLowEdge(minzBin) == 1.0) { // if the low edge == 1.0 excludes the first bin
    minzBin++;
  }

  fHistMatching->GetAxis(jetPt1Axis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(jetPt2Axis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(z1Axis)->SetRange(minzBin, maxzBin);
  fHistMatching->GetAxis(z2Axis)->SetRange(minzBin, maxzBin);
  fHistMatching->GetAxis(dPt1Axis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());
  fHistMatching->GetAxis(dPt2Axis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());

  fHistJets2->GetAxis(jetPtPartAxis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistJets2->GetAxis(zPartAxis)->SetRange(minzBin, maxzBin);
  fHistJets2->GetAxis(dPtPartAxis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());

  hname = Form("%s_ResponseMatrix_JetPt_Z_%d_%d", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100));
  htitle = Form("Response matrix for jet #it{p}_{T}: %.2f < #it{z}_{||}^{part} < %.2f", minZ, maxZ);
  TH2* resp = fHistMatching->Projection(1, 0, "");

  resp->SetName(hname);
  resp->SetTitle(htitle);
  resp->GetXaxis()->SetTitle("#it{p}_{jet}^{det} GeV/#it{c}");
  resp->GetYaxis()->SetTitle("#it{p}_{jet}^{part} GeV/#it{c}");
  resp->Rebin2D(10,10);

  hname = Form("%s_Efficiency_JetPt_Z_%d_%d", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100));
  htitle = Form("Efficiency for jet #it{p}_{T}: %.2f < #it{z}_{||}^{part} < %.2f", minZ, maxZ);
  TH1* eff = resp->ProjectionY(hname);
  
  eff->SetTitle(htitle);
  eff->GetXaxis()->SetTitle("#it{p}_{jet}^{part} GeV/#it{c}");
  eff->GetYaxis()->SetTitle("Efficiency");
  
  TH1* truth_proj = fHistJets2->Projection(2, "");
  truth_proj->Rebin(10);
  eff->Divide(truth_proj);
  delete truth_proj;
  truth_proj = 0;

  fOutputList->Add(resp);
  fOutputList->Add(eff);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseZMatrix(DJetCorrAnalysisParams* params, Double_t minJetPt, Double_t maxJetPt)
{
  TString hname;
  TString htitle;

  if (!fHistMatching || !fHistJets2) return kFALSE;

  Int_t jetPt1Axis = GetAxisIndex("p_{T,1}", fHistMatching, kTRUE);
  Int_t jetPt2Axis = GetAxisIndex("p_{T,2}", fHistMatching, kTRUE);
  Int_t z1Axis = GetAxisIndex("z_{flavour,1}", fHistMatching, kTRUE);
  Int_t z2Axis = GetAxisIndex("z_{flavour,2}", fHistMatching, kTRUE);
  Int_t dPt1Axis = GetAxisIndex("p_{T,1}^{D}", fHistMatching, kTRUE);
  Int_t dPt2Axis = GetAxisIndex("p_{T,2}^{D}", fHistMatching, kTRUE);

  Int_t jetPtPartAxis = GetAxisIndex("p_{T}", fHistJets2, kTRUE);
  Int_t zPartAxis = GetAxisIndex("z_{flavour}", fHistJets2, kTRUE);
  Int_t dPtPartAxis = GetAxisIndex("p_{T}^{D}", fHistJets2, kTRUE);

  if (jetPt1Axis < 0 || jetPt2Axis < 0 || z1Axis < 0 || z2Axis < 0 || dPt1Axis < 0 || dPt2Axis < 0 ||
      jetPtPartAxis < 0 || zPartAxis < 0 || dPtPartAxis < 0) return kFALSE;

  Int_t minzBin = fHistMatching->GetAxis(z2Axis)->FindBin(params->GetMinZ());
  Int_t maxzBin = fHistMatching->GetAxis(z2Axis)->FindBin(params->GetMaxZ());

  if (fHistMatching->GetAxis(z2Axis)->GetBinUpEdge(minzBin) <= params->GetMinZ()) minzBin++;
  if (fHistMatching->GetAxis(z2Axis)->GetBinLowEdge(maxzBin) >= params->GetMaxZ()) maxzBin--;
  
  if (fHistMatching->GetAxis(z2Axis)->GetBinUpEdge(maxzBin) == 1.0) { // if the up edge == 1.0 includes the next bin
    maxzBin++;
  }

  if (fHistMatching->GetAxis(z2Axis)->GetBinLowEdge(minzBin) == 1.0) { // if the low edge == 1.0 excludes the first bin
    minzBin++;
  }

  fHistMatching->GetAxis(jetPt1Axis)->SetRangeUser(minJetPt, maxJetPt);
  fHistMatching->GetAxis(jetPt2Axis)->SetRangeUser(minJetPt, maxJetPt);
  fHistMatching->GetAxis(z1Axis)->SetRange(minzBin, maxzBin);
  fHistMatching->GetAxis(z2Axis)->SetRange(minzBin, maxzBin);
  fHistMatching->GetAxis(dPt1Axis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());
  fHistMatching->GetAxis(dPt2Axis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());

  fHistJets2->GetAxis(jetPtPartAxis)->SetRangeUser(minJetPt, maxJetPt);
  fHistJets2->GetAxis(zPartAxis)->SetRange(minzBin, maxzBin);
  fHistJets2->GetAxis(dPtPartAxis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());

  hname = Form("%s_ResponseMatrix_Z_JetPt_%d_%d", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt));
  htitle = Form("Response matrix for #it{z}_{||}: %.1f < #it{p}_{T,jet}^{part} < %.1f GeV/#it{c}", minJetPt, maxJetPt);
  TH2* resp = fHistMatching->Projection(z2Axis, z1Axis, "");

  resp->SetName(hname);
  resp->SetTitle(htitle);
  resp->GetXaxis()->SetTitle("#it{z}_{||}^{det}");
  resp->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  resp->Rebin2D(5,5);

  hname = Form("%s_Efficiency_Z_JetPt_%d_%d", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt));
  htitle = Form("Efficiency for #it{z}_{||}: %.1f < #it{p}_{T,jet}^{part} < %.1f GeV/#it{c}", minJetPt, maxJetPt);
  TH1* eff = resp->ProjectionY(hname);
  
  eff->SetTitle(htitle);
  eff->GetXaxis()->SetTitle("#it{z}_{||}^{part}");
  eff->GetYaxis()->SetTitle("Efficiency");
  
  TH1* truth_proj = fHistJets2->Projection(zPartAxis, "");
  truth_proj->Rebin(5);
  eff->Divide(truth_proj);
  delete truth_proj;
  truth_proj = 0;

  fOutputList->Add(resp);
  fOutputList->Add(eff);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrices()
{
  // Plot D-jet correlation histograms.

  DJetCorrAnalysisParams* params = 0;
  TIter next(fAnalysisParams);

  while ((params = static_cast<DJetCorrAnalysisParams*>(next()))) PlotResponseMatrices(params);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrices(Int_t i)
{
  // Plot D-jet correlation histograms.

  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(i));

  if (params) return PlotResponseMatrices(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrices(const char* paramName)
{
  // Plot D-jet correlation histograms.

  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->FindObject(paramName));

  if (params) return PlotResponseMatrices(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrices(DJetCorrAnalysisParams* params)
{
  // Plot D-jet correlation histograms.

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  Bool_t result = kFALSE;
  
  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;

    result = LoadOutputHistograms();
    if (!result) return kFALSE;
  }

  PlotResponseMatrix(params);
  PlotResponseMatricesVsJetPt(params);
  PlotResponseMatricesVsZ(params);
  PlotResponseDPtMatrix(params);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatricesVsJetPt(DJetCorrAnalysisParams* params)
{
  for (Int_t i = 0; i < params->GetNJetPtBins(); i++) {
    PlotResponseZMatrix(params, params->GetJetPtBin(i), params->GetJetPtBin(i+1));
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatricesVsZ(DJetCorrAnalysisParams* params)
{
  for (Int_t i = 1; i < params->GetNzBins(); i++) {
    PlotResponseJetPtMatrix(params, params->GetzBin(i), params->GetzBin(i+1));
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrix(DJetCorrAnalysisParams* params)
{
  TString hname;

  hname = Form("%s_Efficiency_JetPt_Z", params->GetName());
  TH2* eff = static_cast<TH2*>(fOutputList->FindObject(hname));
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  TCanvas* canvasEff = SetUpCanvas(hname,
                                   eff->GetXaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                                   eff->GetYaxis()->GetTitle(), 0., 1., kFALSE);
  eff->Draw("colz same");
  gPad->SetLogz();
  if (fSavePlots) SavePlot(canvasEff);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseJetPtMatrix(DJetCorrAnalysisParams* params, Double_t minZ, Double_t maxZ)
{
  TString hname;
  
  hname = Form("%s_ResponseMatrix_JetPt_Z_%d_%d", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100));
  TH2* resp = static_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) {
    Printf("Error-DJetCorrResponse::PlotResponseJetPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }
  
  TCanvas *canvasResp = SetUpCanvas(resp->GetName(),
                                    resp->GetXaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                                    resp->GetYaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE);
  resp->Draw("colz same");
  gPad->SetLogz();
  if (fSavePlots) SavePlot(canvasResp);
  
  hname = Form("%s_Efficiency_JetPt_Z_%d_%d", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100));
  TH1* eff = static_cast<TH1*>(fOutputList->FindObject(hname));
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseJetPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }
  
  TCanvas* canvasEff = SetUpCanvas(eff->GetName(),
                                   eff->GetXaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                                   eff->GetYaxis()->GetTitle(), 0., eff->GetMaximum()*1.5, kFALSE);
  eff->Draw("same");
  if (fSavePlots) SavePlot(canvasEff);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseZMatrix(DJetCorrAnalysisParams* params, Double_t minJetPt, Double_t maxJetPt)
{
  TString hname;

  hname = Form("%s_ResponseMatrix_Z_JetPt_%d_%d", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt));
  TH2* resp = static_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) {
    Printf("Error-DJetCorrResponse::PlotResponseZMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }
  
  TCanvas *canvasResp = SetUpCanvas(resp->GetName(),
                                    resp->GetXaxis()->GetTitle(), 0., 1.0, kFALSE,
                                    resp->GetYaxis()->GetTitle(), 0., 1.0, kFALSE);
  resp->Draw("colz same");
  gPad->SetLogz();
  if (fSavePlots) SavePlot(canvasResp);
  
  hname = Form("%s_Efficiency_Z_JetPt_%d_%d", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt));
  TH1* eff = static_cast<TH1*>(fOutputList->FindObject(hname));
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseZMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }
  
  TCanvas* canvasEff = SetUpCanvas(eff->GetName(),
                                   eff->GetXaxis()->GetTitle(), 0., 1.0, kFALSE,
                                   eff->GetYaxis()->GetTitle(), 0., eff->GetMaximum()*1.5, kFALSE);
  eff->Draw("same");
  if (fSavePlots) SavePlot(canvasEff);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseDPtMatrix(DJetCorrAnalysisParams* params)
{
  TString hname;
  
  hname = Form("%s_ResponseMatrix_DpT", params->GetName());
  TH2* resp = static_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) {
    Printf("Error-DJetCorrResponse::PlotResponseDPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }
  
  TCanvas *canvasResp = SetUpCanvas(resp->GetName(),
                                    resp->GetXaxis()->GetTitle(), params->GetMinDPt(), params->GetMaxDPt(), kFALSE,
                                    resp->GetYaxis()->GetTitle(), params->GetMinDPt(), params->GetMaxDPt(), kFALSE);
  resp->Draw("colz same");
  gPad->SetLogz();
  if (fSavePlots) SavePlot(canvasResp);
  
  hname = Form("%s_Efficiency_DpT", params->GetName());
  TH1* eff = static_cast<TH1*>(fOutputList->FindObject(hname));
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseDPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }
  
  TCanvas* canvasEff = SetUpCanvas(eff->GetName(),
                                   eff->GetXaxis()->GetTitle(), params->GetMinDPt(), params->GetMaxDPt(), kFALSE,
                                   eff->GetYaxis()->GetTitle(), 0., eff->GetMaximum()*1.5, kFALSE);
  eff->Draw("same");
  if (fSavePlots) SavePlot(canvasEff);
  
  return kTRUE;
}
