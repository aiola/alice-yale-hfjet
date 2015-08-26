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
#include <TGraphAsymmErrors.h>

#include "DJetCorrAnalysisParams.h"
#include "HistoStyler.h"

#include "DJetCorrResponse.h"

ClassImp(DJetCorrResponse);

//____________________________________________________________________________________
DJetCorrResponse::DJetCorrResponse() :
  DJetCorrBase(),
  fEfficiencyMode("b(1,1) mode"),
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
  fEfficiencyMode("b(1,1) mode"),
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
  
  ProjectResponseDPtMatrix(params);
  ProjectResponseMatrices4DVsDPt(params);
  ProjectResponseZMatricesVsDPt(params);
  ProjectResponseZMatricesVsJetPt(params);
  ProjectResponseJetPtMatricesVsDPt(params);
  ProjectResponseJetPtMatricesVsZ(params);
  
  TH1::AddDirectory(addDirStatus);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseMatrices4DVsDPt(DJetCorrAnalysisParams* params)
{
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt++) {
    ProjectResponseMatrix4D(params, idpt);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseZMatricesVsDPt(DJetCorrAnalysisParams* params)
{
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt++) {
    ProjectResponseZMatrix(params, -1, idpt);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseZMatricesVsJetPt(DJetCorrAnalysisParams* params)
{
  for (Int_t ijetpt = 0; ijetpt < params->GetNJetPtBins(); ijetpt++) {
    ProjectResponseZMatrix(params, ijetpt, -1);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseJetPtMatricesVsDPt(DJetCorrAnalysisParams* params)
{
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt++) {
    ProjectResponseJetPtMatrix(params, -1, idpt);
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseJetPtMatricesVsZ(DJetCorrAnalysisParams* params)
{
  for (Int_t iz = 0; iz < params->GetNzBins(); iz++) {
    ProjectResponseJetPtMatrix(params, iz, -1);
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseMatrix4D(DJetCorrAnalysisParams* params, Int_t dptBin)
{
  TString hname;
  TString htitle;

  if (!fHistMatching || !fHistJets2) return kFALSE;

  Double_t minDPt = 0;
  Double_t maxDPt = 0;

  if (dptBin >= 0) {
    minDPt = params->GetDPtBin(dptBin);
    maxDPt = params->GetDPtBin(dptBin+1);
  }
  else {
    minDPt = params->GetMinDPt();
    maxDPt = params->GetMaxDPt();
  }

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

  // Apply the axis cuts
  fHistMatching->GetAxis(jetPt1Axis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(jetPt2Axis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(z1Axis)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());
  fHistMatching->GetAxis(z2Axis)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());
  fHistMatching->GetAxis(dPt1Axis)->SetRangeUser(minDPt, maxDPt);
  fHistMatching->GetAxis(dPt2Axis)->SetRangeUser(minDPt, maxDPt);

  fHistJets2->GetAxis(jetPtPartAxis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistJets2->GetAxis(zPartAxis)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());
  fHistJets2->GetAxis(dPtPartAxis)->SetRangeUser(minDPt, maxDPt);

  // Project 4D matrix
  hname = Form("%s_ResponseMatrix_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);
  htitle = Form("Response matrix");
  Int_t dims[4] = {jetPt1Axis, z1Axis, jetPt2Axis, z2Axis};
  THnSparse* resp = fHistMatching->Projection(4, dims, "");
  resp->SetName(hname);
  resp->SetTitle(htitle);

  hname = Form("%s_Measured_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);
  TH2* measured = fHistMatching->Projection(z1Axis, jetPt1Axis);
  measured->SetName(hname);
  measured->SetTitle("Measured");
  measured->GetXaxis()->SetTitle("#it{p}_{T,jet}^{det} GeV/#it{c}");
  measured->GetYaxis()->SetTitle("#it{z}_{||}^{det}");
  measured->GetZaxis()->SetTitle("Counts");

  hname = Form("%s_TruthReco_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);
  TH2* truthReco = fHistMatching->Projection(z2Axis, jetPt2Axis);
  truthReco->SetName(hname);
  truthReco->SetTitle("Truth reconstructed");
  truthReco->GetXaxis()->SetTitle("#it{p}_{T,jet}^{part} GeV/#it{c}");
  truthReco->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  truthReco->GetZaxis()->SetTitle("Counts");

  hname = Form("%s_Truth_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);
  TH2* truth = fHistJets2->Projection(zPartAxis, jetPtPartAxis, "");
  truth->SetName(hname);
  truth->SetTitle("Truth");
  truth->GetXaxis()->SetTitle("#it{p}_{T,jet}^{part} GeV/#it{c}");
  truth->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  truth->GetZaxis()->SetTitle("Counts");

  // Remove cuts at detector level to extract the missed events due to the kinematical restrictions
  fHistMatching->GetAxis(jetPt1Axis)->SetRange(0, fHistMatching->GetAxis(jetPt1Axis)->GetNbins()+1);
  fHistMatching->GetAxis(z1Axis)->SetRangeUser(0, fHistMatching->GetAxis(z1Axis)->GetNbins()+1);
  fHistMatching->GetAxis(dPt1Axis)->SetRangeUser(0, fHistMatching->GetAxis(dPt1Axis)->GetNbins()+1);

  hname = Form("%s_TruthRecoAll_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);
  TH2* truthRecoAll = fHistMatching->Projection(z2Axis, jetPt2Axis);
  truthRecoAll->SetName(hname);
  truthRecoAll->SetTitle("Truth reconstructed (all)");
  truthRecoAll->GetXaxis()->SetTitle("#it{p}_{T,jet}^{part} GeV/#it{c}");
  truthRecoAll->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  truthRecoAll->GetZaxis()->SetTitle("Counts");

  hname = Form("%s_Efficiency_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);  
  TH2* eff = static_cast<TH2*>(truthReco->Clone(hname));
  eff->Divide(truth);
  eff->SetTitle("Efficiecny");
  eff->GetXaxis()->SetTitle("#it{p}_{T,jet}^{part} GeV/#it{c}");
  eff->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  eff->GetZaxis()->SetTitle("Efficiency");

  hname = Form("%s_Misses_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);  
  TH2* misses = static_cast<TH2*>(truth->Clone(hname));
  misses->Add(truthReco, -1);
  misses->SetTitle("Misses = Truth - Truth reconstructed");
  misses->GetXaxis()->SetTitle("#it{p}_{T,jet}^{part} GeV/#it{c}");
  misses->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  misses->GetZaxis()->SetTitle("Misses");

  hname = Form("%s_KinEfficiency_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);  
  TH2* kinEff = static_cast<TH2*>(truthReco->Clone(hname));
  kinEff->Divide(truthRecoAll);
  kinEff->SetTitle("Kinematic efficiency");
  kinEff->GetXaxis()->SetTitle("#it{p}_{T,jet}^{part} GeV/#it{c}");
  kinEff->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  kinEff->GetZaxis()->SetTitle("Efficiency");

  hname = Form("%s_KinMisses_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);  
  TH2* kinMisses = static_cast<TH2*>(truthRecoAll->Clone(hname));
  kinMisses->Add(truthReco, -1);
  kinMisses->SetTitle("Misses = Truth reconstructed (all) - Truth reconstructed");
  kinMisses->GetXaxis()->SetTitle("#it{p}_{T,jet}^{part} GeV/#it{c}");
  kinMisses->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  kinMisses->GetZaxis()->SetTitle("Misses");

  // Rebin to the coarse binning ready for unfolding
  hname = Form("%s_ResponseMatrix_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), minDPt, maxDPt);
  Int_t nbins[4] = {params->GetNJetPtBins(), params->GetNzBins(), params->GetNJetPtBins(), params->GetNzBins()};
  const Double_t* bins[4] = {params->GetJetPtBins(), params->GetzBins(), params->GetJetPtBins(), params->GetzBins()};
  THnSparse* coarseResp = Rebin(resp, hname, nbins, bins);

  hname = Form("%s_Measured_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), minDPt, maxDPt);
  TH2* coarseMeasured = Rebin(measured, hname, params->GetNJetPtBins(), params->GetJetPtBins(), params->GetNzBins(), params->GetzBins());

  hname = Form("%s_TruthReco_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), minDPt, maxDPt);
  TH2* coarseTruthReco = Rebin(truthReco, hname, params->GetNJetPtBins(), params->GetJetPtBins(), params->GetNzBins(), params->GetzBins());
  
  hname = Form("%s_Truth_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), minDPt, maxDPt);
  TH2* coarseTruth = Rebin(truth, hname, params->GetNJetPtBins(), params->GetJetPtBins(), params->GetNzBins(), params->GetzBins());

  hname = Form("%s_TruthRecoAll_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), minDPt, maxDPt);
  TH2* coarseTruthRecoAll = Rebin(truthRecoAll, hname, params->GetNJetPtBins(), params->GetJetPtBins(), params->GetNzBins(), params->GetzBins());
  
  hname = Form("%s_Efficiency_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), minDPt, maxDPt);
  TH2* coarseEff = static_cast<TH2*>(coarseTruthReco->Clone(hname));
  coarseEff->Divide(coarseTruth);

  hname = Form("%s_Misses_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), minDPt, maxDPt);
  TH2* coarseMisses = static_cast<TH2*>(coarseTruth->Clone(hname));
  coarseMisses->Add(coarseTruthReco, -1);

  hname = Form("%s_KinEfficiency_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), minDPt, maxDPt);
  TH2* coarseKinEff = static_cast<TH2*>(coarseTruthReco->Clone(hname));
  coarseKinEff->Divide(coarseTruthRecoAll);

  hname = Form("%s_KinMisses_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), minDPt, maxDPt);
  TH2* coarseKinMisses = static_cast<TH2*>(coarseTruthRecoAll->Clone(hname));
  coarseKinMisses->Add(coarseTruthReco, -1);

  fOutputList->Add(resp);
  fOutputList->Add(measured);
  fOutputList->Add(truthReco);
  fOutputList->Add(truth);
  fOutputList->Add(truthRecoAll);
  fOutputList->Add(eff);
  fOutputList->Add(misses);
  fOutputList->Add(kinEff);
  fOutputList->Add(kinMisses);

  fOutputList->Add(coarseResp);
  fOutputList->Add(coarseMeasured);
  fOutputList->Add(coarseTruthReco);
  fOutputList->Add(coarseTruth);
  fOutputList->Add(coarseTruthRecoAll);
  fOutputList->Add(coarseEff);
  fOutputList->Add(coarseMisses);
  fOutputList->Add(coarseKinEff);
  fOutputList->Add(coarseKinMisses);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseDPtMatrix(DJetCorrAnalysisParams* params)
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

  // Apply the axis cuts
  fHistMatching->GetAxis(jetPt1Axis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(jetPt2Axis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(z1Axis)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());
  fHistMatching->GetAxis(z2Axis)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());
  fHistMatching->GetAxis(dPt1Axis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());
  fHistMatching->GetAxis(dPt2Axis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());

  fHistJets2->GetAxis(jetPtPartAxis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistJets2->GetAxis(zPartAxis)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());
  fHistJets2->GetAxis(dPtPartAxis)->SetRangeUser(params->GetMinDPt(), params->GetMaxDPt());
  
  // Project 2D matrix for D meson pT
  hname = Form("%s_ResponseMatrix_DpT", params->GetName());
  htitle = Form("Response matrix for D meson #it{p}_{T}");
  TH2* respDpT = fHistMatching->Projection(dPt2Axis, dPt1Axis, "");

  respDpT->SetName(hname);
  respDpT->SetTitle(htitle);

  hname = Form("%s_Efficiency_DpT", params->GetName());
  htitle = Form("Efficiency vs D meson #it{p}_{T}");
  TH1* pass = respDpT->ProjectionY(hname);
  
  TH1* total = fHistJets2->Projection(dPtPartAxis, "");

  TGraphAsymmErrors* effDpT = new TGraphAsymmErrors(pass, total, fEfficiencyMode);
  effDpT->SetName(hname);
  effDpT->SetTitle(htitle);
  effDpT->GetXaxis()->SetTitle("#it{p}_{T,D}^{part} GeV/#it{c}");
  effDpT->GetYaxis()->SetTitle("Efficiency");
  
  delete total;
  delete pass;

  fOutputList->Add(respDpT);
  fOutputList->Add(effDpT);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseJetPtMatrix(DJetCorrAnalysisParams* params, Int_t zBin, Int_t dptBin)
{
  TString hname;
  TString htitle;

  if (!fHistMatching || !fHistJets2) return kFALSE;

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  
  if (dptBin >= 0) {
    minDPt = params->GetDPtBin(dptBin);
    maxDPt = params->GetDPtBin(dptBin+1);
  }
  else {
    minDPt = params->GetMinDPt();
    maxDPt = params->GetMaxDPt();
  }

  Double_t minZ = 0;
  Double_t maxZ = 0;
  
  if (zBin >= 0) {
    minZ = params->GetzBin(zBin);
    maxZ = params->GetzBin(zBin+1);
  }
  else {
    minZ = params->GetMinZ();
    maxZ = params->GetMaxZ();
  }

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

  fHistMatching->GetAxis(jetPt1Axis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(jetPt2Axis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(z1Axis)->SetRangeUser(minZ, maxZ);
  fHistMatching->GetAxis(z2Axis)->SetRangeUser(minZ, maxZ);
  fHistMatching->GetAxis(dPt1Axis)->SetRangeUser(minDPt, maxDPt);
  fHistMatching->GetAxis(dPt2Axis)->SetRangeUser(minDPt, maxDPt);

  fHistJets2->GetAxis(jetPtPartAxis)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistJets2->GetAxis(zPartAxis)->SetRangeUser(minZ, maxZ);
  fHistJets2->GetAxis(dPtPartAxis)->SetRangeUser(minDPt, maxDPt);

  hname = Form("%s_ResponseMatrix_JetPt_Z_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100), minDPt, maxDPt);
  htitle = Form("Response matrix for jet #it{p}_{T}: %.2f < #it{z}_{||}^{part} < %.2f", minZ, maxZ);
  TH2* resp = fHistMatching->Projection(1, 0, "");

  resp->SetName(hname);
  resp->SetTitle(htitle);
  resp->GetXaxis()->SetTitle("#it{p}_{T,jet}^{det} GeV/#it{c}");
  resp->GetYaxis()->SetTitle("#it{p}_{T,jet}^{part} GeV/#it{c}");
  resp->Rebin2D(10,10);

  hname = Form("%s_Efficiency_JetPt_Z_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100), minDPt, maxDPt);
  htitle = Form("%.1f < #it{p}_{T,D}^{part} < %.1f GeV/#it{c} and %.1f < #it{z}_{||}^{part} < %.1f", minDPt, maxDPt, minZ, maxZ);
  TH1* pass = resp->ProjectionY(hname);
  
  TH1* total = fHistJets2->Projection(2, "");
  total->Rebin(10);

  TGraphAsymmErrors* eff = new TGraphAsymmErrors(pass, total, fEfficiencyMode);
  eff->SetName(hname);
  eff->SetTitle(htitle);
  eff->GetXaxis()->SetTitle("#it{p}_{T,jet}^{part} GeV/#it{c}");
  eff->GetYaxis()->SetTitle("Efficiency");

  delete pass;
  delete total;

  fOutputList->Add(resp);
  fOutputList->Add(eff);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseZMatrix(DJetCorrAnalysisParams* params, Int_t jetPtBin, Int_t dptBin)
{
  TString hname;
  TString htitle;

  if (!fHistMatching || !fHistJets2) return kFALSE;

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  
  if (dptBin >= 0) {
    minDPt = params->GetDPtBin(dptBin);
    maxDPt = params->GetDPtBin(dptBin+1);
  }
  else {
    minDPt = params->GetMinDPt();
    maxDPt = params->GetMaxDPt();
  }

  Double_t minJetPt = 0;
  Double_t maxJetPt = 0;
  
  if (jetPtBin >= 0) {
    minJetPt = params->GetJetPtBin(jetPtBin);
    maxJetPt = params->GetJetPtBin(jetPtBin+1);
  }
  else {
    minJetPt = params->GetMinJetPt();
    maxJetPt = params->GetMaxJetPt();
  }

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

  fHistMatching->GetAxis(jetPt1Axis)->SetRangeUser(minJetPt, maxJetPt);
  fHistMatching->GetAxis(jetPt2Axis)->SetRangeUser(minJetPt, maxJetPt);
  fHistMatching->GetAxis(z1Axis)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());
  fHistMatching->GetAxis(z2Axis)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());
  fHistMatching->GetAxis(dPt1Axis)->SetRangeUser(minDPt, maxDPt);
  fHistMatching->GetAxis(dPt2Axis)->SetRangeUser(minDPt, maxDPt);

  fHistJets2->GetAxis(jetPtPartAxis)->SetRangeUser(minJetPt, maxJetPt);
  fHistJets2->GetAxis(zPartAxis)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());
  fHistJets2->GetAxis(dPtPartAxis)->SetRangeUser(minDPt, maxDPt);

  hname = Form("%s_ResponseMatrix_Z_JetPt_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt), minDPt, maxDPt);
  htitle = Form("Response matrix for #it{z}_{||}: %.1f < #it{p}_{T,jet}^{part} < %.1f GeV/#it{c}", minJetPt, maxJetPt);
  TH2* resp = fHistMatching->Projection(z2Axis, z1Axis, "");

  resp->SetName(hname);
  resp->SetTitle(htitle);
  resp->GetXaxis()->SetTitle("#it{z}_{||}^{det}");
  resp->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  resp->Rebin2D(5,5);

  hname = Form("%s_Efficiency_Z_JetPt_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt), minDPt, maxDPt);
  htitle = Form("%.1f < #it{p}_{T,D}^{part} < %.1f GeV/#it{c} and %.1f < #it{p}_{T,jet}^{part} < %.1f GeV/#it{c}", minDPt, maxDPt, minJetPt, maxJetPt);
  TH1* pass = resp->ProjectionY(hname);
  
  TH1* total = fHistJets2->Projection(zPartAxis, "");
  total->Rebin(5);

  TGraphAsymmErrors* eff = new TGraphAsymmErrors(pass, total, fEfficiencyMode);
  eff->SetName(hname);
  eff->SetTitle(htitle);
  eff->GetXaxis()->SetTitle("#it{z}_{||}^{part}");
  eff->GetYaxis()->SetTitle("Efficiency");

  delete pass;
  delete total;

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

  PlotResponseDPtMatrix(params);
  PlotResponseMatrices4DVsDPt(params);
  PlotResponseZMatricesVsDPt(params);
  PlotResponseZMatricesVsJetPt(params);
  PlotResponseJetPtMatricesVsDPt(params);
  PlotResponseJetPtMatricesVsZ(params);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrices4DVsDPt(DJetCorrAnalysisParams* params)
{
  TCanvas* canvasEff = 0;
  
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt+=2) {
    PlotResponseMatrix4D(canvasEff, params, idpt);
    if (fSavePlots) SavePlot(canvasEff);
    canvasEff = 0;
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseZMatricesVsDPt(DJetCorrAnalysisParams* params)
{
  TCanvas* canvasResp = 0;
  TCanvas* canvasEff = 0;
  
  TString cEffName = Form("%s_Efficiency_Z_Vs_DPt", params->GetName());
 
  canvasEff = SetUpCanvas(cEffName,
                          "#it{z}_{||}^{part}", params->GetMinZ(), params->GetMaxZ(), kFALSE,
                          "Efficiency", 0., 1e-6, kFALSE);
  canvasEff->cd();
  TLegend* leg = SetUpLegend(0.1, 0.63, 0.50, 0.89, 16);
  leg->Draw();
  
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt+=2) {
    PlotResponseZMatrix(canvasResp, canvasEff, params, -1, idpt);
    if (fSavePlots) SavePlot(canvasResp);
    canvasResp = 0;
  }

  if (fSavePlots) SavePlot(canvasEff);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseZMatricesVsJetPt(DJetCorrAnalysisParams* params)
{
  TCanvas* canvasResp = 0;
  TCanvas* canvasEff = 0;

  TString cEffName = Form("%s_Efficiency_Z_Vs_JetPt", params->GetName());
 
  canvasEff = SetUpCanvas(cEffName,
                          "#it{z}_{||}^{part}", params->GetMinZ(), params->GetMaxZ(), kFALSE,
                          "Efficiency", 0., 1e-6, kFALSE);

  canvasEff->cd();
  TLegend* leg = SetUpLegend(0.1, 0.63, 0.50, 0.89, 16);
  leg->Draw();
  
  for (Int_t ijetpt = -1; ijetpt < params->GetNJetPtBins(); ijetpt++) {
    PlotResponseZMatrix(canvasResp, canvasEff, params, ijetpt, -1);
    if (fSavePlots) SavePlot(canvasResp);
    canvasResp = 0;
  }

  if (fSavePlots) SavePlot(canvasEff);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseJetPtMatricesVsDPt(DJetCorrAnalysisParams* params)
{
  TCanvas* canvasResp = 0;
  TCanvas* canvasEff = 0;

  TString cEffName = Form("%s_Efficiency_JetPt_Vs_DPt", params->GetName());
 
  canvasEff = SetUpCanvas(cEffName,
                          "#it{p}_{T,jet}^{ch}^{part} GeV/#it{c}", params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                          "Efficiency", 0., 1e-6, kFALSE);

  canvasEff->cd();
  TLegend* leg = SetUpLegend(0.1, 0.63, 0.50, 0.89, 16);
  leg->Draw();
  
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt+=2) {
    PlotResponseJetPtMatrix(canvasResp, canvasEff, params, -1, idpt);
    if (fSavePlots) SavePlot(canvasResp);
    canvasResp = 0;
  }

  if (fSavePlots) SavePlot(canvasEff);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseJetPtMatricesVsZ(DJetCorrAnalysisParams* params)
{
  TCanvas* canvasResp = 0;
  TCanvas* canvasEff = 0;

  TString cEffName = Form("%s_Efficiency_JetPt_Vs_Z", params->GetName());
 
  canvasEff = SetUpCanvas(cEffName,
                          "#it{p}_{T,jet}^{ch}^{part} GeV/#it{c}", params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                          "Efficiency", 0., 1e-6, kFALSE);

  canvasEff->cd();
  TLegend* leg = SetUpLegend(0.1, 0.63, 0.50, 0.89, 16);
  leg->Draw();
  
  for (Int_t iz = -1; iz < params->GetNzBins(); iz++) {
    PlotResponseJetPtMatrix(canvasResp, canvasEff, params, iz, -1);
    if (fSavePlots) SavePlot(canvasResp);
    canvasResp = 0;
  }

  if (fSavePlots) SavePlot(canvasEff);
 
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseDPtMatrix(DJetCorrAnalysisParams* params)
{
  TCanvas *c1 = 0;
  TCanvas *c2 = 0;

  return PlotResponseDPtMatrix(c1, c2, params);
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrix4D(TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t dptBin)
{
  TString hname;

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  
  if (dptBin >= 0) {
    minDPt = params->GetDPtBin(dptBin);
    maxDPt = params->GetDPtBin(dptBin+1);
  }
  else {
    minDPt = params->GetMinDPt();
    maxDPt = params->GetMaxDPt();
  }

  hname = Form("%s_Efficiency_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);
  TH2* eff = dynamic_cast<TH2*>(fOutputList->FindObject(hname));
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  if (!canvasEff && !GetCanvas(hname)) {
    canvasEff = SetUpCanvas(hname,
                            eff->GetXaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                            eff->GetYaxis()->GetTitle(), 0., 1., kFALSE);
  }

  if (canvasEff) {
    canvasEff->cd();
    eff->Draw("colz same");
    eff->GetZaxis()->SetRangeUser(0.,1.);
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseJetPtMatrix(TCanvas*& canvasResp, TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t zBin, Int_t dptBin)
{ 
  TString hname;

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  
  if (dptBin >= 0) {
    minDPt = params->GetDPtBin(dptBin);
    maxDPt = params->GetDPtBin(dptBin+1);
  }
  else {
    minDPt = params->GetMinDPt();
    maxDPt = params->GetMaxDPt();
  }

  Double_t minZ = 0;
  Double_t maxZ = 0;
  
  if (zBin >= 0) {
    minZ = params->GetzBin(zBin);
    maxZ = params->GetzBin(zBin+1);
  }
  else {
    minZ = params->GetMinZ();
    maxZ = params->GetMaxZ();
  }
  
  hname = Form("%s_ResponseMatrix_JetPt_Z_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100), minDPt, maxDPt);
  TH2* resp = dynamic_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) {
    Printf("Error-DJetCorrResponse::PlotResponseJetPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  if (!canvasResp && !GetCanvas(resp->GetName())) {
    canvasResp = SetUpCanvas(resp->GetName(),
                             resp->GetXaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                             resp->GetYaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE);
  }

  if (canvasResp) {
    canvasResp->cd();
    resp->Draw("colz same");
    gPad->SetLogz();
  }
  
  hname = Form("%s_Efficiency_JetPt_Z_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100), minDPt, maxDPt);
  TGraphAsymmErrors* eff = dynamic_cast<TGraphAsymmErrors*>(fOutputList->FindObject(hname));
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseJetPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }
  
  if (!canvasEff && !GetCanvas(eff->GetName())) {
    canvasEff = SetUpCanvas(eff->GetName(),
                            eff->GetXaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                            eff->GetYaxis()->GetTitle(), 0., 1e-6, kFALSE);
  }
  
  if (canvasEff) {
    canvasEff->cd();

    TString hcopyName(Form("%s_copy", eff->GetName()));
    TGraphAsymmErrors* effCopy = static_cast<TGraphAsymmErrors*>(eff->Clone(hcopyName));
    FitGraphInPad(effCopy, canvasEff);
    effCopy->Draw("same P");

    Int_t icolor = -1;
    Int_t imarker = -1;

    if (zBin >= 0 && dptBin >= 0) {
      icolor = zBin;
      imarker = (dptBin-1)/2;
    }
    else if (zBin >= 0) {
      icolor = zBin;
    }
    else if (dptBin >= 0) {
      icolor = (dptBin-1)/2;
    }
  
    HistoStyler styler;
    styler.SetMarkerStyle(kFullCircle);
    styler.SetMarkerSize(0.8);
    styler.SetMarkerColor(kBlack);
    styler.SetVariableMarkerStyle();
    styler.SetVariableMarkerColor();
    styler.SetVariableLineColor();
    styler.SetLineWidth(1);
    styler.SetFillStyle(0);
    styler.Apply(effCopy, icolor, imarker);

    TLegend* leg = GetLegend(canvasEff);
    if (leg) {
      leg->AddEntry(effCopy, effCopy->GetTitle(), "pe");
    }
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseZMatrix(TCanvas*& canvasResp, TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t jetPtBin, Int_t dptBin)
{
  TString hname;

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  
  if (dptBin >= 0) {
    minDPt = params->GetDPtBin(dptBin);
    maxDPt = params->GetDPtBin(dptBin+1);
  }
  else {
    minDPt = params->GetMinDPt();
    maxDPt = params->GetMaxDPt();
  }

  Double_t minJetPt = 0;
  Double_t maxJetPt = 0;
  
  if (jetPtBin >= 0) {
    minJetPt = params->GetJetPtBin(jetPtBin);
    maxJetPt = params->GetJetPtBin(jetPtBin+1);
  }
  else {
    minJetPt = params->GetMinJetPt();
    maxJetPt = params->GetMaxJetPt();
  }

  hname = Form("%s_ResponseMatrix_Z_JetPt_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt), minDPt, maxDPt);
  TH2* resp = dynamic_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) {
    Printf("Error-DJetCorrResponse::PlotResponseZMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  if (!canvasResp && !GetCanvas(resp->GetName())) {
    canvasResp = SetUpCanvas(resp->GetName(),
                             resp->GetXaxis()->GetTitle(), 0., 1.0, kFALSE,
                             resp->GetYaxis()->GetTitle(), 0., 1.0, kFALSE);
  }

  if (canvasResp) {
    canvasResp->cd();
    resp->Draw("colz same");
    gPad->SetLogz();
  }
  
  hname = Form("%s_Efficiency_Z_JetPt_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt), minDPt, maxDPt);
  TGraphAsymmErrors* eff = dynamic_cast<TGraphAsymmErrors*>(fOutputList->FindObject(hname));
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseZMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }
  
  if (!canvasEff && !GetCanvas(eff->GetName())) {
    canvasEff = SetUpCanvas(eff->GetName(),
                            eff->GetXaxis()->GetTitle(), 0., 1.0, kFALSE,
                            eff->GetYaxis()->GetTitle(), 0., 1e-6, kFALSE);
  }

  if (canvasEff) {
    canvasEff->cd();

    TString hcopyName(Form("%s_copy", eff->GetName()));
    TGraphAsymmErrors* effCopy = static_cast<TGraphAsymmErrors*>(eff->Clone(hcopyName));
    FitGraphInPad(effCopy, canvasEff); 
    effCopy->Draw("same P");

    Int_t icolor = -1;
    Int_t imarker = -1;
  
    if (jetPtBin >= 0 && dptBin >= 0) {
      icolor = jetPtBin;
      imarker = (dptBin-1)/2;
    }
    else if (jetPtBin >= 0) {
      icolor = jetPtBin;
    }
    else if (dptBin >= 0) {
      icolor = (dptBin-1)/2;
    }
  
    HistoStyler styler;
    styler.SetMarkerStyle(kFullCircle);
    styler.SetMarkerSize(0.8);
    styler.SetMarkerColor(kBlack);
    styler.SetVariableMarkerStyle();
    styler.SetVariableMarkerColor();
    styler.SetVariableLineColor();
    styler.SetLineWidth(1);
    styler.SetFillStyle(0);
    styler.Apply(effCopy, icolor, imarker);

    TLegend* leg = GetLegend(canvasEff);
    if (leg) {
      leg->AddEntry(effCopy, effCopy->GetTitle(), "pe");
    }
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseDPtMatrix(TCanvas*& canvasResp, TCanvas*& canvasEff, DJetCorrAnalysisParams* params)
{
  TString hname;
  
  hname = Form("%s_ResponseMatrix_DpT", params->GetName());
  TH2* resp = dynamic_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) {
    Printf("Error-DJetCorrResponse::PlotResponseDPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }
  
  if (!canvasResp && !GetCanvas(resp->GetName())) {
    canvasResp = SetUpCanvas(resp->GetName(),
                             resp->GetXaxis()->GetTitle(), params->GetMinDPt(), params->GetMaxDPt(), kFALSE,
                             resp->GetYaxis()->GetTitle(), params->GetMinDPt(), params->GetMaxDPt(), kFALSE);
  }

  if (canvasResp) {
    canvasResp->cd();
    resp->Draw("colz same");
    gPad->SetLogz();
  }
  
  hname = Form("%s_Efficiency_DpT", params->GetName());
  TGraphAsymmErrors* eff = dynamic_cast<TGraphAsymmErrors*>(fOutputList->FindObject(hname));
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseDPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }
  
  if (!canvasEff && !GetCanvas(eff->GetName())) {
    canvasEff = SetUpCanvas(eff->GetName(),
                            eff->GetXaxis()->GetTitle(), params->GetMinDPt(), params->GetMaxDPt(), kFALSE,
                            eff->GetYaxis()->GetTitle(), 0., 1., kFALSE);
  }

  if (canvasEff) {
    canvasEff->cd();
    eff->Draw("same P");
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
TString DJetCorrResponse::GetTruthName(Int_t p)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";
  
  TString  hname = Form("%s_Truth_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), params->GetMinDPt(), params->GetMaxDPt());
  return hname;
}

//____________________________________________________________________________________
TString DJetCorrResponse::GetMeasuredName(Int_t p)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";
  
  TString  hname = Form("%s_Measured_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), params->GetMinDPt(), params->GetMaxDPt());
  return hname;
}

//____________________________________________________________________________________
TString DJetCorrResponse::GetResponseName(Int_t p)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";
  
  TString  hname = Form("%s_ResponseMatrix_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), params->GetMinDPt(), params->GetMaxDPt());
  return hname;
}

//____________________________________________________________________________________
TString DJetCorrResponse::GetMissesName(Int_t p)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";
  
  TString  hname = Form("%s_Misses_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), params->GetMinDPt(), params->GetMaxDPt());
  return hname;
}

//____________________________________________________________________________________
TString DJetCorrResponse::GetKinMissesName(Int_t p)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";
  
  TString  hname = Form("%s_KinMisses_JetPt_Z_DPt_%02.0f_%02.0f_Coarse", params->GetName(), params->GetMinDPt(), params->GetMaxDPt());
  return hname;
}

//____________________________________________________________________________________
TH2* DJetCorrResponse::GetTruth(Int_t p, Bool_t copy)
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
TH2* DJetCorrResponse::GetMeasured(Int_t p, Bool_t copy)
{
  TH2* hist = dynamic_cast<TH2*>(GetOutputHistogram(GetMeasuredName(p)));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH2*>(hist->Clone(hname));
  }

  return hist;
}

//____________________________________________________________________________________
THnSparse* DJetCorrResponse::GetResponse(Int_t p, Bool_t copy)
{
  THnSparse* hist = dynamic_cast<THnSparse*>(GetOutputSparseHistogram(GetResponseName(p)));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<THnSparse*>(hist->Clone(hname));
  }

  return hist;
}

//____________________________________________________________________________________
TH2* DJetCorrResponse::GetMisses(Int_t p, Bool_t copy)
{
  TH2* hist = dynamic_cast<TH2*>(GetOutputHistogram(GetMissesName(p)));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH2*>(hist->Clone(hname));
  }

  return hist;
}

//____________________________________________________________________________________
TH2* DJetCorrResponse::GetKinMisses(Int_t p, Bool_t copy)
{
  TH2* hist = dynamic_cast<TH2*>(GetOutputHistogram(GetKinMissesName(p)));

  if (copy && hist) {
    TString hname = hist->GetName();
    hname += "_copy";

    hist = static_cast<TH2*>(hist->Clone(hname));
  }

  return hist;
}
