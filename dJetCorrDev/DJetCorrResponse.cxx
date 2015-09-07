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
  fOutputFileName = "DJetCorrResponse.root";
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
  
  ProjectResponseDPt(params);
  ProjectResponseJetPtZVsDPt(params);
  ProjectResponseZVsDPt(params);
  ProjectResponseZVsJetPt(params);
  ProjectResponseJetPtVsDPt(params);
  ProjectResponseJetPtVsZ(params);
  
  TH1::AddDirectory(addDirStatus);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseJetPtZVsDPt(DJetCorrAnalysisParams* params)
{
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt++) {
    ProjectResponseJetPtZ(params, idpt);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseZVsDPt(DJetCorrAnalysisParams* params)
{
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt++) {
    ProjectResponseZ(params, -1, idpt);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseZVsJetPt(DJetCorrAnalysisParams* params)
{
  for (Int_t ijetpt = 0; ijetpt < params->GetNJetPtBins(); ijetpt++) {
    ProjectResponseZ(params, ijetpt, -1);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseJetPtVsDPt(DJetCorrAnalysisParams* params)
{
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt++) {
    ProjectResponseJetPt(params, -1, idpt);
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseJetPtVsZ(DJetCorrAnalysisParams* params)
{
  for (Int_t iz = 0; iz < params->GetNzBins(); iz++) {
    ProjectResponseJetPt(params, iz, -1);
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseJetPtZ(DJetCorrAnalysisParams* params, Int_t dptBin)
{
  Printf("Entering method DJetCorrResponse::ProjectResponseJetPtZ");

  TString hname;
  TString htitle;

  if (!fHistMatching || !fHistJets2) return kFALSE;

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  params->GetDPtBinRange(minDPt, maxDPt, dptBin);

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
  measured->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{det} (GeV/#it{c})");
  measured->GetYaxis()->SetTitle("#it{z}_{||}^{det}");
  measured->GetZaxis()->SetTitle("Counts");

  hname = Form("%s_TruthReco_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);
  TH2* truthReco = fHistMatching->Projection(z2Axis, jetPt2Axis);
  truthReco->SetName(hname);
  truthReco->SetTitle("Truth reconstructed");
  truthReco->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{part} (GeV/#it{c})");
  truthReco->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  truthReco->GetZaxis()->SetTitle("Counts");

  hname = Form("%s_Truth_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);
  TH2* truth = fHistJets2->Projection(zPartAxis, jetPtPartAxis, "");
  truth->SetName(hname);
  truth->SetTitle("Truth");
  truth->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{part} (GeV/#it{c})");
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
  truthRecoAll->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{part} (GeV/#it{c})");
  truthRecoAll->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  truthRecoAll->GetZaxis()->SetTitle("Counts");

  hname = Form("%s_Efficiency_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);  
  TH2* eff = static_cast<TH2*>(truthReco->Clone(hname));
  eff->Divide(truth);
  eff->SetTitle("Efficiecny");
  eff->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{part} (GeV/#it{c})");
  eff->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  eff->GetZaxis()->SetTitle("Efficiency");

  hname = Form("%s_Misses_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);  
  TH2* misses = static_cast<TH2*>(truth->Clone(hname));
  misses->Add(truthReco, -1);
  misses->SetTitle("Misses = Truth - Truth reconstructed");
  misses->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{part} (GeV/#it{c})");
  misses->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  misses->GetZaxis()->SetTitle("Misses");

  hname = Form("%s_KinEfficiency_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);  
  TH2* kinEff = static_cast<TH2*>(truthReco->Clone(hname));
  kinEff->Divide(truthRecoAll);
  kinEff->SetTitle("Kinematic efficiency");
  kinEff->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{part} (GeV/#it{c})");
  kinEff->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  kinEff->GetZaxis()->SetTitle("Efficiency");

  hname = Form("%s_KinMisses_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);  
  TH2* kinMisses = static_cast<TH2*>(truthRecoAll->Clone(hname));
  kinMisses->Add(truthReco, -1);
  kinMisses->SetTitle("Misses = Truth reconstructed (all) - Truth reconstructed");
  kinMisses->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{part} (GeV/#it{c})");
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
Bool_t DJetCorrResponse::ProjectResponseDPt(DJetCorrAnalysisParams* params)
{
  Printf("Entering method DJetCorrResponse::ProjectResponseDPt");

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
  hname = Form("%s_ResponseMatrix_DpT_Full", params->GetName());
  htitle = Form("Response matrix for D meson #it{p}_{T}");
  TH2* respFull = fHistMatching->Projection(dPt2Axis, dPt1Axis, "A");
  respFull->SetName(hname);
  respFull->SetTitle(htitle);
  respFull->GetXaxis()->SetTitle("#it{p}_{T,D}^{det} (GeV/#it{c})");
  respFull->GetYaxis()->SetTitle("#it{p}_{T,D}^{part} (GeV/#it{c})");
  fOutputList->Add(respFull);

  TH1* totalFull = fHistJets2->Projection(dPtPartAxis, "A");
  totalFull->SetName("totalFull");
  totalFull->Sumw2();

  // Rebinning factor
  Int_t nbins = totalFull->GetXaxis()->FindBin(params->GetMaxDPt()) - totalFull->GetXaxis()->FindBin(params->GetMinDPt());
  if (nbins % 10 == 0) nbins /= 10;
  else if (nbins % 5 == 0) nbins /= 5;
  else if (nbins % 4 == 0) nbins /= 4;
  else if (nbins % 6 == 0) nbins /= 6;
  else if (nbins % 3 == 0) nbins /= 3;
  else if (nbins % 7 == 0) nbins /= 7;
  else if (nbins % 2 == 0) nbins /= 2;

  hname = Form("%s_ResponseMatrix_DpT", params->GetName());

  TH2* resp = Rebin(respFull, hname, nbins, params->GetMinDPt(), params->GetMaxDPt(), nbins, params->GetMinDPt(), params->GetMaxDPt());
  fOutputList->Add(resp);

  hname = Form("%s_ResponseMatrix_DpT_Norm", params->GetName());

  TH2* respNorm = Normalize(resp, hname);
  fOutputList->Add(respNorm);

  TH1* total = Rebin(totalFull, "total", nbins, params->GetMinDPt(), params->GetMaxDPt());

  TH1* pass = resp->ProjectionY("pass");
  //pass->Sumw2();

  hname = Form("%s_Efficiency_DpT", params->GetName());
  htitle = Form("Efficiency vs D meson #it{p}_{T}");
  if (fEfficiencyMode.Contains("TH1")) {
    TH1* eff = static_cast<TH1*>(pass->Clone(hname));
    eff->Divide(total);
    eff->SetTitle(htitle);
    eff->GetXaxis()->SetTitle("#it{p}_{T,D}^{part} (GeV/#it{c})");
    eff->GetYaxis()->SetTitle("Efficiency");

    fOutputList->Add(eff);
  }
  else {
    MakeBinomialConsistent(pass, total);
    TGraphAsymmErrors* eff = new TGraphAsymmErrors(pass, total, fEfficiencyMode);
    eff->SetName(hname);
    eff->SetTitle(htitle);
    eff->GetXaxis()->SetTitle("#it{p}_{T,D}^{part} (GeV/#it{c})");
    eff->GetYaxis()->SetTitle("Efficiency");

    fOutputList->Add(eff);
  }
  
  delete total;
  delete totalFull;
  delete pass;

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseJetPt(DJetCorrAnalysisParams* params, Int_t zBin, Int_t dptBin)
{
  Printf("Entering method DJetCorrResponse::ProjectResponseJetPt");

  TString hname;
  TString htitle;

  if (!fHistMatching || !fHistJets2) return kFALSE;

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  params->GetDPtBinRange(minDPt, maxDPt, dptBin);

  Double_t minZ = 0;
  Double_t maxZ = 0;
  params->GetzBinRange(minZ, maxZ, zBin);

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

  hname = Form("%s_ResponseMatrix_JetPt_Z_%d_%d_DPt_%02.0f_%02.0f_Full", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100), minDPt, maxDPt);
  htitle = Form("Response matrix for jet #it{p}_{T}: %.2f < #it{z}_{||}^{part} < %.2f", minZ, maxZ);
  TH2* respFull = fHistMatching->Projection(jetPt2Axis, jetPt1Axis, "A");
  respFull->SetName(hname);
  respFull->SetTitle(htitle);
  respFull->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{det} (GeV/#it{c})");
  respFull->GetYaxis()->SetTitle("#it{p}_{T,ch.jet}^{part} (GeV/#it{c})");
  fOutputList->Add(respFull);

  TH1* totalFull = fHistJets2->Projection(jetPtPartAxis, "A");
  totalFull->Sumw2();
  totalFull->SetName("totalFull");

  // Rebinning factor
  Int_t nbins = totalFull->GetXaxis()->FindBin(params->GetMaxJetPt()) - totalFull->GetXaxis()->FindBin(params->GetMinJetPt());
  if (nbins % 5 == 0) nbins /= 5;
  else if (nbins % 4 == 0) nbins /= 4;
  else if (nbins % 6 == 0) nbins /= 6;
  else if (nbins % 3 == 0) nbins /= 3;
  else if (nbins % 7 == 0) nbins /= 7;
  else if (nbins % 2 == 0) nbins /= 2;

  hname = Form("%s_ResponseMatrix_JetPt_Z_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100), minDPt, maxDPt);

  TH2* resp = Rebin(respFull, hname, nbins, params->GetMinJetPt(), params->GetMaxJetPt(), nbins, params->GetMinJetPt(), params->GetMaxJetPt());
  fOutputList->Add(resp);

  hname = Form("%s_ResponseMatrix_JetPt_Z_%d_%d_DPt_%02.0f_%02.0f_Norm", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100), minDPt, maxDPt);

  TH2* respNorm = Normalize(resp, hname);
  fOutputList->Add(respNorm);

  TH1* total = Rebin(totalFull, "total", nbins, params->GetMinJetPt(), params->GetMaxJetPt());
  TH1* pass = resp->ProjectionY("pass");
  //pass->Sumw2();

  hname = Form("%s_Efficiency_JetPt_Z_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100), minDPt, maxDPt);
  htitle = Form("%.1f < #it{p}_{T,D}^{part} < %.1f GeV/#it{c} and %.1f < #it{z}_{||}^{part} < %.1f", minDPt, maxDPt, minZ, maxZ);

  if (fEfficiencyMode.Contains("TH1")) {
    TH1* eff = static_cast<TH1*>(pass->Clone(hname));
    eff->Divide(total);
    eff->SetTitle(htitle);
    eff->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{part} (GeV/#it{c})");
    eff->GetYaxis()->SetTitle("Efficiency");

    fOutputList->Add(eff);
  }
  else {
    MakeBinomialConsistent(pass, total);
    TGraphAsymmErrors* eff = new TGraphAsymmErrors(pass, total, fEfficiencyMode);
    eff->SetName(hname);
    eff->SetTitle(htitle);
    eff->GetXaxis()->SetTitle("#it{p}_{T,ch.jet}^{part} (GeV/#it{c})");
    eff->GetYaxis()->SetTitle("Efficiency");
    fOutputList->Add(eff);
  }

  delete pass;
  delete total;
  delete totalFull;

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseZ(DJetCorrAnalysisParams* params, Int_t jetPtBin, Int_t dptBin)
{
  Printf("Entering method DJetCorrResponse::ProjectResponseZ with params=\"%s\", jetPtBin = %d, dptBin = %d",
         params->GetName(), jetPtBin, dptBin);

  TString hname;
  TString htitle;

  if (!fHistMatching || !fHistJets2) return kFALSE;

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  params->GetDPtBinRange(minDPt, maxDPt, dptBin);

  Double_t minJetPt = 0;
  Double_t maxJetPt = 0;
  params->GetJetPtBinRange(minJetPt, maxJetPt, jetPtBin);

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

  hname = Form("%s_ResponseMatrix_Z_JetPt_%d_%d_DPt_%02.0f_%02.0f_Full", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt), minDPt, maxDPt);
  htitle = Form("Response matrix for #it{z}_{||}: %.1f < #it{p}_{T,ch.jet}^{part} < %.1f GeV/#it{c}", minJetPt, maxJetPt);
  TH2* respFull = fHistMatching->Projection(z2Axis, z1Axis, "A");
  respFull->SetName(hname);
  respFull->SetTitle(htitle);
  respFull->GetXaxis()->SetTitle("#it{z}_{||}^{det}");
  respFull->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  fOutputList->Add(respFull);

  hname = Form("%s_ResponseMatrix_Z_JetPt_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt), minDPt, maxDPt);

  TH2* resp = Rebin(respFull, hname, params->GetNzBins(), params->GetzBins(), params->GetNzBins(), params->GetzBins());
  fOutputList->Add(resp);

  hname = Form("%s_ResponseMatrix_Z_JetPt_%d_%d_DPt_%02.0f_%02.0f_Norm", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt), minDPt, maxDPt);

  TH2* respNorm = Normalize(resp, hname);
  fOutputList->Add(respNorm);

  TH1* totalFull = fHistJets2->Projection(zPartAxis, "A");
  totalFull->SetName("totalFull");
  totalFull->Sumw2();

  TH1* total = Rebin(totalFull, "total", params->GetNzBins(), params->GetzBins());

  TH1* pass = resp->ProjectionY("pass");
  //pass->Sumw2();

  hname = Form("%s_Efficiency_Z_JetPt_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt), minDPt, maxDPt);
  htitle = Form("%.1f < #it{p}_{T,D}^{part} < %.1f GeV/#it{c} and %.1f < #it{p}_{T,ch.jet}^{part} < %.1f GeV/#it{c}", minDPt, maxDPt, minJetPt, maxJetPt);
  if (fEfficiencyMode.Contains("TH1")) {
    TH1* eff = static_cast<TH1*>(pass->Clone(hname));
    eff->Divide(total);
    eff->SetTitle(htitle);
    eff->GetXaxis()->SetTitle("#it{z}_{||}^{part}");
    eff->GetYaxis()->SetTitle("Efficiency");

    fOutputList->Add(eff);
  }
  else {
    MakeBinomialConsistent(pass, total);
    TGraphAsymmErrors* eff = new TGraphAsymmErrors(pass, total, fEfficiencyMode);
    eff->SetName(hname);
    eff->SetTitle(htitle);
    eff->GetXaxis()->SetTitle("#it{z}_{||}^{part}");
    eff->GetYaxis()->SetTitle("Efficiency");

    fOutputList->Add(eff);
  }

  delete pass;
  delete total;
  delete totalFull;

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

  PlotResponseJetPtZVsDPt(params);
  PlotResponseZVsDPt(params);
  PlotResponseZVsJetPt(params);
  PlotResponseJetPtVsDPt(params);
  PlotResponseJetPtVsZ(params);
  PlotResponseDPt(params);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseJetPtZVsDPt(DJetCorrAnalysisParams* params)
{
  TCanvas* canvasEff = 0;
  
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt+=1) {
    PlotEfficiencyJetPtZ(canvasEff, params, idpt);

    TPaveText* pave = SetUpPaveText(0.10, 0.86, 0.33, 0.975, 16, "ALICE Simulation, PYTHIA pp #sqrt{#it{s}} = 7 TeV");
    pave->AddText(params->GetLabel());
    pave->AddText(params->GetDPtRangeLabel(idpt));
    pave->Draw();

    if (fSavePlots) SavePlot(canvasEff);
    canvasEff = 0;
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseZVsDPt(DJetCorrAnalysisParams* params)
{
  TCanvas* canvasResp = 0;
  TCanvas* canvasEff = 0;
  
  TString cEffName = Form("%s_Efficiency_Z_Vs_DPt", params->GetName());
 
  canvasEff = SetUpCanvas(cEffName,
                          "#it{z}_{||}^{part}", params->GetMinZ(), params->GetMaxZ(), kFALSE,
                          "Efficiency", 0., 1e-6, kFALSE);
  canvasEff->cd();
  TLegend* leg = SetUpLegend(0.10, 0.83-0.05*params->GetNDPtBins(), 0.50, 0.88, 16);
  leg->Draw();
  
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt+=1) {
    Bool_t r = PlotResponseMatrixZ(canvasResp, params, -2, idpt);

    if (r) {
      canvasResp->cd();
      TPaveText* paveResp = SetUpPaveText(0.10, 0.86, 0.33, 0.975, 16, "ALICE Simulation, PYTHIA pp #sqrt{#it{s}} = 7 TeV");
      paveResp->AddText(params->GetLabel());
      TString respLabel = params->GetJetPtRangeLabel();
      respLabel += ", ";
      respLabel += params->GetDPtRangeLabel(idpt);
      paveResp->AddText(respLabel);
      paveResp->Draw();

      if (fSavePlots) SavePlot(canvasResp);
    }
    canvasResp = 0;

    PlotEfficiencyZ(canvasEff, params, -2, idpt);
  }

  canvasEff->cd();
  TPaveText* pave = SetUpPaveText(0.45, 0.85, 0.90, 0.58, 18, "ALICE Simulation");
  pave->AddText("PYTHIA pp #sqrt{#it{s}} = 7 TeV");
  pave->AddText(params->GetLabelDMeson());
  pave->AddText(params->GetLabelJet());
  pave->AddText(params->GetJetPtRangeLabel());
  pave->Draw();

  if (fSavePlots) SavePlot(canvasEff);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseZVsJetPt(DJetCorrAnalysisParams* params)
{
  TCanvas* canvasResp = 0;
  TCanvas* canvasEff = 0;

  TString cEffName = Form("%s_Efficiency_Z_Vs_JetPt", params->GetName());
 
  canvasEff = SetUpCanvas(cEffName,
                          "#it{z}_{||}^{part}", params->GetMinZ(), params->GetMaxZ(), kFALSE,
                          "Efficiency", 0., 1e-6, kFALSE);

  canvasEff->cd();
  TLegend* leg = SetUpLegend(0.10, 0.82-0.06*params->GetNJetPtBins(), 0.50, 0.86, 16);
  leg->Draw();

  for (Int_t ijetpt = -1; ijetpt < params->GetNJetPtBins(); ijetpt++) {
    Bool_t r = PlotResponseMatrixZ(canvasResp, params, ijetpt, -2);

    if (r) {
      canvasResp->cd();
      TPaveText* paveResp = SetUpPaveText(0.10, 0.86, 0.33, 0.975, 16, "ALICE Simulation, PYTHIA pp #sqrt{#it{s}} = 7 TeV");
      paveResp->AddText(params->GetLabel());
      TString respLabel = params->GetJetPtRangeLabel(ijetpt);
      respLabel += ", ";
      respLabel += params->GetDPtRangeLabel();
      paveResp->AddText(respLabel);
      paveResp->Draw();

      if (fSavePlots) SavePlot(canvasResp);
    }
    canvasResp = 0;

    PlotEfficiencyZ(canvasEff, params, ijetpt, -2);
  }

  canvasEff->cd();
  TPaveText* pave = SetUpPaveText(0.45, 0.85, 0.90, 0.58, 18, "ALICE Simulation");
  pave->AddText("PYTHIA pp #sqrt{#it{s}} = 7 TeV");
  pave->AddText(params->GetLabelDMeson());
  pave->AddText(params->GetLabelJet());
  pave->AddText(params->GetDPtRangeLabel());
  pave->Draw();

  if (fSavePlots) SavePlot(canvasEff);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseJetPtVsDPt(DJetCorrAnalysisParams* params)
{
  TCanvas* canvasResp = 0;
  TCanvas* canvasEff = 0;

  TString cEffName = Form("%s_Efficiency_JetPt_Vs_DPt", params->GetName());
 
  canvasEff = SetUpCanvas(cEffName,
                          "#it{p}_{T,ch.jet}^{ch}^{part} (GeV/#it{c})", params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                          "Efficiency", 0., 1e-6, kFALSE);

  canvasEff->cd();
  TLegend* leg = SetUpLegend(0.10, 0.83-0.05*params->GetNDPtBins(), 0.50, 0.88, 16);
  leg->Draw();
  
  for (Int_t idpt = -1; idpt < params->GetNDPtBins(); idpt+=1) {
    Bool_t r = PlotResponseMatrixJetPt(canvasResp, params, -2, idpt);

    if (r) {
      canvasResp->cd();
      TPaveText* paveResp = SetUpPaveText(0.10, 0.86, 0.33, 0.975, 16, "ALICE Simulation, PYTHIA pp #sqrt{#it{s}} = 7 TeV");
      paveResp->AddText(params->GetLabel());
      TString respLabel = params->GetzRangeLabel();
      respLabel += ", ";
      respLabel += params->GetDPtRangeLabel(idpt);
      paveResp->AddText(respLabel);
      paveResp->Draw();

      if (fSavePlots) SavePlot(canvasResp);
    }
    canvasResp = 0;

    PlotEfficiencyJetPt(canvasEff, params, -2, idpt);
  }

  canvasEff->cd();
  TPaveText* pave = SetUpPaveText(0.45, 0.85, 0.90, 0.58, 18, "ALICE Simulation");
  pave->AddText("PYTHIA pp #sqrt{#it{s}} = 7 TeV");
  pave->AddText(params->GetLabelDMeson());
  pave->AddText(params->GetLabelJet());
  pave->AddText(params->GetzRangeLabel());
  pave->Draw();

  if (fSavePlots) SavePlot(canvasEff);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseJetPtVsZ(DJetCorrAnalysisParams* params)
{

  TCanvas* canvasResp = 0;
  TCanvas* canvasEff = 0;

  TString cEffName = Form("%s_Efficiency_JetPt_Vs_Z", params->GetName());
 
  canvasEff = SetUpCanvas(cEffName,
                          "#it{p}_{T,ch.jet}^{ch}^{part} (GeV/#it{c})", params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                          "Efficiency", 0., 1e-6, kFALSE);

  canvasEff->cd();
  TLegend* leg = SetUpLegend(0.10, 0.83-0.06*params->GetNzBins(), 0.50, 0.88, 16);
  leg->Draw();

  for (Int_t iz = -1; iz < params->GetNzBins(); iz++) {
    Bool_t r = PlotResponseMatrixJetPt(canvasResp, params, iz, -2);

    if (r) {
      canvasResp->cd();
      TPaveText* paveResp = SetUpPaveText(0.10, 0.86, 0.33, 0.975, 16, "ALICE Simulation, PYTHIA pp #sqrt{#it{s}} = 7 TeV");
      paveResp->AddText(params->GetLabel());
      TString respLabel = params->GetzRangeLabel(iz);
      respLabel += ", ";
      respLabel += params->GetDPtRangeLabel();
      paveResp->AddText(respLabel);
      paveResp->Draw();

      if (fSavePlots) SavePlot(canvasResp);
    }
    canvasResp = 0;

    PlotEfficiencyJetPt(canvasEff, params, iz, -2);
  }

  canvasEff->cd();
  TPaveText* pave = SetUpPaveText(0.45, 0.85, 0.90, 0.58, 18, "ALICE Simulation");
  pave->AddText("PYTHIA pp #sqrt{#it{s}} = 7 TeV");
  pave->AddText(params->GetLabelDMeson());
  pave->AddText(params->GetLabelJet());
  pave->AddText(params->GetDPtRangeLabel());
  pave->Draw();

  if (fSavePlots) SavePlot(canvasEff);
 
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseDPt(DJetCorrAnalysisParams* params)
{
  TCanvas* canvasEff = 0;
  TCanvas* canvasResp = 0;

  Bool_t res1 = PlotResponseMatrixDPt(canvasResp, params);
  if (res1) {
    canvasResp->cd();
    TPaveText* paveResp = SetUpPaveText(0.10, 0.90, 0.33, 0.97, 16, "ALICE Simulation, PYTHIA pp #sqrt{#it{s}} = 7 TeV");
    paveResp->AddText(params->GetLabelDMeson());
    paveResp->Draw();

    if (fSavePlots) SavePlot(canvasResp);
  }

  Bool_t res2 = PlotEfficiencyDPt(canvasEff, params);
  if (res2) {
    canvasEff->cd();
    TPaveText* pave = SetUpPaveText(0.35, 0.85, 0.90, 0.66, 20, "ALICE Simulation");
    pave->AddText("PYTHIA pp #sqrt{#it{s}} = 7 TeV");
    pave->AddText(params->GetLabelDMeson());
    pave->Draw();

    if (fSavePlots) SavePlot(canvasEff);
  }

  return res1 && res2;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotEfficiencyJetPtZ(TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t dptBin)
{
  TString hname;

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  params->GetDPtBinRange(minDPt, maxDPt, dptBin);

  hname = Form("%s_Efficiency_JetPt_Z_DPt_%02.0f_%02.0f", params->GetName(), minDPt, maxDPt);
  TH2* eff = dynamic_cast<TH2*>(fOutputList->FindObject(hname));
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  Bool_t doPlotting = kFALSE;
  if (canvasEff) {
    doPlotting = kTRUE;
  }
  else {
    canvasEff = GetCanvas(eff->GetName());
    if (!canvasEff) {
      canvasEff = SetUpCanvas(hname,
                              eff->GetXaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                              eff->GetYaxis()->GetTitle(), params->GetMinZ(), params->GetMaxZ(), kFALSE,
                              700, 530, 1, 1,
                              0.10, 0.10, 0.10, 0.15);
      doPlotting = kTRUE;
    }
  }

  if (doPlotting) {
    canvasEff->cd();
    eff->Draw("colz same");
    eff->GetZaxis()->SetRangeUser(0.,1.);
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrixJetPt(TCanvas*& canvasResp, DJetCorrAnalysisParams* params, Int_t zBin, Int_t dptBin)
{
  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  params->GetDPtBinRange(minDPt, maxDPt, dptBin);

  Double_t minZ = 0;
  Double_t maxZ = 0;
  params->GetzBinRange(minZ, maxZ, zBin);

  TString hname = Form("%s_ResponseMatrix_JetPt_Z_%d_%d_DPt_%02.0f_%02.0f_Norm", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100), minDPt, maxDPt);
  TH2* resp = dynamic_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) {
    Printf("Error-DJetCorrResponse::PlotResponseJetPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  Bool_t doPlotting = kFALSE;
  if (canvasResp) {
    doPlotting = kTRUE;
  }
  else {
    canvasResp = GetCanvas(resp->GetName());
    if (!canvasResp) {
      canvasResp = SetUpCanvas(resp->GetName(),
                               resp->GetXaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                               resp->GetYaxis()->GetTitle(), params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                               700, 530, 1, 1,
                               0.10, 0.12, 0.10, 0.15);
      doPlotting = kTRUE;
    }
  }

  if (doPlotting) {
    canvasResp->cd();
    resp->Draw("colz same");
    gPad->SetLogz();
  }

  return doPlotting;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotEfficiencyJetPt(TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t zBin, Int_t dptBin)
{
  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  params->GetDPtBinRange(minDPt, maxDPt, dptBin);

  Double_t minZ = 0;
  Double_t maxZ = 0;
  params->GetzBinRange(minZ, maxZ, zBin);

  TString hname = Form("%s_Efficiency_JetPt_Z_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100), minDPt, maxDPt);
  TObject* eff = fOutputList->FindObject(hname);
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseJetPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  Bool_t doPlotting = kFALSE;
  if (canvasEff) {
    doPlotting = kTRUE;
  }
  else {
    canvasEff = GetCanvas(eff->GetName());
    if (!canvasEff) {
      canvasEff = SetUpCanvas(eff->GetName(),
                              "", params->GetMinJetPt(), params->GetMaxJetPt(), kFALSE,
                              "", 0., 1e-6, kFALSE);
      doPlotting = kTRUE;
    }
  }

  if (doPlotting) {
    TString hcopyName(Form("%s_copy", eff->GetName()));
    TString hcopyTitle;
    Int_t icolor = -1;
    Int_t imarker = -1;

    if (zBin >= -1 && dptBin >= -1) {
      icolor = zBin;
      imarker = dptBin;
      hcopyTitle = Form("%.1f < #it{p}_{T,D}^{part} < %.1f GeV/#it{c} and %.1f < #it{z}_{||}^{part} < %.1f",
                        minDPt, maxDPt, minZ, maxZ);
    }
    else if (zBin >= -1) {
      icolor = zBin;
      hcopyTitle = Form("%.1f < #it{z}_{||}^{part} < %.1f",
                        minZ, maxZ);
    }
    else if (dptBin >= -1) {
      icolor = dptBin;
      hcopyTitle = Form("%.0f < #it{p}_{T,D}^{part} < %.0f GeV/#it{c}",
                        minDPt, maxDPt);
    }

    TObject* effCopy = eff->Clone(hcopyName);

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

    FitObjectInPad(effCopy, canvasEff, "P E", kTRUE);

    TLegend* leg = GetLegend(canvasEff);
    if (leg) {
      leg->AddEntry(effCopy, hcopyTitle, "pe");
    }
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrixZ(TCanvas*& canvasResp, DJetCorrAnalysisParams* params, Int_t jetPtBin, Int_t dptBin)
{
  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  params->GetDPtBinRange(minDPt, maxDPt, dptBin);

  Double_t minJetPt = 0;
  Double_t maxJetPt = 0;
  params->GetJetPtBinRange(minJetPt, maxJetPt, jetPtBin);

  TString hname = Form("%s_ResponseMatrix_Z_JetPt_%d_%d_DPt_%02.0f_%02.0f_Norm",
                       params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt), minDPt, maxDPt);
  TH2* resp = dynamic_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) {
    Printf("Error-DJetCorrResponse::PlotResponseZMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  Bool_t doPlotting = kFALSE;
  if (canvasResp) {
    doPlotting = kTRUE;
  }
  else {
    canvasResp = GetCanvas(resp->GetName());
    if (!canvasResp) {
      canvasResp = SetUpCanvas(resp->GetName(),
                               resp->GetXaxis()->GetTitle(), params->GetMinZ(), params->GetMaxZ(), kFALSE,
                               resp->GetYaxis()->GetTitle(), params->GetMinZ(), params->GetMaxZ(), kFALSE,
                               700, 530, 1, 1,
                               0.10, 0.12, 0.10, 0.15);
      doPlotting = kTRUE;
    }
  }

  if (doPlotting) {
    canvasResp->cd();
    resp->Draw("colz same");
    gPad->SetLogz();
  }

  return doPlotting;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotEfficiencyZ(TCanvas*& canvasEff, DJetCorrAnalysisParams* params, Int_t jetPtBin, Int_t dptBin)
{
  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  params->GetDPtBinRange(minDPt, maxDPt, dptBin);

  Double_t minJetPt = 0;
  Double_t maxJetPt = 0;
  params->GetJetPtBinRange(minJetPt, maxJetPt, jetPtBin);

  TString hname = Form("%s_Efficiency_Z_JetPt_%d_%d_DPt_%02.0f_%02.0f", params->GetName(), TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt), minDPt, maxDPt);
  TObject* eff = fOutputList->FindObject(hname);
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseZMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  Bool_t doPlotting = kFALSE;
  if (canvasEff) {
    doPlotting = kTRUE;
  }
  else {
    canvasEff = GetCanvas(eff->GetName());
    if (!canvasEff) {
      canvasEff = SetUpCanvas(eff->GetName(),
                              "", 0., 1.0, kFALSE,
                              "", 0., 1e-6, kFALSE);
      doPlotting = kTRUE;
    }
  }

  if (doPlotting) {
    TString hcopyName(Form("%s_copy", eff->GetName()));
    TString hcopyTitle;
    Int_t icolor = -1;
    Int_t imarker = -1;

    if (jetPtBin >= -1 && dptBin >= -1) {
      icolor = jetPtBin;
      imarker = dptBin;
      hcopyTitle = Form("%.1f < #it{p}_{T,D}^{part} < %.1f GeV/#it{c} and %.1f < #it{p}_{T,ch.jet}^{part} < %.1f GeV/#it{c}",
                        minDPt, maxDPt, minJetPt, maxJetPt);
    }
    else if (jetPtBin >= -1) {
      icolor = jetPtBin;
      hcopyTitle = Form("%.0f < #it{p}_{T,ch.jet}^{part} < %.0f GeV/#it{c}",
                        minJetPt, maxJetPt);
    }
    else if (dptBin >= -1) {
      icolor = dptBin;
      hcopyTitle = Form("%.0f < #it{p}_{T,D}^{part} < %.0f GeV/#it{c}",
                        minDPt, maxDPt);
    }

    TObject* effCopy = eff->Clone(hcopyName);

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

    FitObjectInPad(effCopy, canvasEff, "P E", kTRUE, 1.9);

    TLegend* leg = GetLegend(canvasEff);
    if (leg) {
      leg->AddEntry(effCopy, hcopyTitle, "pe");
    }
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrixDPt(DJetCorrAnalysisParams* params)
{
  TCanvas *c = 0;

  return PlotResponseMatrixDPt(c, params);
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrixDPt(TCanvas*& canvasResp, DJetCorrAnalysisParams* params)
{
  TString hname = Form("%s_ResponseMatrix_DpT_Norm", params->GetName());
  TH2* resp = dynamic_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) {
    Printf("Error-DJetCorrResponse::PlotResponseDPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  Bool_t doPlotting = kFALSE;
  if (canvasResp) {
    doPlotting = kTRUE;
  }
  else {
    canvasResp = GetCanvas(resp->GetName());
    if (!canvasResp) {
      canvasResp = SetUpCanvas(resp->GetName(),
                               resp->GetXaxis()->GetTitle(), params->GetMinDPt(), params->GetMaxDPt(), kFALSE,
                               resp->GetYaxis()->GetTitle(), params->GetMinDPt(), params->GetMaxDPt(), kFALSE,
                               700, 500, 1, 1,
                               0.10, 0.12, 0.10, 0.10);
      doPlotting = kTRUE;
    }
  }

  if (doPlotting) {
    canvasResp->cd();
    resp->Draw("colz same");
    gPad->SetLogz();
  }

  return doPlotting;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotEfficiencyDPt(TCanvas*& canvasEff, DJetCorrAnalysisParams* params)
{
  TString hname = Form("%s_Efficiency_DpT", params->GetName());
  TObject* eff = fOutputList->FindObject(hname);
  if (!eff) {
    Printf("Error-DJetCorrResponse::PlotResponseDPtMatrix : Could not find histogram '%s'!", hname.Data());
    return kFALSE;
  }

  Bool_t doPlotting = kFALSE;
  if (canvasEff) {
    doPlotting = kTRUE;
  }
  else {
    canvasEff = GetCanvas(eff->GetName());
    if (!canvasEff) {
      canvasEff = SetUpCanvas(eff->GetName(),
                              "", params->GetMinDPt(), params->GetMaxDPt(), kFALSE,
                              "", 0., 0., kFALSE);
      doPlotting = kTRUE;
    }
  }

  if (doPlotting) {
    canvasEff->cd();
    FitObjectInPad(eff, canvasEff, "P E", kTRUE);
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::Regenerate()
{
  return ProjectResponseMatrices();
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
