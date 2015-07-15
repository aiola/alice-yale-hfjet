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

  fHistMatching->GetAxis(0)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(10)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());

  fHistJets2->GetAxis(2)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistJets2->GetAxis(4)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());

  hname = Form("ResponseMatrix_JetPt_Z");
  htitle = Form("Response matrix");
  Int_t dims[4] = {0, 9, 1, 10};
  THnSparse* resp = fHistMatching->Projection(4, dims, "O");

  resp->SetName(hname);
  resp->SetTitle(htitle);

  hname = Form("Efficiency_JetPt_Z");
  htitle = Form("Efficiency");
  TH2* eff = resp->Projection(3, 2);

  eff->SetName(hname);
  eff->SetTitle(htitle);
  eff->GetXaxis()->SetTitle("#it{p}_{jet}^{part} GeV/#it{c}");
  eff->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  eff->GetZaxis()->SetTitle("Efficiency");
  
  TH2* truth_proj = fHistJets2->Projection(4, 2, "O");

  eff->Divide(truth_proj);
  delete truth_proj;
  truth_proj = 0;

  fOutputList->Add(resp);
  fOutputList->Add(eff);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::ProjectResponseJetPtMatrix(DJetCorrAnalysisParams* params, Double_t minZ, Double_t maxZ)
{
  TString hname;
  TString htitle;

  if (!fHistMatching || !fHistJets2) return kFALSE;

  Int_t maxzBin = fHistMatching->GetAxis(10)->FindBin(maxZ);
  Int_t minzBin = fHistMatching->GetAxis(10)->FindBin(minZ);
      
  if (fHistMatching->GetAxis(10)->GetBinUpEdge(maxzBin) == 1.0) { // if the up edge == 1.0 includes the next bin
    maxzBin++;
    Printf("maxz == %.10f", fHistMatching->GetAxis(10)->GetBinUpEdge(maxzBin));
  }

  if (fHistMatching->GetAxis(10)->GetBinLowEdge(minzBin) == 1.0) { // if the low edge == 1.0 excludes the first bin
    minzBin++;
    Printf("minz == %.10f", fHistMatching->GetAxis(10)->GetBinLowEdge(minzBin));
  }
  
  fHistMatching->GetAxis(0)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistMatching->GetAxis(10)->SetRange(minzBin, maxzBin);

  fHistJets2->GetAxis(2)->SetRangeUser(params->GetMinJetPt(), params->GetMaxJetPt());
  fHistJets2->GetAxis(4)->SetRangeUser(minzBin, maxzBin);

  hname = Form("ResponseMatrix_JetPt_Z_%d_%d", TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100));
  htitle = Form("Response matrix for jet #it{p}_{T}: %.2f < #it{z}_{||}^{part} < %.2f", minZ, maxZ);
  TH2* resp = fHistMatching->Projection(1, 0, "O");

  resp->SetName(hname);
  resp->SetTitle(htitle);
  resp->GetXaxis()->SetTitle("#it{p}_{jet}^{det} GeV/#it{c}");
  resp->GetYaxis()->SetTitle("#it{p}_{jet}^{part} GeV/#it{c}");
  //resp->Rebin2D(5,5);

  hname = Form("Efficiency_JetPt_Z_%d_%d", TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100));
  htitle = Form("Efficiency for jet #it{p}_{T}: %.2f < #it{z}_{||}^{part} < %.2f", minZ, maxZ);
  TH1* eff = resp->ProjectionY(hname);
  
  eff->SetTitle(htitle);
  eff->GetXaxis()->SetTitle("#it{p}_{jet}^{part} GeV/#it{c}");
  eff->GetYaxis()->SetTitle("Efficiency");
  
  TH1* truth_proj = fHistJets2->Projection(2, "O");
  //truth_proj->Rebin(5);
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
  
  fHistMatching->GetAxis(1)->SetRangeUser(minJetPt, maxJetPt);
  fHistMatching->GetAxis(10)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());

  fHistJets2->GetAxis(2)->SetRangeUser(minJetPt, maxJetPt);
  fHistJets2->GetAxis(4)->SetRangeUser(params->GetMinZ(), params->GetMaxZ());

  hname = Form("ResponseMatrix_Z_JetPt_%d_%d", TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt));
  htitle = Form("Response matrix for #it{z}_{||}: %.1f < #it{p}_{T,jet}^{part} < %.1f GeV/#it{c}", minJetPt, maxJetPt);
  TH2* resp = fHistMatching->Projection(10, 9, "O");

  resp->SetName(hname);
  resp->SetTitle(htitle);
  resp->GetXaxis()->SetTitle("#it{z}_{||}^{det}");
  resp->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
  //resp->Rebin2D(5,5);

  hname = Form("Efficiency_Z_JetPt_%d_%d", TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt));
  htitle = Form("Efficiency for #it{z}_{||}: %.1f < #it{p}_{T,jet}^{part} < %.1f GeV/#it{c}", minJetPt, maxJetPt);
  TH1* eff = resp->ProjectionY(hname);
  
  eff->SetTitle(htitle);
  eff->GetXaxis()->SetTitle("#it{z}_{||}^{part}");
  eff->GetYaxis()->SetTitle("Efficiency");
  
  TH1* truth_proj = fHistJets2->Projection(4, "O");
  //truth_proj->Rebin(5);
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

  PlotResponseMatrix(params);
  PlotResponseMatricesVsJetPt(params);
  PlotResponseMatricesVsZ(params);
  
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
  for (Int_t i = 0; i < params->GetNzBins(); i++) {
    PlotResponseJetPtMatrix(params, params->GetzBin(i), params->GetzBin(i+1));
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseMatrix(DJetCorrAnalysisParams* /*params*/)
{
  TString hname;
  
  hname = Form("Efficiency_JetPt_Z");
  TH2* eff = static_cast<TH2*>(fOutputList->FindObject(hname));
  if (!eff) return kFALSE;

  SetUpCanvas(eff->GetName(),
              eff->GetXaxis()->GetTitle(), 0., 30., kFALSE,
              eff->GetYaxis()->GetTitle(), 0., 1.2, kFALSE);
  eff->Draw("colz same");

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseJetPtMatrix(DJetCorrAnalysisParams* /*params*/, Double_t minZ, Double_t maxZ)
{
  TString hname;

  hname = Form("ResponseMatrix_JetPt_Z_%d_%d", TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100));
  TH2* resp = static_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) return kFALSE;
  SetUpCanvas(resp, kFALSE, kFALSE);
  resp->Draw("colz same");
  
  hname = Form("Efficiency_JetPt_Z_%d_%d", TMath::CeilNint(minZ*100), TMath::CeilNint(maxZ*100));
  TH1* eff = static_cast<TH1*>(fOutputList->FindObject(hname));
  if (!eff) return kFALSE;
  SetUpCanvas(eff->GetName(),
              eff->GetXaxis()->GetTitle(), 0., eff->GetXaxis()->GetXmax(), kFALSE,
              eff->GetYaxis()->GetTitle(), 0., eff->GetMaximum()*1.5, kFALSE);
  eff->Draw("same");
 
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrResponse::PlotResponseZMatrix(DJetCorrAnalysisParams* /*params*/, Double_t minJetPt, Double_t maxJetPt)
{
  TString hname;

  hname = Form("ResponseMatrix_Z_JetPt_%d_%d", TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt));
  TH2* resp = static_cast<TH2*>(fOutputList->FindObject(hname));
  if (!resp) return kFALSE;
  SetUpCanvas(resp->GetName(),
              resp->GetXaxis()->GetTitle(), 0., 1.2, kFALSE,
              resp->GetYaxis()->GetTitle(), 0., 1.2, kFALSE);
  resp->Draw("colz same");
  
  hname = Form("Efficiency_Z_JetPt_%d_%d", TMath::CeilNint(minJetPt), TMath::CeilNint(maxJetPt));
  TH1* eff = static_cast<TH1*>(fOutputList->FindObject(hname));
  if (!eff) return kFALSE;
  SetUpCanvas(eff->GetName(),
              eff->GetXaxis()->GetTitle(), 0., 1.2, kFALSE,
              eff->GetYaxis()->GetTitle(), 0., eff->GetMaximum()*1.5, kFALSE);
  eff->Draw("same");

  return kTRUE;
}
