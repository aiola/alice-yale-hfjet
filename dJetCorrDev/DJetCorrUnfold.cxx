// Class to unfold results from DJetCorrAnalysis using response matrix from DJetCorrResponse
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TList.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THnSparse.h>
#include <TF1.h>
#include <Riostream.h>
#include <TObjString.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TLine.h>

#include "/opt/alicesw/RooUnfold/src/RooUnfold.h"
#include "/opt/alicesw/RooUnfold/src/RooUnfoldBayes.h"
#include "/opt/alicesw/RooUnfold/src/RooUnfoldResponse.h"

#include "MassFitter.h"
#include "HistoStyler.h"
#include "DJetCorrAnalysisParams.h"
#include "DJetCorrAnalysis.h"
#include "DJetCorrResponse.h"

#include "DJetCorrUnfold.h"

ClassImp(DJetCorrUnfold);

//____________________________________________________________________________________
DJetCorrUnfold::DJetCorrUnfold() :
  TNamed(),
  fDataParamIndex(0),
  fRespParamIndex(0),
  fForceRegeneration(kFALSE),
  fUseEfficiency(kTRUE),
  fUseKinEfficiency(kTRUE),
  fTruth(0),
  fMeasured(0),
  fUnfolded(0),
  fRefolded(0),
  fResponseTruth(0),
  fResponseMeasured(0),
  fResponseMatrix(0),
  fResponseMisses(0),
  fResponseKinMisses(0),
  fAnalysis(0),
  fResponse(0)
{
  // Default constructor.
}

//____________________________________________________________________________________
DJetCorrUnfold::DJetCorrUnfold(DJetCorrAnalysis* ana, DJetCorrResponse* resp) :
  TNamed(),
  fDataParamIndex(0),
  fRespParamIndex(0),
  fForceRegeneration(kFALSE),
  fUseEfficiency(kTRUE),
  fUseKinEfficiency(kTRUE),
  fTruth(0),
  fMeasured(0),
  fUnfolded(0),
  fRefolded(0),
  fResponseTruth(0),
  fResponseMeasured(0),
  fResponseMatrix(0),
  fResponseMisses(0),
  fResponseKinMisses(0),
  fAnalysis(ana),
  fResponse(resp)
{
  // Constructor.

  TString name(Form("%s_%s", ana->GetName(), resp->GetName()));
  SetName(name);
  SetTitle(name);
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::PrepareData()
{
  if (!fAnalysis) return kFALSE;
  
  fAnalysis->SetSavePlots(kFALSE);
  fAnalysis->SetAddTrainToCanvasName(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fForceRegeneration) {
    result = fAnalysis->LoadOutputHistograms();
  }
  
  if (!result) {
    result = fAnalysis->Regenerate();
    if (!result) return kFALSE;
  }
  
  fTruth = fAnalysis->GetTruth(fDataParamIndex, kTRUE);
  fMeasured = fAnalysis->GetMeasured(fDataParamIndex, kTRUE);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::PrepareResponse()
{
  if (!fResponse) return kFALSE;
  
  fResponse->SetSavePlots(kFALSE);
  fResponse->SetAddTrainToCanvasName(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fForceRegeneration) {
    result = fResponse->LoadOutputHistograms();
  }
  
  if (!result) {
    result = fResponse->Regenerate();
    if (!result) return kFALSE;
  }

  fResponseTruth = fResponse->GetTruth(fRespParamIndex, kTRUE);
  fResponseMeasured = fResponse->GetMeasured(fRespParamIndex, kTRUE);
  fResponseMatrix = fResponse->GetResponse(fRespParamIndex, kTRUE);
  fResponseMisses = fResponse->GetMisses(fRespParamIndex, kTRUE);
  fResponseKinMisses = fResponse->GetKinMisses(fRespParamIndex, kTRUE);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::Start()
{
  Bool_t result = kFALSE;

  result = PrepareData();
  if (!result) return kFALSE;

  result = PrepareResponse();
  if (!result) return kFALSE;

  result = Unfold();
  return result;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::Unfold()
{
  RooUnfoldResponse* resp = GenerateRooUnfoldResponse();

  if (!resp) {
    Printf("ERROR - Response matrix was not correctly generated!");
    return kFALSE;
  }

  if (!fMeasured) {
    Printf("ERROR - Measured histogram was not loaded!");
    return kFALSE;
  }
  
  RooUnfoldBayes bayes(resp, fMeasured);
  bayes.SetRegParm(2);
  fUnfolded = static_cast<TH2*>(bayes.Hreco());
  fUnfolded->GetXaxis()->SetName("xaxis");
  fUnfolded->GetYaxis()->SetName("yaxis");
  fUnfolded->GetZaxis()->SetName("zaxis");

  new TCanvas;
  fMeasured->Draw("colz");
  
  new TCanvas;
  fUnfolded->Draw("colz");

  new TCanvas;
  fTruth->Draw("colz");

  return kTRUE;
}

//____________________________________________________________________________________
RooUnfoldResponse* DJetCorrUnfold::GenerateRooUnfoldResponse()
{
  // Preparing the response matrix.

  if (!fResponseMeasured || !fResponseTruth || !fResponseMatrix) return 0;
  
  RooUnfoldResponse* resp = new RooUnfoldResponse(fResponseMeasured, fResponseTruth);

  Int_t* coord_ind = new Int_t[fResponseMatrix->GetNdimensions()];
  Double_t* coord = new Double_t[fResponseMatrix->GetNdimensions()];
    
  for (Int_t ibin = 0; ibin < fResponseMatrix->GetNbins(); ibin++) {
    Double_t content = fResponseMatrix->GetBinContent(ibin, coord_ind);
    DJetCorrBase::GetBinCenter(fResponseMatrix, coord_ind, coord);
    resp->Fill(coord[0], coord[1], coord[2], coord[3], content);
  }

  delete[] coord_ind;
  delete[] coord;

  if (fUseEfficiency) {
    AddEfficiency(resp, fResponseMisses);
  }

  if (fUseKinEfficiency) {
    AddEfficiency(resp, fResponseKinMisses);
  }
  
  return resp;
}

//____________________________________________________________________________________
void DJetCorrUnfold::AddEfficiency(RooUnfoldResponse* resp, TH2* misses)
{
  if (!misses || !resp) return;
  
  for (Int_t jbin = 1; jbin <= misses->GetNbinsX(); jbin++) {
    Double_t jetpt = misses->GetXaxis()->GetBinCenter(jbin);
    for (Int_t zbin = 1; zbin <= misses->GetNbinsY(); zbin++) {
      Double_t z = misses->GetYaxis()->GetBinCenter(zbin);

      resp->Miss(jetpt, z, misses->GetBinContent(jbin, zbin));
    }
  }
}
