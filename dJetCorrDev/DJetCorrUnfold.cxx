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
  fTruth(0),
  fMeasured(0),
  fUnfolded(0),
  fRefolded(0),
  fResponseTruth(0),
  fResponseMeasured(0),
  fResponseMatrix(0),
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
  fTruth(0),
  fMeasured(0),
  fUnfolded(0),
  fRefolded(0),
  fResponseTruth(0),
  fResponseMeasured(0),
  fResponseMatrix(0),
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
  fAnalysis->SetInvMassPlotNorm(DJetCorrAnalysis::kNormalizeBackground);
  
  Bool_t result = kFALSE;

  if (!fForceRegeneration) {
    Printf("Now loading histograms"); 
    result = fAnalysis->LoadOutputHistograms();
  }
  
  if (!result) {
    result = fAnalysis->GenerateDJetCorrHistograms();
    if (!result) return kFALSE;
    
    result = fAnalysis->PlotDJetCorrHistograms(kTRUE);
    if (!result) return kFALSE;
  }

  Printf("Now getting histograms");
  
  fTruth = fAnalysis->GetTruth(fRespParamIndex, kTRUE);
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
    result = fResponse->ProjectResponseMatrices();
    if (!result) return kFALSE;
  }

  fResponseTruth = fResponse->GetTruth(fRespParamIndex, kTRUE);
  fResponseMeasured = fResponse->GetMeasured(fRespParamIndex, kTRUE);
  fResponseMatrix = fResponse->GetResponse(fRespParamIndex, kTRUE);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::Start()
{
  Bool_t result = kFALSE;

  Printf("now preparing data");
  result = PrepareData();
  if (!result) return kFALSE;

  Printf("now preparing response");
  result = PrepareResponse();
  if (!result) return kFALSE;

  Printf("now unfolding");
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
  fUnfolded = static_cast<TH2*>(bayes.Hreco());
  fUnfolded->GetXaxis()->SetName("xaxis");
  fUnfolded->GetYaxis()->SetName("yaxis");
  fUnfolded->GetZaxis()->SetName("zaxis");

  fUnfolded->Draw("colz");

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
  
  return resp;
}
