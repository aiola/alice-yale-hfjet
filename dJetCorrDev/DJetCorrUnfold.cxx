// Class to unfold results from DJetCorrAnalysis using response matrix from DJetCorrResponse
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <vector>

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
#include <TSystem.h>
#include <TFile.h>
#include <TPaletteAxis.h>

#include "/opt/alicesw/RooUnfold/src/RooUnfold.h"
#include "/opt/alicesw/RooUnfold/src/RooUnfoldBayes.h"
#include "/opt/alicesw/RooUnfold/src/RooUnfoldResponse.h"

#include "MassFitter.h"
#include "HistoStyler.h"
#include "DJetCorrAnalysisParams.h"
#include "DJetCorrAnalysis.h"
#include "DJetCorrResponse.h"

#include "DJetCorrUnfold.h"

//TString DJetCorrUnfold::UnfoldingParams::fgkUnfoldingMethodNames[4] = { "NoUnfolding", "Bayesian", "SVD", "Chi2" };
const TString DJetCorrUnfold::fgkAxisLabels[2] = { "JetPt", "DmesonZ" };

ClassImp(DJetCorrUnfold);
ClassImp(DJetCorrUnfold::UnfoldingParams);

//____________________________________________________________________________________
DJetCorrUnfold::DJetCorrUnfold() :
  DJetCorrBase(),
  fDataParamIndex(0),
  fRespParamIndex(0),
  fForceRegeneration(kFALSE),
  fUseEfficiency(kTRUE),
  fUseKinEfficiency(kTRUE),
  fUnfoldingParams(UnfoldingParams::kBayesian),
  fAnalysis(0),
  fResponse(0),
  fRooUnfoldResponse(0)
{
  // Default constructor.

  fRooUnfoldResponse1D[kJetPtAxis] = 0;
  fRooUnfoldResponse1D[kDmesonZAxis] = 0;
}

//____________________________________________________________________________________
DJetCorrUnfold::DJetCorrUnfold(DJetCorrAnalysis* ana, DJetCorrResponse* resp) :
  DJetCorrBase(),
  fDataParamIndex(0),
  fRespParamIndex(0),
  fForceRegeneration(kFALSE),
  fUseEfficiency(kTRUE),
  fUseKinEfficiency(kTRUE),
  fUnfoldingParams(UnfoldingParams::kBayesian, 2, 10, 2),
  fAnalysis(ana),
  fResponse(resp),
  fRooUnfoldResponse(0)
{
  // Constructor.

  TString name(Form("%s_%s", ana->GetName(), resp->GetName()));
  SetName(name);
  SetTitle(name);

  fRooUnfoldResponse1D[kJetPtAxis] = 0;
  fRooUnfoldResponse1D[kDmesonZAxis] = 0;

  fUnfoldingParams1D[kJetPtAxis].Set(UnfoldingParams::kBayesian, 2, 10, 2);
  fUnfoldingParams1D[kDmesonZAxis].Set(UnfoldingParams::kBayesian, 2, 10, 2);
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
  
  TH2* truth = fAnalysis->GetTruth(fDataParamIndex, kTRUE);
  truth->SetName("histTruth");
  truth->SetTitle("Truth");
  fOutputList->Add(truth);

  TH2* measured = fAnalysis->GetMeasured(fDataParamIndex, kTRUE);
  measured->SetName("histMeasured");
  measured->SetTitle("Measured");
  fOutputList->Add(measured);

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
    Printf("Trying to load response matrix...");
    result = fResponse->LoadOutputHistograms();
  }
  
  if (result) {
    Printf("Response matrix successfully loaded!");
  }
  else {
    Printf("Generating response matrix...");
    result = fResponse->Regenerate();
    if (!result) {
      Printf("Failed to generate response matrix...");
      return kFALSE;
    }
  }

  TH2* responseTruth = fResponse->GetTruth(fRespParamIndex, kTRUE);
  if (!responseTruth) {
    Printf("Could not load response truth!");
    result = kFALSE;
  }
  responseTruth->SetName("histResponseTruth");
  responseTruth->SetTitle("Response Truth");
  fOutputList->Add(responseTruth);

  TH2* responseMeasured = fResponse->GetMeasured(fRespParamIndex, kTRUE);
  if (!responseMeasured) {
    Printf("Could not load response measured!");
    result = kFALSE;
  }
  responseMeasured->SetName("histResponseMeasured");
  responseMeasured->SetTitle("Response Measured");
  fOutputList->Add(responseMeasured);

  THnSparse* responseMatrix = fResponse->GetResponse(fRespParamIndex, kTRUE);
  if (!responseMatrix) {
    Printf("Could not load response matrix!");
    result = kFALSE;
  }
  responseMatrix->SetName("histResponseMatrix");
  responseMatrix->SetTitle("Response Matrix");
  fOutputList->Add(responseMatrix);

  TH2* responseMisses = fResponse->GetMisses(fRespParamIndex, kTRUE);
  if (!responseMisses) {
    Printf("Could not load response misses!");
    result = kFALSE;
  }
  responseMisses->SetName("histResponseMisses");
  responseMisses->SetTitle("Response Misses");
  fOutputList->Add(responseMisses);

  TH2* responseKinMisses = fResponse->GetKinMisses(fRespParamIndex, kTRUE);
  if (!responseKinMisses) {
    Printf("Could not load response kinematic misses!");
    result = kFALSE;
  }
  responseKinMisses->SetName("histResponseKinMisses");
  responseKinMisses->SetTitle("Response Kinematic Misses");
  fOutputList->Add(responseKinMisses);

  if (result) {
    Printf("Prepare response successful!");
  }

  return result;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::Start()
{
  Bool_t result = kFALSE;

  result = Init();
  if (!result) return kFALSE;

  result = PrepareData();
  if (!result) return kFALSE;

  result = PrepareResponse();
  if (!result) return kFALSE;

  result = MakeProjectionsBeforeUnfolding();
  if (!result) return kFALSE;

  result = GenerateRooUnfoldResponse();
  if (!result) return kFALSE;

  result = GenerateRooUnfoldResponse1D(kJetPtAxis);
  if (!result) return kFALSE;

  result = GenerateRooUnfoldResponse1D(kDmesonZAxis);
  if (!result) return kFALSE;

  result = Unfold();
  if (!result) return kFALSE;

  result = Unfold1D(kJetPtAxis);
  if (!result) return kFALSE;

  result = Unfold1D(kDmesonZAxis);
  if (!result) return kFALSE;

  result = MakeProjectionsAfterUnfolding();
  if (!result) return kFALSE;

  result = MakePlots();
  if (!result) return kFALSE;

  if (fSaveOutputFile) {
    result = SaveOutputFile();
    if (!result) return kFALSE;
  }

  return result;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::MakePlots()
{
  gStyle->SetOptStat(0);

  PlotResponse();
  PlotUnfolded();
  PlotRefolded();
  PlotFolded();

  PlotMeasured2D();
  PlotFolded2D();
  PlotTruth2D();

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::Unfold()
{
  if (!fRooUnfoldResponse) {
    Printf("ERROR - Response matrix was not correctly generated!");
    return kFALSE;
  }

  TH2* measured = GetMeasured();

  if (!measured) {
    Printf("ERROR - Measured histogram was not loaded!");
    return kFALSE;
  }
  
  TH2* truth = GetTruth();

  if (!truth) {
    Printf("ERROR - Truth histogram was not loaded!");
    return kFALSE;
  }

  TString hname;

  hname = GetFoldedName();
  Fold(truth, fRooUnfoldResponse, hname, hname);

  RooUnfoldBayes bayes(fRooUnfoldResponse, measured);

  Int_t regParam = 0;
  fUnfoldingParams.ResetCurrentRegIndex();
  while ((regParam = fUnfoldingParams.NextRegParam()) > 0) {
    TString name = Form("Unfold_RegParam%d", regParam);
    RooUnfold* unf = bayes.Clone(name);
    unf->SetRegParm(regParam);
    TH2* unfolded = static_cast<TH2*>(unf->Hreco());
    unfolded->GetXaxis()->SetName("xaxis");
    unfolded->GetYaxis()->SetName("yaxis");
    unfolded->GetZaxis()->SetName("zaxis");
    hname = GetUnfoldedName(regParam);
    unfolded->SetName(hname);
    unfolded->SetTitle(hname);

    fOutputList->Add(unfolded);

    hname = GetRefoldedName(regParam);
    Fold(unfolded, fRooUnfoldResponse, hname, hname);

    delete unf;
    unf = 0;
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::Unfold1D(EAxisType_t axis)
{
  if (!fRooUnfoldResponse1D[axis]) {
    Printf("ERROR - 1D Response matrix for axis %d was not correctly generated!", axis);
    return kFALSE;
  }

  TH1* measured = GetMeasuredProj(axis, 0);

  if (!measured) {
    Printf("ERROR - Measured histogram for axis %d was not loaded!", axis);
    return kFALSE;
  }

  switch (fUnfoldingParams1D[axis].GetUnfoldingMethod()) {
  case UnfoldingParams::kBayesian:
    Unfold1DBayes(measured, fRooUnfoldResponse1D[axis], fUnfoldingParams1D[axis]);
    break;
  default:
    Printf("Unfolding method %d not implemented.", fUnfoldingParams1D[axis].GetUnfoldingMethod());
    return kFALSE;
  }

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrUnfold::Unfold1DBayes(TH1* measured, RooUnfoldResponse* resp, UnfoldingParams params)
{
  RooUnfoldBayes bayes(resp, measured);

  TString hname;

  Int_t regParam = 0;
  params.ResetCurrentRegIndex();
  while ((regParam = params.NextRegParam()) > 0) {
    TString name = Form("Unfold_RegParam%d", regParam);
    RooUnfold* unf = bayes.Clone(name);
    unf->SetRegParm(regParam);
    TH2* unfolded = static_cast<TH2*>(unf->Hreco());
    unfolded->GetXaxis()->SetName("xaxis");
    unfolded->GetYaxis()->SetName("yaxis");
    unfolded->GetZaxis()->SetName("zaxis");
    hname = GetUnfoldedName(regParam);
    unfolded->SetName(hname);
    unfolded->SetTitle(hname);

    fOutputList->Add(unfolded);

    hname = GetRefoldedName(regParam);
    Fold(unfolded, resp, hname, hname);

    delete unf;
    unf = 0;
  }
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::SaveOutputFile()
{
  TString path(fAnalysis->GetOutpuFileName());
  path.Remove(path.Last('/'));

  path += "/Unfolding/";
  path += fResponse->GetName();

  if (gSystem->AccessPathName(path)) {
    Printf("Info-DJetCorrUnfold::SaveOutputFile : creating directory '%s'.", path.Data());
    gSystem->mkdir(path.Data(), kTRUE);
  }

  TString fname(path);
  fname += "/UnfoldingResults.root";

  TFile* outputFile = TFile::Open(fname, "recreate");

  if (!outputFile || outputFile->IsZombie()) {
    Printf("Error-DJetCorrUnfold::SaveOutputFile : could not open file '%s' to write.", fname.Data());
    outputFile = 0;
    return kFALSE;
  }

  outputFile->cd();

  Printf("Info-DJetCorrUnfold::SaveOutputFile : Now streaming results.");
  TIter next(fOutputList);
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

  Printf("Info-DJetCorrUnfold::SaveOutputFile : Closing the output file.");
  outputFile->Close();
  delete outputFile;
  outputFile = 0;

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::Fold(TH1* hist, RooUnfoldResponse* resp, TString name, TString title)
{
  if (!resp) {
    Printf("ERROR - Response matrix was not correctly generated!");
    return kFALSE;
  }

  TH1* histCopy = static_cast<TH2*>(hist->Clone("temp"));

  TH1* folded = static_cast<TH1*>(resp->ApplyToTruth(histCopy, name));
  folded->SetTitle(title);
  fOutputList->Add(folded);

  delete histCopy;

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::GenerateRooUnfoldResponse()
{
  // Preparing the response matrix.

  TH2* responseMeasured = GetResponseMeasured();
  TH2* responseTruth = GetResponseTruth();
  THnSparse* responseMatrix = GetResponseMatrix();
  TH2* responseMisses = GetResponseMisses();
  TH2* responseKinMisses = GetResponseKinMisses();

  if (!responseMeasured || !responseTruth || !responseMatrix) return kFALSE;
  
  if (fRooUnfoldResponse) delete fRooUnfoldResponse;

  fRooUnfoldResponse = new RooUnfoldResponse(responseMeasured, responseTruth);

  Int_t* coord_ind = new Int_t[responseMatrix->GetNdimensions()];
  Double_t* coord = new Double_t[responseMatrix->GetNdimensions()];
    
  for (Int_t ibin = 0; ibin < responseMatrix->GetNbins(); ibin++) {
    Double_t content = responseMatrix->GetBinContent(ibin, coord_ind);
    DJetCorrBase::GetBinCenter(responseMatrix, coord_ind, coord);
    fRooUnfoldResponse->Fill(coord[0], coord[1], coord[2], coord[3], content);
  }

  delete[] coord_ind;
  delete[] coord;

  if (fUseEfficiency) {
    AddEfficiency(fRooUnfoldResponse, responseMisses);
  }

  if (fUseKinEfficiency) {
    AddEfficiency(fRooUnfoldResponse, responseKinMisses);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::GenerateRooUnfoldResponse1D(EAxisType_t axis)
{
  // Preparing the response matrix.

  TH1* responseMeasured = GetResponseMeasuredProj(axis, 0);
  TH1* responseTruth = GetResponseTruthProj(axis, 0);
  TH2* responseMatrix = GetResponseMatrixProj(axis, 0);
  TH1* responseMisses = GetResponseMissesProj(axis, 0);
  TH1* responseKinMisses = GetResponseKinMissesProj(axis, 0);

  if (!responseMeasured) {
    Printf("Could not get response matrix measured projection for axis %d!", axis);
    return kFALSE;
  }

  if (!responseTruth) {
    Printf("Could not get response matrix truth projection for axis %d!", axis);
    return kFALSE;
  }

  if (!responseMatrix) {
    Printf("Could not get response matrix for axis %d!", axis);
    return kFALSE;
  }

  if (fRooUnfoldResponse1D[axis]) delete fRooUnfoldResponse1D[axis];

  fRooUnfoldResponse1D[axis] = new RooUnfoldResponse(responseMeasured, responseTruth);

  for (Int_t xbin = 0; xbin < responseMatrix->GetXaxis()->GetNbins(); xbin++) {
    Double_t x = responseMatrix->GetXaxis()->GetBinCenter(xbin);
    for (Int_t ybin = 0; ybin < responseMatrix->GetYaxis()->GetNbins(); ybin++) {
      Double_t y = responseMatrix->GetYaxis()->GetBinCenter(ybin);
      Double_t content = responseMatrix->GetBinContent(xbin, ybin);
      fRooUnfoldResponse1D[axis]->Fill(x, y, content);
    }
  }

  if (fUseEfficiency) {
    AddEfficiency1D(fRooUnfoldResponse1D[axis], responseMisses);
  }

  if (fUseKinEfficiency) {
    AddEfficiency1D(fRooUnfoldResponse1D[axis], responseKinMisses);
  }

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrUnfold::AddEfficiency1D(RooUnfoldResponse* resp, TH1* misses)
{
  if (!misses || !resp) return;

  for (Int_t xbin = 1; xbin <= misses->GetNbinsX(); xbin++) {
    Double_t x = misses->GetXaxis()->GetBinCenter(xbin);
    resp->Miss(x, misses->GetBinContent(xbin));
  }
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

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::MakeProjectionsBeforeUnfolding()
{
  std::vector<TString> hnames;
  hnames.push_back("histResponseTruth");
  hnames.push_back("histResponseMeasured");
  hnames.push_back("histResponseMisses");
  hnames.push_back("histResponseKinMisses");
  hnames.push_back("histTruth");
  hnames.push_back("histMeasured");

  for (std::vector<TString>::iterator it = hnames.begin() ; it != hnames.end(); ++it) {
    MakeProjections(*it);
  }

  MakeProjectionsTHnSparse(GetResponseMatrixName());

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::MakeProjectionsAfterUnfolding()
{
  std::vector<TString> hnames;
  hnames.push_back("histFolded");
  Int_t regParam = 0;
  fUnfoldingParams.ResetCurrentRegIndex();
  while ((regParam = fUnfoldingParams.NextRegParam()) > 0) {
    hnames.push_back(TString::Format("histUnfolded_RegParam%d", regParam));
    hnames.push_back(TString::Format("histRefolded_RegParam%d", regParam));
  }

  for (std::vector<TString>::iterator it = hnames.begin() ; it != hnames.end(); ++it) {
    MakeProjections(*it);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::MakeProjectionsTHnSparse(TString hname)
{
  THnSparse* hnsparse = static_cast<THnSparse*>(fOutputList->FindObject(hname));

  for (Int_t i = 0; i < 2; i++) {
    for (Int_t x = 0; x <= hnsparse->GetAxis((i+1)%2)->GetNbins(); x++) {
      TString pname;

      for (Int_t j = 0; j < 4; j++) {
        hnsparse->GetAxis(j)->SetRange(0, 0);
      }

      pname = TString::Format("%s_Proj%s_Bin%d", hnsparse->GetName(), fgkAxisLabels[i].Data(), x);
      if (x > 0)  hnsparse->GetAxis((i+1)%2)->SetRange(x, x);

      TH2* h2d = hnsparse->Projection(i+2, i);
      h2d->SetName(pname);
      h2d->SetTitle(pname);
      fOutputList->Add(h2d);
    }
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrUnfold::MakeProjections(TString hname)
{
  TH2* hist = static_cast<TH2*>(fOutputList->FindObject(hname));

  // X projection
  for (Int_t y = 0; y <= hist->GetNbinsY(); y++) {
    TString pname = TString::Format("%s_Proj%s_Bin%d", hist->GetName(), fgkAxisLabels[0].Data(), y);
    TString ptitle = TString::Format("Bin = %d", y);
    TH1* proj = 0;

    if (y > 0) {
      proj = hist->ProjectionX(pname, y, y);
    }
    else {
      proj = hist->ProjectionX(pname);
    }

    proj->SetTitle(ptitle);
    fOutputList->Add(proj);
  }

  // Y projection
  for (Int_t x = 0; x <= hist->GetNbinsX(); x++) {
    TString pname = TString::Format("%s_Proj%s_Bin%d", hname.Data(), fgkAxisLabels[1].Data(), x);
    TString ptitle = TString::Format("Bin = %d", x);

    TH1* proj = 0;

    if (x > 0) {
      proj = hist->ProjectionY(pname, x, x);
    }
    else {
      proj = hist->ProjectionY(pname);
    }

    proj->SetTitle(ptitle);
    fOutputList->Add(proj);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Int_t DJetCorrUnfold::GetNbinsX()
{
  TH2* truth = GetTruth();
  if (!truth) return 0;

  return truth->GetNbinsX();
}

//____________________________________________________________________________________
Int_t DJetCorrUnfold::GetNbinsY()
{
  TH2* truth = GetTruth();
  if (!truth) return 0;

  return truth->GetNbinsY();
}

//____________________________________________________________________________________
Int_t DJetCorrUnfold::GetNbins(Int_t axis)
{
  if (axis == kJetPtAxis) {
    return GetNbinsX();
  }
  else {
    return GetNbinsY();
  }
}

//____________________________________________________________________________________
void DJetCorrUnfold::PlotResponse()
{
  for (Int_t i = 0; i < 2; i++){
    TH2* resp = GetResponseMatrixProj((EAxisType_t)i, 2);
    if (!resp) return;

    TCanvas *canvas = SetUpCanvas(resp, kFALSE, kFALSE);
    canvas->SetLogz();
    resp->Draw("colz");

    canvas->Update();

    TPaletteAxis* palette = static_cast<TPaletteAxis*>(resp->GetListOfFunctions()->FindObject("palette"));
    if (palette) {
      palette->SetX2NDC(0.95);
    }

    if (fSavePlots) SavePlot(canvas);
  }
}

//____________________________________________________________________________________
void DJetCorrUnfold::PlotUnfolded(EAxisType_t axis)
{
  Int_t nbins = GetNbins((axis+1)%2);

  Printf("Axis %s, nbins = %d", fgkAxisLabels[axis].Data(), nbins);

  TObjArray histos;
  histos.SetOwner(kFALSE);

  TString title;

  for (Int_t x = 0; x <= nbins; x++) {
    histos.Clear();

    TH1* truth = GetTruthProj(axis, x);
    if (!truth) {
      Printf("Could not get truth for axis %s and bin %d!", fgkAxisLabels[axis].Data(), x);
      continue;
    }
    title = "Truth";
    truth->SetTitle(title);
    truth->SetOption("hist");
    truth->SetLineWidth(2);
    histos.Add(truth);

    TH1* prior = GetResponseTruthProj(axis, x);
    if (!prior) {
      Printf("Could not get prior for axis %s and bin %d!", fgkAxisLabels[axis].Data(), x);
      continue;
    }
    title = "Prior";
    prior->SetTitle(title);
    prior->SetOption("hist");
    prior->SetLineWidth(2);
    prior->SetLineStyle(2);
    prior->Scale(truth->Integral() / prior->Integral());
    histos.Add(prior);

    Int_t regParam = 0;
    fUnfoldingParams.ResetCurrentRegIndex();
    while ((regParam = fUnfoldingParams.NextRegParam()) > 0) {
      TH1* unfolded = GetUnfoldedProj(regParam, axis, x);
      if (!unfolded) {
        Printf("Could not get unfolded for axis %s, bin %d, reg param %d!", fgkAxisLabels[axis].Data(), x, regParam);
        continue;
      }
      title = TString::Format("Reg param = %d", regParam);
      unfolded->SetTitle(title);
      histos.Add(unfolded);
    }
    TString cname(Form("Unfolded_%s_Bin%d", fgkAxisLabels[axis].Data(), x));

    Plot1DHistos(cname, histos, 0, 0);
  }
}

//____________________________________________________________________________________
void DJetCorrUnfold::SavePlot(TCanvas* canvas)
{
  TString path(fAnalysis->GetOutpuFileName());
  path.Remove(path.Last('/'));

  path += "/Unfolding/";
  path += fResponse->GetName();
  path += "/plots";

  if (gSystem->AccessPathName(path)) {
    Printf("Info-DJetCorrUnfold::SavePlot : creating directory '%s'.", path.Data());
    gSystem->mkdir(path.Data(), kTRUE);
  }

  TString fname(path);
  fname += "/";
  fname += canvas->GetName();
  fname += ".";
  fname += fPlotFormat;

  canvas->SaveAs(fname);
}

//____________________________________________________________________________________
void DJetCorrUnfold::PlotUnfolded()
{
  PlotUnfolded(kJetPtAxis);
  PlotUnfolded(kDmesonZAxis);
}

//____________________________________________________________________________________
void DJetCorrUnfold::PlotRefolded()
{
  PlotRefolded(kJetPtAxis);
  PlotRefolded(kDmesonZAxis);
}

//____________________________________________________________________________________
void DJetCorrUnfold::PlotFolded()
{
  PlotFolded(kJetPtAxis);
  PlotFolded(kDmesonZAxis);
}

//____________________________________________________________________________________
void DJetCorrUnfold::PlotRefolded(EAxisType_t axis)
{
  Int_t nbins = GetNbins((axis+1)%2);

  Printf("Axis %s, nbins = %d", fgkAxisLabels[axis].Data(), nbins);

  TObjArray histos;
  histos.SetOwner(kFALSE);

  TString title;

  for (Int_t x = 0; x <= nbins; x++) {
    histos.Clear();

    TH1* meas = GetMeasuredProj(axis, x);
    if (!meas) {
      Printf("Could not get measured for axis %s and bin %d!", fgkAxisLabels[axis].Data(), x);
      continue;
    }
    title = "Measured";
    meas->SetTitle(title);
    histos.Add(meas);
    Int_t regParam = 0;
    fUnfoldingParams.ResetCurrentRegIndex();
    while ((regParam = fUnfoldingParams.NextRegParam()) > 0) {
      TH1* refolded = GetRefoldedProj(regParam, axis, x);
      if (!refolded) {
        Printf("Could not get refolded for axis %s, bin %d, reg param %d!", fgkAxisLabels[axis].Data(), x, regParam);
        continue;
      }
      title = TString::Format("Reg param = %d", regParam);
      refolded->SetTitle(title);
      histos.Add(refolded);
    }
    TString cname(Form("Refolded_%s_Bin%d", fgkAxisLabels[axis].Data(), x));

    Plot1DHistos(cname, histos, 0, 0);
  }
}

//____________________________________________________________________________________
void DJetCorrUnfold::PlotFolded(EAxisType_t axis)
{
  Int_t nbins = GetNbins((axis+1)%2);

  Printf("Axis %s, nbins = %d", fgkAxisLabels[axis].Data(), nbins);

  TObjArray histos;
  histos.SetOwner(kFALSE);

  TString title;

  for (Int_t x = 0; x <= nbins; x++) {
    histos.Clear();

    TH1* meas = GetMeasuredProj(axis, x);
    if (!meas) {
      Printf("Could not get measured for axis %s and bin %d!", fgkAxisLabels[axis].Data(), x);
      continue;
    }
    title = "Measured";
    meas->SetTitle(title);
    histos.Add(meas);
/*
    TH1* respMeas = GetResponseMeasuredProj(axis, x);
    if (!respMeas) {
      Printf("Could not get measured for axis %s and bin %d!", fgkAxisLabels[axis].Data(), x);
      continue;
    }
    title = "MC measured";
    respMeas->SetTitle(title);
    respMeas->Scale(meas->Integral() / respMeas->Integral());
    histos.Add(respMeas);
*/
    TH1* folded = GetFoldedProj(axis, x);
    if (!folded) {
      Printf("Could not get folded for axis %s and bin %d!", fgkAxisLabels[axis].Data(), x);
      continue;
    }
    title = "Folded";
    folded->SetTitle(title);
    histos.Add(folded);

    TString cname(Form("Folded_%s_Bin%d", fgkAxisLabels[axis].Data(), x));

    Plot1DHistos(cname, histos, 0, 0);
  }
}

//____________________________________________________________________________________
void DJetCorrUnfold::PlotMeasured2D()
{
  TH2* measured = GetMeasured();
  if (!measured) return;

  TCanvas *canvas = SetUpCanvas(measured, kFALSE, kFALSE);
  canvas->SetLogz();
  measured->Draw("colz");

  measured->SetMinimum(5e1);
  measured->SetMaximum(1e4);

  canvas->Update();

  TPaletteAxis* palette = static_cast<TPaletteAxis*>(measured->GetListOfFunctions()->FindObject("palette"));
  if (palette) {
    palette->SetX2NDC(0.95);
  }

  if (fSavePlots) SavePlot(canvas);
}

//____________________________________________________________________________________
void DJetCorrUnfold::PlotFolded2D()
{
  TH2* folded = GetFolded();
  if (!folded) return;

  TCanvas *canvas = SetUpCanvas(folded, kFALSE, kFALSE);
  canvas->SetLogz();
  folded->Draw("colz");

  folded->SetMinimum(5e1);
  folded->SetMaximum(1e4);

  canvas->Update();

  TPaletteAxis* palette = static_cast<TPaletteAxis*>(folded->GetListOfFunctions()->FindObject("palette"));
  if (palette) {
    palette->SetX2NDC(0.95);
  }

  if (fSavePlots) SavePlot(canvas);
}

//____________________________________________________________________________________
void DJetCorrUnfold::PlotTruth2D()
{
  TH2* truth = GetTruth();
  if (!truth) return;

  TCanvas *canvas = SetUpCanvas(truth, kFALSE, kFALSE);
  canvas->SetLogz();
  truth->Draw("colz");

  truth->SetMinimum(5e1);
  truth->SetMaximum(1e4);

  canvas->Update();

  TPaletteAxis* palette = static_cast<TPaletteAxis*>(truth->GetListOfFunctions()->FindObject("palette"));
  if (palette) {
    palette->SetX2NDC(0.95);
  }

  if (fSavePlots) SavePlot(canvas);
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetResponseTruthName()
{
  TString name;

  name = TString::Format("histResponseTruth");

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetResponseMeasuredName()
{
  TString name;

  name = TString::Format("histResponseMeasured");

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetResponseMatrixName()
{
  TString name;

  name = TString::Format("histResponseMatrix");

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetResponseMissesName()
{
  TString name;

  name = TString::Format("histResponseMisses");

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetResponseKinMissesName()
{
  TString name;

  name = TString::Format("histResponseKinMisses");

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetTruthName(Int_t)
{
  TString name;

  name = TString::Format("histTruth");

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetMeasuredName(Int_t)
{
  TString name;

  name = TString::Format("histMeasured");

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetUnfoldedName(Int_t regParam)
{
  TString name;

  name = TString::Format("histUnfolded_RegParam%d", regParam);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetRefoldedName(Int_t regParam)
{
  TString name;

  name = TString::Format("histRefolded_RegParam%d", regParam);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetFoldedName()
{
  TString name;

  name = TString::Format("histFolded");

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetTruthProjName(EAxisType_t axis, Int_t bin)
{
  TString name;

  name = TString::Format("histTruth_Proj%s_Bin%d", fgkAxisLabels[axis].Data(), bin);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetMeasuredProjName(EAxisType_t axis, Int_t bin)
{
  TString name;

  name = TString::Format("histMeasured_Proj%s_Bin%d", fgkAxisLabels[axis].Data(), bin);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetUnfoldedProjName(Int_t regParam, EAxisType_t axis, Int_t bin)
{
  TString name;

  name = TString::Format("histUnfolded_RegParam%d_Proj%s_Bin%d", regParam, fgkAxisLabels[axis].Data(), bin);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetRefoldedProjName(Int_t regParam, EAxisType_t axis, Int_t bin)
{
  TString name;

  name = TString::Format("histRefolded_RegParam%d_Proj%s_Bin%d", regParam, fgkAxisLabels[axis].Data(), bin);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetFoldedProjName(EAxisType_t axis, Int_t bin)
{
  TString name;

  name = TString::Format("histFolded_Proj%s_Bin%d", fgkAxisLabels[axis].Data(), bin);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetResponseTruthProjName(EAxisType_t axis, Int_t bin)
{
  TString name;

  name = TString::Format("histResponseTruth_Proj%s_Bin%d", fgkAxisLabels[axis].Data(), bin);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetResponseMeasuredProjName(EAxisType_t axis, Int_t bin)
{
  TString name;

  name = TString::Format("histResponseMeasured_Proj%s_Bin%d", fgkAxisLabels[axis].Data(), bin);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetResponseMatrixProjName(EAxisType_t axis, Int_t bin)
{
  TString name;

  name = TString::Format("histResponseMatrix_Proj%s_Bin%d", fgkAxisLabels[axis].Data(), bin);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetResponseMissesProjName(EAxisType_t axis, Int_t bin)
{
  TString name;

  name = TString::Format("histResponseMisses_Proj%s_Bin%d", fgkAxisLabels[axis].Data(), bin);

  return name;
}

//____________________________________________________________________________________
TString DJetCorrUnfold::GetResponseKinMissesProjName(EAxisType_t axis, Int_t bin)
{
  TString name;

  name = TString::Format("histResponseKinMisses_Proj%s_Bin%d", fgkAxisLabels[axis].Data(), bin);

  return name;
}

//____________________________________________________________________________________
DJetCorrUnfold::UnfoldingParams::UnfoldingParams() :
  TNamed(),
  fUnfoldingMethod(kNoUnfolding),
  fMinRegParam(0),
  fRegParamStep(0),
  fNRegParams(0),
  fCurrentRegIndex(0)
{
  fgkUnfoldingMethodNames[kNoUnfolding] = "NoUnfolding";
  fgkUnfoldingMethodNames[kBayesian] = "Bayesian";
  fgkUnfoldingMethodNames[kSVD] = "SVD";
  fgkUnfoldingMethodNames[kChi2] = "Chi2";

  SetName(fgkUnfoldingMethodNames[fUnfoldingMethod]);
  SetTitle(fgkUnfoldingMethodNames[fUnfoldingMethod]);
}

//____________________________________________________________________________________
DJetCorrUnfold::UnfoldingParams::UnfoldingParams(EUnfoldingMethod_t m, Int_t min, Int_t max, Int_t step) :
  TNamed(),
  fUnfoldingMethod(m),
  fMinRegParam(0),
  fRegParamStep(0),
  fNRegParams(0),
  fCurrentRegIndex(0)
{
  fgkUnfoldingMethodNames[kNoUnfolding] = "NoUnfolding";
  fgkUnfoldingMethodNames[kBayesian] = "Bayesian";
  fgkUnfoldingMethodNames[kSVD] = "SVD";
  fgkUnfoldingMethodNames[kChi2] = "Chi2";

  SetRegParam(min, max, step);
  SetName(fgkUnfoldingMethodNames[fUnfoldingMethod]);
  SetTitle(fgkUnfoldingMethodNames[fUnfoldingMethod]);
}


//____________________________________________________________________________________
DJetCorrUnfold::UnfoldingParams::UnfoldingParams(const UnfoldingParams &source) :
  TNamed(source),
  fUnfoldingMethod(source.fUnfoldingMethod),
  fMinRegParam(source.fMinRegParam),
  fRegParamStep(source.fRegParamStep),
  fNRegParams(source.fNRegParams),
  fCurrentRegIndex(source.fCurrentRegIndex)
{
  fgkUnfoldingMethodNames[kNoUnfolding] = "NoUnfolding";
  fgkUnfoldingMethodNames[kBayesian] = "Bayesian";
  fgkUnfoldingMethodNames[kSVD] = "SVD";
  fgkUnfoldingMethodNames[kChi2] = "Chi2";
}

//____________________________________________________________________________________
DJetCorrUnfold::UnfoldingParams& DJetCorrUnfold::UnfoldingParams::operator=(const UnfoldingParams& source)
{
  fUnfoldingMethod = source.fUnfoldingMethod;
  fMinRegParam = source.fMinRegParam;
  fRegParamStep = source.fRegParamStep;
  fNRegParams = source.fNRegParams;
  fCurrentRegIndex = source.fCurrentRegIndex;

  return *this;
}

//____________________________________________________________________________________
void DJetCorrUnfold::UnfoldingParams::Set(EUnfoldingMethod_t m, Int_t min, Int_t max, Int_t step)
{
  fUnfoldingMethod = m;
  SetRegParam(min, max, step);
  SetName(fgkUnfoldingMethodNames[fUnfoldingMethod]);
  SetTitle(fgkUnfoldingMethodNames[fUnfoldingMethod]);
}

//____________________________________________________________________________________
void DJetCorrUnfold::UnfoldingParams::SetRegParam(Int_t min, Int_t max, Int_t step)
{
  fMinRegParam = min;
  fRegParamStep = step;
  fNRegParams = TMath::CeilNint(Double_t(max - min) / step) + 1;
}

//____________________________________________________________________________________
Int_t DJetCorrUnfold::UnfoldingParams::NextRegParam()
{
  Int_t reg = RegParam(fCurrentRegIndex);
  if (reg >= 0) fCurrentRegIndex++;
  return reg;
}
