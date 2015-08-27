// Helper class to compare results from DJetCorrAnalysis
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TList.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
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

#include "MassFitter.h"
#include "HistoStyler.h"
#include "DJetCorrAnalysisParams.h"
#include "DJetCorrAnalysis.h"

#include "DJetCorrAnalysisComparer.h"

ClassImp(DJetCorrAnalysisComparer);

//____________________________________________________________________________________
DJetCorrAnalysisComparer::DJetCorrAnalysisComparer() :
  DJetCorrBase(),
  fParamIndexes(),
  fForceRegeneration(kFALSE),
  fMakeRatios(kTRUE),
  fCompareTask(kNone),
  fAnalysisArray(0)
{
  // Default constructor.
}

//____________________________________________________________________________________
DJetCorrAnalysisComparer::DJetCorrAnalysisComparer(UInt_t task) :
  DJetCorrBase("DJetCorrAnalysisComparer", ""),
  fParamIndexes(),
  fForceRegeneration(kFALSE),
  fMakeRatios(kTRUE),
  fCompareTask(task),
  fAnalysisArray(0)
{
  // Default constructor.

  fAnalysisArray = new TObjArray();
}

//____________________________________________________________________________________
DJetCorrAnalysisComparer::DJetCorrAnalysisComparer(UInt_t task, DJetCorrBase* ana1, DJetCorrBase* ana2, Int_t ipar1, Int_t ipar2) :
  DJetCorrBase("DJetCorrAnalysisComparer", ""),
  fParamIndexes(10),
  fForceRegeneration(kFALSE),
  fMakeRatios(kTRUE),
  fCompareTask(task),
  fAnalysisArray(0)
{
  // Constructor for the very common case in which two analysis need to be compared.

  fAnalysisArray = new TObjArray();

  TString name = Form("%s_%s", ana1->GetName(), ana2->GetName());

  SetInputTrain(name);
  SetName(name);
  SetTitle(name);
  
  AddAnalysis(ana1, ipar1);
  AddAnalysis(ana2, ipar2);
}

//____________________________________________________________________________________
void DJetCorrAnalysisComparer::Start()
{
  if (!fAnalysisArray) fAnalysisArray = new TObjArray();

  Init();

  Prepare();
  ExecuteTasks();
}

//____________________________________________________________________________________
DJetCorrBase* DJetCorrAnalysisComparer::AddAnalysis(DJetCorrBase* ana, Int_t ipar)
{
  fAnalysisArray->Add(ana);
  if (fParamIndexes.GetSize() < fAnalysisArray->GetEntriesFast()) {
    fParamIndexes.Set(fAnalysisArray->GetEntriesFast()*2+1);
  }
  fParamIndexes[fAnalysisArray->GetEntriesFast()-1] = ipar;

  return ana;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysisComparer::Prepare()
{
  TIter next(fAnalysisArray);

  DJetCorrBase* ana = 0;

  while ((ana = static_cast<DJetCorrBase*>(next()))) {
    Bool_t result = kFALSE;

    if (!fForceRegeneration) {
      result = ana->LoadOutputHistograms();
    }

    if (!result) {
      ana->SetSavePlots(kFALSE);
      ana->SetAddTrainToCanvasName(kTRUE);

      result = ana->Regenerate();
      if (!result) return kFALSE;
    }
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysisComparer::ExecuteTasks()
{
  if ((fCompareTask & kCompareTruth) == kCompareTruth) CompareTruth();
  if ((fCompareTask & kCompareMeasured) == kCompareMeasured) CompareMeasured();

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrAnalysisComparer::CompareMeasured()
{
  TObjArray array(fAnalysisArray->GetEntriesFast());
  array.SetOwner(kTRUE);

  TIter next(fAnalysisArray);

  DJetCorrBase* ana = 0;
  Int_t i = 0;

  while ((ana = static_cast<DJetCorrBase*>(next()))) {
    TH2* measured = ana->GetMeasured(fParamIndexes[i], kTRUE);
    if (!measured) continue;
    measured->Scale(1. / measured->Integral());
    TString hname = Form("%s_%s", ana->GetName(), measured->GetName());
    measured->SetName(hname);
    measured->SetTitle(ana->GetName());
    array.Add(measured);
    i++;
  }

  CompareZvsJetPt("Measured", array);
}

//____________________________________________________________________________________
void DJetCorrAnalysisComparer::CompareTruth()
{
  TObjArray array(fAnalysisArray->GetEntriesFast());
  array.SetOwner(kTRUE);

  TIter next(fAnalysisArray);

  DJetCorrBase* ana = 0;
  Int_t i = 0;

  while ((ana = static_cast<DJetCorrBase*>(next()))) {
    TH2* truth = ana->GetTruth(fParamIndexes[i], kTRUE);
    if (!truth) continue;
    truth->Scale(1. / truth->Integral());
    TString hname = Form("%s_%s", ana->GetName(), truth->GetName());
    truth->SetName(hname);
    truth->SetTitle(ana->GetName());
    array.Add(truth);
    i++;
  }

  CompareZvsJetPt("Truth", array);
}

//____________________________________________________________________________________
void DJetCorrAnalysisComparer::CompareZvsJetPt(const char* name, TObjArray& array)
{
  gStyle->SetOptStat(0);

  Printf("Entering method DJetCorrAnalysisComparer::CompareZvsJetPt");

  TIter next(&array);

  TCanvas* canvas = 0;
  TH2* hist = 0;
  Int_t i = 0;

  TObjArray baselineHistos;

  HistoStyler styler;
  styler.SetMarkerStyle(kFullCircle);
  styler.SetMarkerSize(0.8);
  styler.SetVariableMarkerColor();
  styler.SetVariableLineColor();
  styler.SetLineWidth(1);

  TLegend* leg = SetUpLegend(0.10, 0.70, 0.46, 0.86, 12);

  while ((hist = static_cast<TH2*>(next()))) {

    if (!canvas) {
      Int_t nrows = 1;
      if (fMakeRatios) nrows = 2;
      canvas = SetUpCanvas(name,
                           "#it{z}_{||}", 0, 1., kFALSE,
                           "arb. units", 0, 0, kFALSE,
                           1200, 400*nrows, nrows, hist->GetNbinsX());
    }

    for (Int_t ibin = 1; ibin <= hist->GetNbinsX(); ibin++) {
      TString hname = Form("%s_JetPt_%d_%d", hist->GetName(), TMath::FloorNint(hist->GetXaxis()->GetBinLowEdge(ibin)*10), TMath::FloorNint(hist->GetXaxis()->GetBinUpEdge(ibin)*10));
      Printf("Now projecting %s", hname.Data());
      TH1* proj = hist->ProjectionY(hname, ibin, ibin);
      proj->SetTitle(hist->GetTitle());
      styler.Apply(proj, i, 0);
      fOutputList->Add(proj);
      FitHistogramInPad(proj, canvas->cd(ibin));
      if (ibin == 1) leg->AddEntry(proj, proj->GetTitle(), "pe");

      if (fMakeRatios) {
        TVirtualPad* pad = canvas->cd(ibin+hist->GetNbinsX());
        if (baselineHistos.GetEntriesFast() < hist->GetNbinsX()) {
          baselineHistos.Add(proj);

          TH1* blankHist = dynamic_cast<TH1*>(pad->GetListOfPrimitives()->At(0));
          if (blankHist) blankHist->GetYaxis()->SetTitle("ratio");
        }
        else {
          hname += "Ratio";
          TH1* ratio = static_cast<TH1*>(proj->Clone(hname));
          ratio->Divide(static_cast<TH1*>(baselineHistos.At(ibin-1)));
          styler.Apply(ratio, i, 0);
          FitHistogramInPad(ratio, pad);
          fOutputList->Add(ratio);
        }
      }
    }
    i++;
  }

  canvas->cd(1);
  leg->Draw();
}
