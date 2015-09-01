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
  fOutputFileName = "DJetCorrAnalysisComparer.root";
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
  fOutputFileName = "DJetCorrAnalysisComparer.root";

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
  TObjArray arrayDz(fAnalysisArray->GetEntriesFast());
  arrayDz.SetOwner(kTRUE);

  TObjArray arrayDPt(fAnalysisArray->GetEntriesFast());
  arrayDPt.SetOwner(kFALSE);

  TIter next(fAnalysisArray);

  DJetCorrBase* ana = 0;
  Int_t i = 0;

  TString compareDzName("DzSpectraComparisonMeasured");
  TString compareDPtName("DPtSpectraComparisonMeasured");

  while ((ana = static_cast<DJetCorrBase*>(next()))) {
    compareDzName += "_";
    compareDzName += ana->GetParamName(fParamIndexes[i]);

    compareDPtName += "_";
    compareDPtName += ana->GetParamName(fParamIndexes[i]);

    TH2* measuredDz = ana->GetDzMeasured(fParamIndexes[i], kTRUE);
    if (measuredDz) {
      measuredDz->Scale(1. / measuredDz->Integral(), "width");
      TString hname = Form("%s_%s", ana->GetName(), measuredDz->GetName());
      measuredDz->SetName(hname);
      measuredDz->SetTitle(ana->GetName());
      arrayDz.Add(measuredDz);
    }

    TH1* measuredDPt = ana->GetDPtMeasured(fParamIndexes[i], kTRUE);
    if (measuredDPt) {
      measuredDPt->Scale(1. / measuredDPt->Integral(), "width");
      TString hname = Form("%s_%s", ana->GetName(), measuredDPt->GetName());
      measuredDPt->SetName(hname);
      measuredDPt->SetTitle(ana->GetName());
      arrayDPt.Add(measuredDPt);
    }

    i++;
  }

  CompareZvsJetPt(compareDzName, arrayDz);
  CompareDPt(compareDPtName, arrayDPt);
}

//____________________________________________________________________________________
void DJetCorrAnalysisComparer::CompareTruth()
{
  TObjArray arrayDz(fAnalysisArray->GetEntriesFast());
  arrayDz.SetOwner(kTRUE);

  TObjArray arrayDPt(fAnalysisArray->GetEntriesFast());
  arrayDPt.SetOwner(kFALSE);

  TIter next(fAnalysisArray);

  DJetCorrBase* ana = 0;
  Int_t i = 0;

  TString compareDzName("DzSpectraComparisonTruth");
  TString compareDPtName("DPtSpectraComparisonTruth");

  while ((ana = static_cast<DJetCorrBase*>(next()))) {
    compareDzName += "_";
    compareDzName += ana->GetParamName(fParamIndexes[i]);

    compareDPtName += "_";
    compareDPtName += ana->GetParamName(fParamIndexes[i]);

    TH2* truthDz = ana->GetDzTruth(fParamIndexes[i], kTRUE);
    if (truthDz) {
      truthDz->Scale(1. / truthDz->Integral(), "width");
      TString hname = Form("%s_%s", ana->GetName(), truthDz->GetName());
      truthDz->SetName(hname);
      truthDz->SetTitle(ana->GetName());
      arrayDz.Add(truthDz);
    }

    TH1* truthDPt = ana->GetDPtTruth(fParamIndexes[i], kTRUE);
    if (truthDPt) {
      truthDPt->Scale(1. / truthDPt->Integral(), "width");
      TString hname = Form("%s_%s", ana->GetName(), truthDPt->GetName());
      truthDPt->SetName(hname);
      truthDPt->SetTitle(ana->GetName());
      arrayDPt.Add(truthDPt);
    }

    i++;
  }

  CompareZvsJetPt(compareDzName, arrayDz);
  CompareDPt(compareDPtName, arrayDPt);
}

//____________________________________________________________________________________
void DJetCorrAnalysisComparer::CompareZvsJetPt(const char* name, TObjArray& array)
{
  gStyle->SetOptStat(0);

  Printf("Entering method DJetCorrAnalysisComparer::CompareZvsJetPt: %d histograms to compare", array.GetEntriesFast());

  if (array.GetEntriesFast() == 0) return;

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
                           400*hist->GetNbinsX(), 400*nrows, nrows, hist->GetNbinsX());
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
        pad->SetLogy(kFALSE);
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

  if (fSavePlots) SavePlot(canvas);
}

//____________________________________________________________________________________
void DJetCorrAnalysisComparer::CompareDPt(const char* name, TObjArray& array)
{
  gStyle->SetOptStat(0);

  Printf("Entering method DJetCorrAnalysisComparer::CompareDPt: %d histograms to compare", array.GetEntriesFast());

  if (array.GetEntriesFast() == 0) return;

  TIter next(&array);

  TCanvas* canvas = 0;
  TH1* hist = 0;
  Int_t i = 0;

  TH1* baselineHisto = 0;

  HistoStyler styler;
  styler.SetMarkerStyle(kFullCircle);
  styler.SetMarkerSize(0.8);
  styler.SetVariableMarkerColor();
  styler.SetVariableLineColor();
  styler.SetLineWidth(1);

  TLegend* leg = SetUpLegend(0.40, 0.70, 0.76, 0.86, 12);

  while ((hist = static_cast<TH1*>(next()))) {
    if (!canvas) {
      Int_t nrows = 1;
      if (fMakeRatios) nrows = 2;
      canvas = SetUpCanvas(name,
                           "#it{p}_{T,D}", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax(), kFALSE,
                           "arb. units", 1, 1, kTRUE,
                           500, 500*nrows, nrows, 1);
    }

    TVirtualPad* pad = canvas->cd();
    if (fMakeRatios) pad = canvas->cd(1);

    styler.Apply(hist, i, 0);
    fOutputList->Add(hist);
    FitHistogramInPad(hist, pad);
    leg->AddEntry(hist, hist->GetTitle(), "pe");

    if (fMakeRatios) {
      pad = canvas->cd(2);
      pad->SetLogy(kFALSE);
      if (!baselineHisto) {
        baselineHisto = hist;

        TH1* blankHist = dynamic_cast<TH1*>(pad->GetListOfPrimitives()->At(0));
        if (blankHist) blankHist->GetYaxis()->SetTitle("ratio");
      }
      else {
        TString hname = Form("%sRatio", hist->GetName());
        TH1* ratio = static_cast<TH1*>(hist->Clone(hname));
        ratio->Divide(baselineHisto);
        styler.Apply(ratio, i, 0);
        FitHistogramInPad(ratio, pad);
        fOutputList->Add(ratio);
      }
    }
    i++;
  }

  canvas->cd(1);
  leg->Draw();

  if (fSavePlots) SavePlot(canvas);
}

