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
  TNamed(),
  fAnalysis1(0),
  fAnalysis2(0)
{
  // Default constructor.
}

//____________________________________________________________________________________
DJetCorrAnalysisComparer::DJetCorrAnalysisComparer(DJetCorrAnalysis* ana1, DJetCorrAnalysis* ana2) :
  TNamed(),
  fAnalysis1(ana1),
  fAnalysis2(ana2)
{
  // Constructor.

  TString name(Form("%s_%s", ana1->GetName(), ana2->GetName()));
  SetName(name);
  SetTitle(name);
}

//____________________________________________________________________________________
void DJetCorrAnalysisComparer::Start()
{
  fAnalysis1->SetSavePlots(kFALSE);
  fAnalysis1->SetAddTrainToCanvasName(kTRUE);
  fAnalysis2->SetSavePlots(kFALSE);
  fAnalysis2->SetAddTrainToCanvasName(kTRUE);
  
  fAnalysis1->PlotDJetCorrHistograms();
  fAnalysis2->PlotDJetCorrHistograms();
}
