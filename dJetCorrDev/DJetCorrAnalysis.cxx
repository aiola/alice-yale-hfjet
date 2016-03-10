// This class projects the THnSparse results of the AliAnalysisTaskDmesonJetCorrelations into TH1/TH2/TH3
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
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TLine.h>

#include "MassFitter.h"
#include "HistoStyler.h"

#include "DJetCorrAnalysisParams.h"

#include "DJetCorrAnalysis.h"

ClassImp(DJetCorrAnalysis);

//____________________________________________________________________________________
DJetCorrAnalysis::DJetCorrAnalysis() :
  DJetCorrBase(),
  fQAListName(),
  fInvMassPlotNorm(),
  fForceRefit(kFALSE),
  fDPtAxisTitle(),
  fDEtaAxisTitle(),
  fDPhiAxisTitle(),
  fDInvMassAxisTitle(),
  fD2prongInvMassAxisTitle(),
  fDDeltaInvMassAxisTitle(),
  fDSoftPionPtAxisTitle(),
  fDzAxisTitle(),
  fDeltaRAxisTitle(),
  fDeltaEtaAxisTitle(),
  fDeltaPhiAxisTitle(),
  fJetPtAxisTitle(),
  fJetEtaAxisTitle(),
  fJetPhiAxisTitle(),
  fJetLeadPtAxisTitle(),
  fJetAreaAxisTitle(),
  fJetConstAxisTitle(),
  fMatchingStatusAxisTitle(),
  fInputQAList(0),
  fMassFitters(0)
{
  // Default ctr.

  fAnaType = DJetCorrAnalysisParams::kInvMassAna;
}

//____________________________________________________________________________________
DJetCorrAnalysis::DJetCorrAnalysis(const char* train, const char* path) :
  DJetCorrBase(train, path),
  fQAListName("AliAnalysisTaskSAQA_AODFilterTracks_TPC_histos"),
  fInvMassPlotNorm(kPureCounts),
  fForceRefit(kFALSE),
  fDPtAxisTitle("#it{p}_{T,D} (GeV/#it{c})"),
  fDEtaAxisTitle("#eta_{D}"),
  fDPhiAxisTitle("#phi_{D} (rad)"),
  fDInvMassAxisTitle("#it{M}_{D} (GeV/#it{c}^{2})"),
  fD2prongInvMassAxisTitle("#it{M}_{2-prong} (GeV/#it{c}^{2})"),
  fDDeltaInvMassAxisTitle("#it{M}_{D*} - #it{M}_{D_{0}} (GeV/#it{c}^{2})"),
  fDSoftPionPtAxisTitle("#it{p}_{T,#pi} (GeV/#it{c})"),
  fDzAxisTitle("#it{z}_{D}"),
  fDeltaRAxisTitle("#Delta R_{D-jet}"),
  fDeltaEtaAxisTitle("#eta_{D} - #eta_{jet}"),
  fDeltaPhiAxisTitle("#phi_{D} - #phi_{jet} (rad)"),
  fJetPtAxisTitle("#it{p}_{T,jet} (GeV/#it{c})"),
  fJetEtaAxisTitle("#eta_{jet}"),
  fJetPhiAxisTitle("#phi_{jet} (rad)"),
  fJetLeadPtAxisTitle("#it{p}_{T,particle}^{leading} (GeV/#it{c})"),
  fJetAreaAxisTitle("#it{A}_{jet}"),
  fJetConstAxisTitle("No. of constituents"),
  fMatchingStatusAxisTitle("Matching status"),
  fInputQAList(0),
  fMassFitters(0)
{
  // Standard ctr.

  fAnaType = DJetCorrAnalysisParams::kInvMassAna;
  fInputFileName = "AnalysisResults.root";
  fInputDirFileName = "SA_DmesonJetCorr";
  fOutputFileName = "DJetCorrAnalysis.root";
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::Init()
{
  // Init class.

  Bool_t ret = DJetCorrBase::Init();

  if (!ret) return kFALSE;
  
  Printf("Info-DJetCorrAnalysis::Init : Now initializing.");
 
  if (fMassFitters) {
    delete fMassFitters;
    fMassFitters = 0;
  }
  
  fMassFitters = new TList();
  fMassFitters->SetOwner(kTRUE);
    
  Printf("Info-DJetCorrAnalysis::Init : Initialization done.");

  return kTRUE;
}

/*
//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::LoadTruthList(DJetCorrAnalysisParams* params)
{
  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;
  }
  
  result = params->LoadInputList(params->GetTruthInputListName());
  if (!result) return kFALSE;

  return kTRUE;
}
  

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectTruthSpectrum()
{
  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  DJetCorrAnalysisParams* params = 0;
  TIter next(fAnalysisParams);

  while ((params = static_cast<DJetCorrAnalysisParams*>(next()))) ProjectTruthSpectrum(params);
  
  TH1::AddDirectory(addDirStatus);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectTruthSpectrum(DJetCorrAnalysisParams* params)
{
  // Project histograms related to the D meson correlated to a jet.

  Bool_t result = kFALSE;

  result = LoadTruthList(params);
  if (!result) return kFALSE;
  
  result = LoadTHnSparse();
  if (!result) return kFALSE;
  
  TString prefix(params->GetName());

  Printf("Info-DJetCorrAnalysis::ProjectTruthSpectrum : Entering method.");

  ProjectTruthSpectrum(prefix, "AnyMatchingStatus_Truth", params, kAnyMatchingStatus, -1, -1, -1);
  ProjectTruthSpectrum(prefix, "NotMatched_Truth", params, kNotMatched, -1, -1, -1);
  ProjectTruthSpectrum(prefix, "Matched_Truth", params, kMatched, -1, -1, -1);

  for (Int_t i = 0; i < params->GetNDPtBins(); i++) {
    ProjectTruthSpectrum(prefix, "AnyMatchingStatus_Truth", params, kAnyMatchingStatus, i, -1, -1, 0);
    ProjectTruthSpectrum(prefix, "NotMatched_Truth", params, kNotMatched, i, -1, -1, 0);
    ProjectTruthSpectrum(prefix, "Matched_Truth", params, kMatched, i, -1, -1, 0);
  }

  for (Int_t i = 0; i < params->GetNJetPtBins(); i++) {
    ProjectTruthSpectrum(prefix, "Matched_Truth", params, kMatched, -1, i);
  }

  for (Int_t i = 0; i < params->GetNzBins(); i++) {
    ProjectTruthSpectrum(prefix, "Matched_Truth", params, kMatched, -1, -1, i);
    for (Int_t j = 0; j < params->GetNJetPtBins(); j++) {
      ProjectTruthSpectrum(prefix, "Matched_Truth", params, kMatched, -1, j, i);
    }
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectTruthSpectrum(TString prefix, TString suffix, DJetCorrAnalysisParams* params, EMatchingStatus st, Int_t dptBin, Int_t jetptBin, Int_t dzBin, Int_t minJetConst)
{
  TString cuts(params->GetCutString(st, dptBin, jetptBin, dzBin));

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  params->GetDPtBinRange(minDPt, maxDPt, dptBin);

  Double_t minz = 0;
  Double_t maxz = 0;
  params->GetzBinRange(minz, maxz, dzBin);

  Double_t minJetPt = 0;
  Double_t maxJetPt = 0;
  params->GetJetPtBinRange(minJetPt, maxJetPt, jetptBin);

  if (!suffix.IsNull()) suffix.Prepend("_");

  Int_t dPtAxis = GetAxisIndex(fDPtAxisTitle);
  if (dPtAxis < 0) return kFALSE;
  
  Int_t dEtaAxis = GetAxisIndex(fDEtaAxisTitle);
  if (dEtaAxis < 0) return kFALSE;
  
  Int_t dPhiAxis = GetAxisIndex(fDPhiAxisTitle);
  if (dPhiAxis < 0) return kFALSE;
  
  Int_t jetPtAxis = GetAxisIndex(fJetPtAxisTitle);
  if (jetPtAxis < 0) return kFALSE;
    
  Int_t jetEtaAxis = GetAxisIndex(fJetEtaAxisTitle);
  if (jetEtaAxis < 0) return kFALSE;
  
  Int_t jetPhiAxis = GetAxisIndex(fJetPhiAxisTitle);
  if (jetPhiAxis < 0) return kFALSE;

  Int_t dSoftPionPtAxis = GetAxisIndex(fDSoftPionPtAxisTitle);
  
  Int_t deltaRAxis = GetAxisIndex(fDeltaRAxisTitle);
  Int_t deltaEtaAxis = GetAxisIndex(fDeltaEtaAxisTitle);
  Int_t deltaPhiAxis = GetAxisIndex(fDeltaPhiAxisTitle);
  Int_t dzAxis = GetAxisIndex(fDzAxisTitle);
  Int_t jetConstAxis = GetAxisIndex(fJetConstAxisTitle);
  
  Int_t matchingStatusAxis = GetAxisIndex(fMatchingStatusAxisTitle);

  if (matchingStatusAxis >= 0) {
    if (st == kMatched) {
      fDmesons->GetAxis(matchingStatusAxis)->SetRange(1,2);
    }
    else if (st == kNotMatched) {
      fDmesons->GetAxis(matchingStatusAxis)->SetRange(3,4);
    }
    else {
      fDmesons->GetAxis(matchingStatusAxis)->SetRange(1,4);
    }
  }

  fDmesons->GetAxis(dPtAxis)->SetRangeUser(minDPt, maxDPt);

  if (st == kMatched) { // apply cuts on z and jet pt only if matched requested
    fDmesons->GetAxis(jetPtAxis)->SetRangeUser(minJetPt, maxJetPt);

    if (dzAxis >= 0) fDmesons->GetAxis(dzAxis)->SetRangeUser(minz, maxz);
    
    if (jetConstAxis >= 0) {
      if (minJetConst > 0) {
        fDmesons->GetAxis(jetConstAxis)->SetRangeUser(params->GetMinJetConstituents(), fDmesons->GetAxis(jetConstAxis)->GetXmax()+fgkEpsilon);
      }
      else if (minJetConst < 0) {
        fDmesons->GetAxis(jetConstAxis)->SetRangeUser(0, params->GetMinJetConstituents());
      }
      else {
        fDmesons->GetAxis(jetConstAxis)->SetRange(0, fDmesons->GetAxis(jetConstAxis)->GetNbins()+1);
      }
    }
  }
  else {
    fDmesons->GetAxis(dEtaAxis)->SetRangeUser(params->GetMinDEta(), params->GetMaxDEta()); // this is for non-matched D mesons: look only for D mesons in a defined eta region (usually the eta region is defined by the jet axis)
     
    if (jetPtAxis >= 0) fDmesons->GetAxis(jetPtAxis)->SetRange(0, fDmesons->GetAxis(jetPtAxis)->GetNbins()+1);
    if (dzAxis >= 0) fDmesons->GetAxis(dzAxis)->SetRange(0, fDmesons->GetAxis(dzAxis)->GetNbins()+1);
    if (jetConstAxis >= 0) fDmesons->GetAxis(jetConstAxis)->SetRange(0, fDmesons->GetAxis(jetConstAxis)->GetNbins()+1);
  }

  TH2* hDpos = fDmesons->Projection(dPhiAxis, dEtaAxis, "EO");
  hDpos->SetName(Form("h%s_MesonPhiVsEta_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
  hDpos->Rebin2D(2,2);
  fOutputList->Add(hDpos);

  TH1* hDEta = fDmesons->Projection(dEtaAxis, "EO");
  hDEta->SetName(Form("h%s_MesonEta_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
  hDEta->Rebin(2);
  fOutputList->Add(hDEta);

  if (dSoftPionPtAxis >= 0) {
    TH1* hdsoftpionpt = fDmesons->Projection(dSoftPionPtAxis, "EO");
    hdsoftpionpt->SetName(Form("h%s_SoftPion_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
    fOutputList->Add(hdsoftpionpt);
  }

  TString hname;
  
  if (dPtAxis >=0 && dptBin == -1) {
    hname = Form("h%s_MesonPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data());
    if (!fOutputList->Contains(hname)) {
      TH1* hDpt = fDmesons->Projection(dPtAxis, "EO");
      hDpt->SetName(hname);
      fOutputList->Add(hDpt);

      // Rebin to the coarse binning ready for unfolding
      TString hname = Form("%s_Coarse", hDpt->GetName());
      TH1* hDptCoarse = Rebin(hDpt, hname, params->GetNDPtBins(), params->GetDPtBins());
      fOutputList->Add(hDptCoarse);

      if (dSoftPionPtAxis >= 0) {
        TH2* hdsoftpionptVsDpt = fDmesons->Projection(dSoftPionPtAxis, dPtAxis, "EO");
        hdsoftpionptVsDpt->SetName(Form("h%s_SoftPionVsDPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hdsoftpionptVsDpt);
      }
    }
  }
  
  if (st == kMatched) {
    if (deltaRAxis >= 0) {
      TH1* hDeltaR = fDmesons->Projection(deltaRAxis, "EO");
      hDeltaR->SetName(Form("h%s_DeltaR_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDeltaR);

      if (dPtAxis >=0 && dptBin == -1) {
        TH2* hDeltaRVsDpt = fDmesons->Projection(deltaRAxis, dPtAxis, "EO");
        hDeltaRVsDpt->SetName(Form("h%s_DeltaRvsDPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hDeltaRVsDpt);
      }
    }

    if (jetConstAxis >= 0) {
      TH1* hJetConst = fDmesons->Projection(jetConstAxis, "EO");
      hJetConst->SetName(Form("h%s_JetConstituents_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hJetConst);

      if (dzAxis >= 0 && dzBin == -1) {
        TH2* hJetConstVsDz = fDmesons->Projection(jetConstAxis, dzAxis, "EO");
        hJetConstVsDz->SetName(Form("h%s_JetConstituentsVsDz_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hJetConstVsDz);
      }

      if (jetPtAxis >= 0 && jetptBin == -1) {
        TH2* hJetConstVsJetPt = fDmesons->Projection(jetConstAxis, jetPtAxis, "EO");
        hJetConstVsJetPt->SetName(Form("h%s_JetConstituentsVsJetPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hJetConstVsJetPt);
      }
    }

    if (deltaEtaAxis >= 0) {
      TH1* hDeltaEta = fDmesons->Projection(deltaEtaAxis, "EO");
      hDeltaEta->SetName(Form("h%s_DeltaEta_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDeltaEta);
    }

    if (deltaPhiAxis >= 0) {
      TH1* hDeltaPhi = fDmesons->Projection(deltaPhiAxis, "EO");
      hDeltaPhi->SetName(Form("h%s_DeltaPhi_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDeltaPhi);
    }
    
    if (dzAxis >= 0 && dzBin == -1) {
      TH1* hdz = fDmesons->Projection(dzAxis, "EO");
      hdz->SetName(Form("h%s_MesonZ_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      hdz->GetXaxis()->SetTitle("#it{z}_{||}^{part}");
      fOutputList->Add(hdz);

      // Rebin to the coarse binning ready for unfolding
      TString hname = Form("%s_Coarse", hdz->GetName());
      TH1* hdzCoarse = Rebin(hdz, hname, params->GetNzBins(), params->GetzBins());
      fOutputList->Add(hdzCoarse);

      if (dPtAxis >=0 && dptBin == -1) {
        TH2* hdzVsDPt = fDmesons->Projection(dzAxis, dPtAxis, "EO");
        hdzVsDPt->SetName(Form("h%s_MesonZvsDPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        hdzVsDPt->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
        fOutputList->Add(hdzVsDPt);
      }
      
      if (jetPtAxis >= 0 && jetptBin == -1) {
        TH2* hdzVsJetPt = fDmesons->Projection(dzAxis, jetPtAxis, "EO");
        hdzVsJetPt->SetName(Form("h%s_MesonZvsJetPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        hdzVsJetPt->GetXaxis()->SetTitle("#it{p}_{T,jet}^{part}");
        hdzVsJetPt->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
        fOutputList->Add(hdzVsJetPt);

        // Rebin to the coarse binning ready for unfolding
        TString hname = Form("%s_Coarse", hdzVsJetPt->GetName());
        TH2* hdzVsJetPtCoarse = Rebin(hdzVsJetPt, hname, params->GetNJetPtBins(), params->GetJetPtBins(), params->GetNzBins(), params->GetzBins());
        fOutputList->Add(hdzVsJetPtCoarse);
      }

      if (deltaRAxis >= 0) {
        TH2* hdzVsDeltaR = fDmesons->Projection(dzAxis, deltaRAxis, "EO");
        hdzVsDeltaR->SetName(Form("h%s_MesonZvsDeltaR_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        hdzVsDeltaR->GetYaxis()->SetTitle("#it{z}_{||}^{part}");
        fOutputList->Add(hdzVsDeltaR);
      }
    }

    if (jetPtAxis >= 0 && jetptBin == -1) {
      TH1* hJetPt = fDmesons->Projection(jetPtAxis, "EO");
      hJetPt->SetName(Form("h%s_JetPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      hJetPt->GetXaxis()->SetTitle("#it{p}_{T,jet}^{part}");
      fOutputList->Add(hJetPt);

      // Rebin to the coarse binning ready for unfolding
      TString hname = Form("%s_Coarse", hJetPt->GetName());
      TH1* hJetPtCoarse = Rebin(hJetPt, hname, params->GetNJetPtBins(), params->GetJetPtBins());
      fOutputList->Add(hJetPtCoarse);

      if (dPtAxis >=0 && dptBin == -1) {
        TH2* hDPtVsJetPt = fDmesons->Projection(dPtAxis, jetPtAxis, "EO");
        hDPtVsJetPt->SetName(Form("h%s_MesonPtvsJetPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hDPtVsJetPt);
      }
    }

    if (jetEtaAxis >= 0 && jetPhiAxis >= 0) {
      TH2* hJetpos = fDmesons->Projection(jetPhiAxis, jetEtaAxis, "EO");
      hJetpos->SetName(Form("h%s_JetPhiVsEta_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hJetpos);
    }
  }

  return kTRUE;
}
*/
//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateQAHistograms()
{
  // Generate QA histograms.

  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;
  }
  
  result = LoadQAList();
  if (!result) return kFALSE;
  
  result = ProjectQA();
  if (!result) return kFALSE;

  TH1::AddDirectory(addDirStatus);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateDJetCorrHistograms()
{
  // Generate D-jet correlation histograms.

  DJetCorrAnalysisParams* params = 0;
  TIter next(fAnalysisParams);

  while ((params = static_cast<DJetCorrAnalysisParams*>(next()))) GenerateDJetCorrHistograms(params);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateDJetCorrHistograms(Int_t i)
{
  // Generate D-jet correlation histograms.

  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(i));

  if (params) return GenerateDJetCorrHistograms(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateDJetCorrHistograms(const char* paramName)
{
  // Generate D-jet correlation histograms.

  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->FindObject(paramName));

  if (params) return GenerateDJetCorrHistograms(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateDJetCorrHistograms(DJetCorrAnalysisParams* params)
{
  // Generate D-jet correlation histograms.
  
  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;
  }
  
  result = OpenInputFile();
  if (!result) return kFALSE;

  result = params->LoadInputList(fInputFile);
  if (!result) return kFALSE;

  result = ProjectCorrD(params);
  if (!result) return kFALSE;

  TH1::AddDirectory(addDirStatus);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDJetCorrHistograms(Bool_t forceRefit)
{
  fForceRefit = forceRefit;
  
  DJetCorrAnalysisParams* params = 0;
  TIter next(fAnalysisParams);

  while ((params = static_cast<DJetCorrAnalysisParams*>(next()))) PlotDJetCorrHistograms(params);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDJetCorrHistograms(Int_t i)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(i));

  if (params) return PlotDJetCorrHistograms(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDJetCorrHistograms(const char* paramName)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->FindObject(paramName));

  if (params) return PlotDJetCorrHistograms(params);

  return kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDJetCorrHistograms(DJetCorrAnalysisParams* params)
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;

    result = LoadOutputHistograms();
    if (!result) result = GenerateDJetCorrHistograms();

    if (!result) return kFALSE;
  }

  if (fOutputList->GetEntries() == 0) return kFALSE;
  
  // Invariant mass
  Printf("Info-DJetCorrAnalysis::PlotDJetCorrHistograms : Plotting invariant mass vs D pt for all D mesons.");
  PlotInvMassHistogramsVsDPt(params, kAnyMatchingStatus);
  Printf("Info-DJetCorrAnalysis::PlotDJetCorrHistograms : Plotting invariant mass vs D pt for uncorrelated D mesons.");
  PlotInvMassHistogramsVsDPt(params, kNotMatched);
  Printf("Info-DJetCorrAnalysis::PlotDJetCorrHistograms : Plotting invariant mass vs D pt for correlated D mesons.");
  PlotInvMassHistogramsVsDPt(params, kMatched);
  Printf("Info-DJetCorrAnalysis::PlotDJetCorrHistograms : Plotting invariant mass vs z for correlated D mesons.");
  PlotInvMassHistogramsVsDz(params);
  for (Int_t j = 0; j < params->GetNJetPtBins(); j++) {
    PlotInvMassHistogramsVsDz(params, -1, j);
  }

  GenerateMeaduredSpectrum(params);

  // Delta R
  PlotObservable(params, "DeltaR", 0, 0,
                 -1, -1,
                 params->GetMinJetPt(), params->GetMaxJetPt(),
                 params->GetMinZ(), params->GetMaxZ(),
                 2, 4, 1);

  // Delta R
  PlotObservable(params, "DeltaR", 0, 0,
                 params->GetMinDPt(), params->GetMaxDPt(), 
                 -1, -1,
                 params->GetMinZ(), params->GetMaxZ(),
                 1, 1, 1);

  // Delta Phi
  PlotObservable(params, "DeltaPhi", -0.7, 0.7,
                 params->GetMinDPt(), params->GetMaxDPt(), 
                 -1, -1,
                 params->GetMinZ(), params->GetMaxZ(),
                 1, 1, 1);

  // Delta Eta
  PlotObservable(params, "DeltaEta", -0.7, 0.7,
                 params->GetMinDPt(), params->GetMaxDPt(), 
                 -1, -1,
                 params->GetMinZ(), params->GetMaxZ(),
                 1, 1, 1);

  PlotDPtSpectraVsJetPt(params, kTRUE);
  PlotDPtSpectraVsDz(params, kTRUE);
  PlotDPtSpectraVsMatchingStatus(params, kTRUE);
  PlotDzSpectraVsJetPt(params, kTRUE);
  
  TH1::AddDirectory(addDirStatus);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::LoadOutputHistograms()
{
  TFile* outputFile = OpenOutputFile();

  if (!outputFile) return kFALSE;

  fOutputList = static_cast<THashList*>(outputFile->Get("fOutputList"));
  fMassFitters = static_cast<TList*>(outputFile->Get("fMassFitters"));

  outputFile->Close();
  delete outputFile;
  outputFile = 0;

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::SaveOutputFile()
{
  // Save the output in a file.

  TObjArray arr;
  arr.SetOwner(kFALSE);
  
  fOutputList->SetName("fOutputList");
  fMassFitters->SetName("fMassFitters");
  
  arr.Add(fOutputList);
  arr.Add(fMassFitters);

  return DJetCorrBase::SaveOutputFile(arr);
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotTrackHistograms()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;

    result = LoadOutputHistograms();
    if (!result) result = GenerateQAHistograms();
    
    if (!result) return kFALSE;
  }

  if (fOutputList->GetEntries() == 0) return kFALSE;

  TString labels[4] = {"Hybrid", "Global", "Constr", "ConstrNoITS"};
  Color_t colors[4] = {kBlack, kBlue+2, kRed+2, kGreen+2};
  TString titles[4] = {"Hybrid tracks", "Global tracks", "Constrained, w/ ITS refit", "Constrained, w/o ITS refit"};

  TCanvas* cTrackPt = 0;
  TCanvas* cTrackEta = 0;
  TCanvas* cTrackPhi = 0;

  TLegend* legPt = SetUpLegend(0.58, 0.68, 0.88, 0.88, 14);
  TLegend* legEta = SetUpLegend(0.19, 0.30, 0.49, 0.50, 14);
  TLegend* legPhi = SetUpLegend(0.58, 0.68, 0.88, 0.88, 14);
  
  for (Int_t i = 0; i < 4; i++) {
    TH1* hTrEta = dynamic_cast<TH1*>(fOutputList->FindObject(Form("hTracks_%s_Eta", labels[i].Data())));
    TH1* hTrPhi = dynamic_cast<TH1*>(fOutputList->FindObject(Form("hTracks_%s_Phi", labels[i].Data())));
    TH1* hTrPt = dynamic_cast<TH1*>(fOutputList->FindObject(Form("hTracks_%s_Pt", labels[i].Data())));

    if (!hTrEta) continue;
    if (!hTrPhi) continue;
    if (!hTrPt) continue;
    
    if (!cTrackPt) cTrackPt = SetUpCanvas("fig_TrackPt", "#it{p}_{T} (GeV/#it{c})", 0, 50, kFALSE, "counts", 1e-1, hTrPt->GetMaximum()*3, kTRUE);
    if (!cTrackEta) cTrackEta = SetUpCanvas("fig_TrackEta", "#eta", -1, 1, kFALSE, "counts", 0, hTrEta->GetMaximum()*1.2, kFALSE);
    if (!cTrackPhi) cTrackPhi = SetUpCanvas("fig_TrackPhi", "#phi", 0, TMath::TwoPi(), kFALSE, "counts", 0, hTrPhi->GetMaximum()*1.5, kFALSE);  

    cTrackPt->cd();
    TH1* hTrPt_copy = hTrPt->DrawCopy("same hist");
    hTrPt_copy->SetLineColor(colors[i]);
    hTrPt_copy->SetTitle(titles[i]);

    cTrackEta->cd();
    TH1* hTrEta_copy = hTrEta->DrawCopy("same hist");
    hTrEta_copy->SetLineColor(colors[i]);
    hTrEta_copy->SetTitle(titles[i]);

    cTrackPhi->cd();
    TH1* hTrPhi_copy = hTrPhi->DrawCopy("same hist");
    hTrPhi_copy->SetLineColor(colors[i]);
    hTrPhi_copy->SetTitle(titles[i]);

    legPt->AddEntry(hTrPt_copy, titles[i], "l");
    legEta->AddEntry(hTrEta_copy, titles[i], "l");
    legPhi->AddEntry(hTrPhi_copy, titles[i], "l");
  }

  if (cTrackPt) {
    cTrackPt->cd();
    legPt->Draw();
    cTrackPt->Update();
  }

  if (cTrackEta) {
    cTrackEta->cd();
    legEta->Draw();
    cTrackEta->Update();
  }

  if (cTrackPhi) {
    cTrackPhi->cd();
    legPhi->Draw();
    cTrackPhi->Update();
  }
  
  if (fSavePlots) {
    SavePlot(cTrackPt);
    SavePlot(cTrackEta);
    SavePlot(cTrackPhi);
  }

  TH1::AddDirectory(addDirStatus);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDPtSpectraVsJetPt(DJetCorrAnalysisParams* params, Bool_t eventScaling)
{
  TH1** histos = new TH1*[params->GetNJetPtBins()];
  
  for (Int_t i = 0; i < params->GetNJetPtBins(); i++) {
    histos[i] = 0;
    
    TString spectrumCuts(params->GetCutString(kMatched, -1, i, -1));
    TString spectrumName(Form("h%s_DPtSpectrum_%s_Matched", params->GetName(), spectrumCuts.Data()));
    TH1* hs = dynamic_cast<TH1*>(fOutputList->FindObject(spectrumName));
    if (!hs) {
      Printf("Error-DJetCorrAnalysis::PlotDPtSpectraVsJetPt : Histogram '%s' not found!", spectrumName.Data());
      continue;
    }
    spectrumName += "_copy";
    histos[i] = static_cast<TH1*>(hs->Clone(spectrumName));

    if (eventScaling && params->GetEvents() > 0) {
      histos[i]->Scale(1. / params->GetEvents(), "width");
      histos[i]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d#it{z}_{D}}");
    }
  }

  TString cname(Form("fig_%s_DPtSpectraVsJetPt", params->GetName()));
  Bool_t res = PlotSpectra(params->GetNJetPtBins(), histos, cname, kTRUE);

  delete[] histos;

  return res;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDPtSpectraVsDz(DJetCorrAnalysisParams* params, Bool_t eventScaling)
{
  TH1** histos = new TH1*[params->GetNzBins()];
  
  for (Int_t i = 0; i < params->GetNzBins(); i++) {
    histos[i] = 0;
    
    TString spectrumCuts(params->GetCutString(kMatched, -1, -1, i));
    TString spectrumName(Form("h%s_DPtSpectrum_%s_Matched", params->GetName(), spectrumCuts.Data()));
    TH1* hs = dynamic_cast<TH1*>(fOutputList->FindObject(spectrumName));
    if (!hs) {
      Printf("Error-DJetCorrAnalysis::PlotDPtSpectraVsJetPt : Histogram '%s' not found!", spectrumName.Data());
      continue;
    }
    spectrumName += "_copy";
    histos[i] = static_cast<TH1*>(hs->Clone(spectrumName));
    if (eventScaling && params->GetEvents() > 0) {
      histos[i]->Scale(1. / params->GetEvents(), "width");
      histos[i]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d#it{z}_{D}}");
    }
  }

  TString cname(Form("fig_%s_DPtSpectraVsZ", params->GetName()));
  Bool_t res = PlotSpectra(params->GetNzBins(), histos, cname, kTRUE);

  delete[] histos;

  return res;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDPtSpectraVsMatchingStatus(DJetCorrAnalysisParams* params, Bool_t eventScaling)
{
  TH1* histos[3] = {0};
  
  TString spectrumCuts = params->GetCutString(kMatched, -1, -1, -1);
  
  TString spectrumNameMatched(Form("h%s_DPtSpectrum_%s_Matched", params->GetName(), spectrumCuts.Data()));
  TH1* hs0 = dynamic_cast<TH1*>(fOutputList->FindObject(spectrumNameMatched));
  if (hs0) {
    spectrumNameMatched += "_copy";
    histos[0] = static_cast<TH1*>(hs0->Clone(spectrumNameMatched));

    if (eventScaling && params->GetEvents() > 0) {
      histos[0]->Scale(1. / params->GetEvents(), "width");
      histos[0]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d#it{z}_{D}}");
    }

    histos[0]->SetTitle("Matched");
  }
  else {
    Printf("Error-DJetCorrAnalysis::PlotDPtSpectraVsMatchingStatus : Histogram '%s' not found!", spectrumNameMatched.Data());
  }

  spectrumCuts = params->GetCutString(kNotMatched, -1, -1, -1);

  TString spectrumNameNotMatched(Form("h%s_DPtSpectrum_%s_NotMatched", params->GetName(), spectrumCuts.Data()));
  TH1* hs1 = dynamic_cast<TH1*>(fOutputList->FindObject(spectrumNameNotMatched));
  if (hs1) {
    spectrumNameNotMatched += "_copy";
    histos[1] = static_cast<TH1*>(hs1->Clone(spectrumNameNotMatched));

    if (eventScaling && params->GetEvents() > 0) {
      histos[1]->Scale(1. / params->GetEvents(), "width");
      histos[1]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d#it{z}_{D}}");
    }

    histos[1]->SetTitle("Not Matched");
  }
  else {
    Printf("Error-DJetCorrAnalysis::PlotDPtSpectraVsMatchingStatus : Histogram '%s' not found!", spectrumNameNotMatched.Data());
  }

  spectrumCuts = params->GetCutString(kAnyMatchingStatus, -1, -1, -1);

  TString spectrumNameAnyMatchingStatus(Form("h%s_DPtSpectrum_%s_AnyMatchingStatus", params->GetName(), spectrumCuts.Data()));
  TH1* hs2 = dynamic_cast<TH1*>(fOutputList->FindObject(spectrumNameAnyMatchingStatus));
  if (hs1) {
    spectrumNameAnyMatchingStatus += "_copy";
    histos[2] = static_cast<TH1*>(hs2->Clone(spectrumNameAnyMatchingStatus));

    if (eventScaling && params->GetEvents() > 0) {
      histos[2]->Scale(1. / params->GetEvents(), "width");
      histos[2]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d#it{z}_{D}}");
    }

    histos[2]->SetTitle("All");
  }
  else {
    Printf("Error-DJetCorrAnalysis::PlotDPtSpectraVsMatchingStatus : Histogram '%s' not found!", spectrumNameNotMatched.Data());
  }
  
  TString cname(Form("fig_%s_DPtSpectraVsMatchingStatus", params->GetName()));
  return PlotSpectra(3, histos, cname, kTRUE);
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotDzSpectraVsJetPt(DJetCorrAnalysisParams* params, Bool_t eventScaling)
{
  TH1** histos = new TH1*[params->GetNJetPtBins()];
  
  for (Int_t i = 0; i < params->GetNJetPtBins(); i++) {
    histos[i] = 0;
    
    TString spectrumCuts(params->GetCutString(kMatched, -1, i, -1));
    TString spectrumName(Form("h%s_DzSpectrum_%s_Matched", params->GetName(), spectrumCuts.Data()));
    TH1* hs = dynamic_cast<TH1*>(fOutputList->FindObject(spectrumName));
    if (!hs) {
      Printf("Error-DJetCorrAnalysis::PlotDzSpectraVsJetPt : Histogram '%s' not found!", spectrumName.Data());
      continue;
    }
    spectrumName += "_copy";
    histos[i] = static_cast<TH1*>(hs->Clone(spectrumName));
    if (eventScaling && params->GetEvents() > 0) {
      histos[i]->Scale(1. / params->GetEvents(), "width");
      histos[i]->GetYaxis()->SetTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d#it{z}_{D}}");
    }
  }

  TString cname(Form("fig_%s_DzSpectraVsJetPt", params->GetName()));
  Bool_t res = PlotSpectra(params->GetNJetPtBins(), histos, cname, kFALSE);

  delete[] histos;

  return res;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotInvMassHistogramsVsDPt(DJetCorrAnalysisParams* params, EMatchingStatus st, Int_t jetptBin, Int_t dzBin)
{
  if (!fOutputList) return kFALSE;
  
  Double_t D0mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(421))->Mass();
  Double_t Dstarmass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(413))->Mass();
  
  TString fixedCuts(params->GetCutString(st, -1, jetptBin, dzBin));
  
  TString cname("fig_InvMassVsDPt");
  cname += params->GetName();
  cname += fixedCuts;

  TString matchString;
  TString jetCuts;
  if (st == kAnyMatchingStatus) {
    matchString = "AnyMatchingStatus";
    jetCuts = "";
  }
  else if (st == kNotMatched) {
    matchString = "NotMatched";
    jetCuts = "No matching jet";
  }
  else if (st == kMatched) {
    matchString = "Matched";
    if (jetptBin >= 0) {
      jetCuts = Form("%.1f < #it{p}_{T,jet}^{ch} < %.1f GeV/#it{c}", params->GetJetPtBin(jetptBin), params->GetJetPtBin(jetptBin+1));
    }
    else {
      jetCuts = Form("%.1f < #it{p}_{T,jet}^{ch} < %.1f GeV/#it{c}", params->GetMinJetPt(), params->GetMaxJetPt());
    }
  }

  cname += "_";
  cname += matchString;
  
  TString hname;
  TString xTitle;
  Double_t minMass = 0;
  Double_t maxMass = 0;
  Double_t pdgMass = -1;
  if (params->IsD0()) {
    xTitle = "#it{m}(K#pi) (GeV/#it{c}^{2})";
    hname = "InvMass";
    minMass = params->GetInvMinMass();
    maxMass = params->GetInvMaxMass();
    pdgMass = D0mass;
  }
  else if (params->IsDStar()) {
    xTitle = "#it{m}(K#pi#pi) - #it{m}(K#pi) (GeV/#it{c}^{2})";
    hname = "DeltaInvMass";
    pdgMass = Dstarmass - D0mass;
    minMass = params->GetDeltaInvMinMass();
    maxMass = params->GetDeltaInvMaxMass();
  }
  else {
    Printf("Error-DJetCorrAnalysis::PlotInvMassHistogramsVsDPt : Meson type '%s' not recognized!", params->GetDmesonName());
    return kFALSE;
  }
  
  TString prefix(params->GetName());

  TString spectrumCuts(params->GetCutString(st, -1, jetptBin, dzBin));
  TString spectrumName(Form("h%s_DPtSpectrum_%s_%s", prefix.Data(), spectrumCuts.Data(), matchString.Data()));

  TH1* histSpectrum = GetOutputHistogram(spectrumName);
  if ((histSpectrum && fForceRefit) || !histSpectrum) {
    if (histSpectrum) {
      fOutputList->Remove(histSpectrum);
      delete histSpectrum;
    }

    histSpectrum = new TH1D(spectrumName, jetCuts, params->GetNDPtBins(), params->GetDPtBins());
    histSpectrum->GetXaxis()->SetTitle("#it{p}_{T,D}");
    histSpectrum->GetYaxis()->SetTitle("counts");
    fOutputList->Add(histSpectrum);
  }
  else {
    histSpectrum = 0;
  }
  
  TH1** histos = new TH1*[params->GetNDPtBins()];
  Int_t n = 0;
  for (Int_t i = 0; i < params->GetNDPtBins(); i++) {
    TString cuts(params->GetCutString(st, i, jetptBin, dzBin));
    TString objname(Form("h%s_%s_%s_%s", prefix.Data(), hname.Data(), cuts.Data(), matchString.Data()));
    //Printf("Info-DJetCorrAnalysis::PlotInvMassHistogramsVsDPt : Retrieving histogram '%s'", objname.Data());
    TH1* hist = dynamic_cast<TH1*>(fOutputList->FindObject(objname));
    if (!hist) {
      Printf("Error-DJetCorrAnalysis::PlotInvMassHistogramsVsDPt : Histogram '%s' not found!", objname.Data());
      continue;
    }
    TString newName(objname);
    newName += "_copy";
    //Printf("Info-DJetCorrAnalysis::PlotInvMassHistogramsVsDPt : Cloning histogram '%s'", objname.Data());
    histos[n] = static_cast<TH1*>(hist->Clone(newName));
    //Printf("Info-DJetCorrAnalysis::PlotInvMassHistogramsVsDPt : Setting title of histogram '%s'", histos[i]->GetName());
    TString htitle(Form("%.1f < #it{p}_{T,D} < %.1f GeV/#it{c}", params->GetDPtBin(i), params->GetDPtBin(i+1)));
    histos[n]->SetTitle(htitle);
    n++;
  }

  TObjArray extraInfo;
  extraInfo.SetOwner(kTRUE);
  extraInfo.Add(new TObjString("ALICE Work in progress"));
  extraInfo.Add(new TObjString(jetCuts));
  
  Bool_t result = PlotInvMassHistogramArray(n, histos, cname, xTitle, minMass, maxMass, pdgMass, 0x0, 0x0, params, kTRUE, &extraInfo, histSpectrum);
  
  delete[] histos;
  
  return result;
}

//____________________________________________________________________________________
void DJetCorrAnalysis::GenerateMeaduredSpectrum(DJetCorrAnalysisParams* params)
{
  TString spectrum2Dcuts(params->GetCutString(kMatched, -1, -1, -1));
  TString spectrum2Dname(Form("h%s_DzSpectrum2D_%s_Matched", params->GetName(), spectrum2Dcuts.Data()));

  TH2* spectrum2D = static_cast<TH2*>(GetOutputHistogram(spectrum2Dname));
  if (spectrum2D) {
    if (fForceRefit) {
      fOutputList->Remove(spectrum2D);
      delete spectrum2D;
    }
    else {
      return;
    }
  }
  
  spectrum2D = new TH2D(spectrum2Dname, spectrum2Dname, params->GetNJetPtBins(), params->GetJetPtBins(), params->GetNzBins(), params->GetzBins());
  spectrum2D->GetXaxis()->SetTitle("#it{p}_{T,jet}^{det} GeV/#it{c}");
  spectrum2D->GetYaxis()->SetTitle("#it{z}_{||}^{det}");
  spectrum2D->GetZaxis()->SetTitle("Counts");
  
  for (Int_t j = 0; j < params->GetNJetPtBins(); j++) {
    TString spectrumCuts(params->GetCutString(kMatched, -1, j, -1));
    TString spectrumName(Form("h%s_DzSpectrum_%s_Matched", params->GetName(), spectrumCuts.Data()));

    TH1* spectrum = GetOutputHistogram(spectrumName);

    if (!spectrum) continue;
    
    Int_t jetbin = j+1;
    for (Int_t zbin = 0; zbin <= spectrum->GetNbinsX()+1; zbin++) {
      spectrum2D->SetBinContent(jetbin, zbin, spectrum->GetBinContent(zbin));
      spectrum2D->SetBinError(jetbin, zbin, spectrum->GetBinError(zbin));
    }
  }

  fOutputList->Add(spectrum2D);
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotInvMassHistogramsVsDz(DJetCorrAnalysisParams* params, Int_t dptBin, Int_t jetptBin)
{
  if (!fOutputList) return kFALSE;
  
  Double_t D0mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(421))->Mass();
  Double_t Dstarmass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(413))->Mass();
  
  TString fixedCuts(params->GetCutString(kMatched, dptBin, jetptBin, -1));
  fixedCuts += "_Matched";
  
  TString cname;
  TString hname;
  TString xTitle;
  Double_t minMass = 0;
  Double_t maxMass = 0;
  Double_t pdgMass = -1;
  TH1** histos = 0;
    
  TString cname2;
  TString hname2;
  TString xTitle2;
  Double_t minMass2 = 0;
  Double_t maxMass2 = 0;
  Double_t pdgMass2 = -1;
  TH1** histos2 = 0;
  
  if (params->IsD0()) {
    cname = "fig_InvMassVsDz";
    xTitle = "#it{m}(K#pi) (GeV/#it{c}^{2})";
    hname = "InvMass";
    pdgMass = D0mass;
    minMass = params->GetInvMinMass();
    maxMass = params->GetInvMaxMass();
    histos = new TH1*[params->GetNzBins()];
  }
  else if (params->IsDStar()) {
    cname = "fig_DeltaInvMassVsDz";
    xTitle = "#it{m}(K#pi#pi) - #it{m}(K#pi) (GeV/#it{c}^{2})";
    hname = "DeltaInvMass";
    pdgMass = Dstarmass - D0mass;
    minMass = params->GetDeltaInvMinMass();
    maxMass = params->GetDeltaInvMaxMass();
    histos = new TH1*[params->GetNzBins()];

    cname2 = "fig_D0InvMassVsDz";
    xTitle2 = "#it{m}(K#pi) (GeV/#it{c}^{2})";
    hname2 = "D0InvMass";
    pdgMass2 = D0mass;
    minMass2 = params->Get2ProngMinMass();
    maxMass2 = params->Get2ProngMaxMass();
    histos2 = new TH1*[params->GetNzBins()];
    cname2 += params->GetName();
    cname2 += fixedCuts;
  }
  else {
    Printf("Error-DJetCorrAnalysis::PlotInvMassHistogramsVsDz : Meson type '%s' not recognized!", params->GetDmesonName());
    return kFALSE;
  }

  cname += params->GetName();
  cname += fixedCuts;

  TString prefix(params->GetName());

  TString jetCuts;
  if (jetptBin >= 0) {
    jetCuts = Form("%.1f < #it{p}_{T,jet}^{ch} < %.1f GeV/#it{c}", params->GetJetPtBin(jetptBin), params->GetJetPtBin(jetptBin+1));
  }
  else {
    jetCuts = Form("%.1f < #it{p}_{T,jet}^{ch} < %.1f GeV/#it{c}", params->GetMinJetPt(), params->GetMaxJetPt());
  }

  TString spectrumCuts(params->GetCutString(kMatched, dptBin, jetptBin, -1));
  TString spectrumName(Form("h%s_DzSpectrum_%s_Matched", prefix.Data(), spectrumCuts.Data()));

  TH1* histSpectrum = GetOutputHistogram(spectrumName);
  if ((histSpectrum && fForceRefit) || !histSpectrum) {
    if (histSpectrum) {
      fOutputList->Remove(histSpectrum);
      delete histSpectrum;
    }

    histSpectrum = new TH1D(spectrumName, jetCuts, params->GetNzBins(), params->GetzBins());
    histSpectrum->GetXaxis()->SetTitle("#it{z}_{||}");
    histSpectrum->GetYaxis()->SetTitle("counts");
    fOutputList->Add(histSpectrum);
  }
  else {
    histSpectrum = 0;
  }

  Double_t* minMassSel = new Double_t[params->GetNzBins()];
  Double_t* maxMassSel = new Double_t[params->GetNzBins()];

  Double_t meanJetPt = params->GetMeanJetPt(jetptBin);

  Int_t n = 0;
  for (Int_t i = 0; i < params->GetNzBins(); i++) {
    TString cuts(params->GetCutString(kMatched, dptBin, jetptBin, i));

    Double_t meanDPt = meanJetPt * params->GetMeanZ(i);

    minMassSel[i] = params->Get2ProngMinMass(meanDPt);
    maxMassSel[i] = params->Get2ProngMaxMass(meanDPt);
    
    TString objname(Form("h%s_%s_%s_Matched", prefix.Data(), hname.Data(), cuts.Data()));
    TH1* hist = dynamic_cast<TH1*>(fOutputList->FindObject(objname));
    if (!hist) {
      Printf("Error-DJetCorrAnalysis::PlotInvMassHistogramsVsDz : Histogram '%s' not found!", objname.Data());
      continue;
    }
    TString newName(objname);
    newName += "_copy";
    histos[n] = static_cast<TH1*>(hist->Clone(newName));
    TString htitle(Form("%.1f < #it{z}_{D} < %.1f", params->GetzBin(i), params->GetzBin(i+1)));
    histos[n]->SetTitle(htitle);

    if (histos2) {
      TString objname2(Form("h%s_%s_%s_Matched", prefix.Data(), hname2.Data(), cuts.Data()));
      TH1* hist2 = dynamic_cast<TH1*>(fOutputList->FindObject(objname2));
      if (!hist2) {
        Printf("Error-DJetCorrAnalysis::PlotInvMassHistogramsVsDz : Histogram '%s' not found!", objname2.Data());
        continue;
      }
      TString newName2(objname2);
      newName2 += "_copy";
      histos2[n] = static_cast<TH1*>(hist2->Clone(newName2));
      TString htitle2(Form("%.1f < #it{z}_{D} < %.1f", params->GetzBin(i), params->GetzBin(i+1)));
      histos2[n]->SetTitle(htitle2);
    }
    
    n++;
  }
  
  Bool_t result = kTRUE;
  
  TObjArray extraInfo;
  extraInfo.SetOwner(kTRUE);
  extraInfo.Add(new TObjString("ALICE Work in progress"));
  extraInfo.Add(new TObjString(jetCuts));

  result = PlotInvMassHistogramArray(n, histos, cname, xTitle, minMass, maxMass, pdgMass, 0x0, 0x0, params, kTRUE, &extraInfo, histSpectrum);
  delete[] histos;

  if (histos2) {
    result = PlotInvMassHistogramArray(n, histos2, cname2, xTitle2, minMass2, maxMass2, pdgMass2, minMassSel, maxMassSel, params, kFALSE, &extraInfo, 0x0) && result;
    delete[] histos2;
  }

  delete[] minMassSel;
  delete[] maxMassSel;

  return result;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotSpectra(Int_t n, TH1** histSpectra, const char* name, Bool_t logY)
{
  // Plot spectra.

  TString cname(name);

  if (n <= 0) return kFALSE;

  if (!histSpectra[0]) return kFALSE;

  HistoStyler styler;
  styler.SetMarkerStyle(kFullCircle);
  styler.SetMarkerSize(0.8);
  styler.SetVariableMarkerColor();
  styler.SetVariableLineColor();
  styler.SetLineWidth(1);
  styler.Apply(n, histSpectra);
  
  Double_t min = histSpectra[0]->GetMinimum();
  Double_t max = histSpectra[0]->GetMaximum();

  if (logY && min <= 0) min = fgkEpsilon;
  if (max < min) max = min;

  TCanvas* canvas = SetUpCanvas(cname,
                                histSpectra[0]->GetXaxis()->GetTitle(), histSpectra[0]->GetXaxis()->GetXmin(), histSpectra[0]->GetXaxis()->GetXmax(), kFALSE,
                                histSpectra[0]->GetYaxis()->GetTitle(), min, max, logY);

  TLegend* leg = SetUpLegend(0.58, 0.68, 0.88, 0.88, 14);

  Double_t extraFactor = 1.5;
  if (logY) extraFactor = 2.0;
  
  for (Int_t i = 0; i < n; i++) {
    if (!histSpectra[i]) continue;

    FitHistogramInPad(histSpectra[i], canvas, "", kFALSE, extraFactor);

    leg->AddEntry(histSpectra[i], histSpectra[i]->GetTitle(), "pe");
  }


  leg->Draw();

  if (fSavePlots) SavePlot(canvas);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::PlotInvMassHistogramArray(Int_t n, TH1** histos,
                                                   const char* name, const char* xTitle,
                                                   Double_t minMass, Double_t maxMass, Double_t pdgMass, Double_t* minMassSel, Double_t* maxMassSel,
                                                   DJetCorrAnalysisParams* params, Bool_t doFit, TObjArray* extraInfo, TH1* histSpectrum)
{
  // Plot invariant mass histograms contained in histos.

  if (n == 0) return kFALSE;
  
  Printf("Info-DJetCorrAnalysis::PlotInvMassHistograms : Plotting invariant mass histograms '%s'", name);
  
  Int_t cols = TMath::FloorNint(TMath::Sqrt(n));
  Int_t rows = TMath::CeilNint(1. * n / cols);

  TString cname(name);

  Double_t w = cols*350;
  Double_t h = rows*320;

  TString yaxisTitle;
  yaxisTitle = "counts";
  
  TCanvas* canvas = SetUpCanvas(cname, xTitle, minMass, maxMass, kFALSE, yaxisTitle, 0, 1, kFALSE, h, w, cols, rows, 0.18, 0.02, 0.12, 0.08);
  for (Int_t i = 0; i < n; i++) {    
    TVirtualPad* pad = canvas->cd(i+1);
    TH1* blankHist = dynamic_cast<TH1*>(pad->GetListOfPrimitives()->At(0));
    if (!blankHist) {
      Printf("Error-DJetCorrAnalysis::PlotInvMassHistograms : Could not find blank histogram!");
      continue;
    }

    blankHist->GetXaxis()->SetTitleFont(43);
    blankHist->GetXaxis()->SetTitleSize(15);
    blankHist->GetXaxis()->SetTitleOffset(2.5);
    blankHist->GetXaxis()->SetLabelFont(43);
    blankHist->GetXaxis()->SetLabelSize(15);
    blankHist->GetXaxis()->SetNdivisions(404, kTRUE);

    blankHist->GetYaxis()->SetTitleFont(43);
    blankHist->GetYaxis()->SetTitleSize(15);
    blankHist->GetYaxis()->SetTitleOffset(3.8);
    blankHist->GetYaxis()->SetLabelFont(43);
    blankHist->GetYaxis()->SetLabelSize(15);

    if (histos[i]->GetSumw2N() == 0) histos[i]->Sumw2();
    Printf("Info-DJetCorrAnalysis::PlotInvMassHistograms : Now plotting '%s'", histos[i]->GetName());
    histos[i]->SetMarkerStyle(kFullCircle);
    histos[i]->SetMarkerSize(0.6);
    histos[i]->SetMarkerColor(kBlue+3);
    histos[i]->SetLineColor(kBlue+3);
    TH1* hcopy = histos[i]->DrawCopy("same p");

    MassFitter* fitter = 0;
    
    if (doFit) {      
      TString fitterName(Form("%s_fitter", histos[i]->GetName()));

      fitter = static_cast<MassFitter*>(fMassFitters->FindObject(fitterName));
      if (fForceRefit) {
        fMassFitters->Remove(fitter);
        delete fitter;
        fitter = 0;
      }
      
      if (!fitter) {
        fitter = params->CreateMassFitter(fitterName);
        fMassFitters->Add(fitter);

        Double_t integral = histos[i]->Integral(histos[i]->GetXaxis()->FindBin(minMass), histos[i]->GetXaxis()->FindBin(maxMass));
        //if (!params->IsBackgroundOnly()) {
        fitter->GetFitFunction()->FixParameter(0, integral); // total integral is fixed
        fitter->GetFitFunction()->SetParameter(2, integral / 100); // signal integral (start with very small signal)
        fitter->GetFitFunction()->SetParLimits(2, 0, integral); // signal integral has to be contained in the total integral
        fitter->GetFitFunction()->SetParameter(3, pdgMass); // start fitting using PDG mass
        //}
        //else {
        //  fitter->GetFitFunction()->SetParameter(0, integral); // total integral is NOT fixed
        // }

        fitter->Fit(histos[i], "0 E S");
      }

      TFitResultPtr r = fitter->GetFitStatus();
      Int_t fitStatus = r;

      if (histSpectrum) {
        //if (params->IsBkgSub() && !params->IsBackgroundOnly()) {
        //histSpectrum->SetBinContent(i+1, fitter->GetTotalEntries());
        // histSpectrum->SetBinError(i+1, fitter->GetTotalEntriesError());
        //}
        // else {
        if (fitStatus == 0) {
          histSpectrum->SetBinContent(i+1, fitter->GetSignal());
          histSpectrum->SetBinError(i+1, fitter->GetSignalError());
        }
        //}
      }

      fitter->Draw("same");

      TPaveText* paveSig = SetUpPaveText(0.22, 0.41, 0.51, 0.84, 13, fitter->GetSignalString());
      paveSig->AddText(fitter->GetBackgroundString());
      paveSig->AddText(fitter->GetSignalOverBackgroundString());
      paveSig->AddText(fitter->GetSignalOverSqrtSignalBackgroundString());
      if (fitStatus == 0) {
        paveSig->AddText(fitter->GetChisquareString());
      }
      else {
        paveSig->AddText("Fit failed");
      }
      paveSig->Draw();

      TPaveText* paveFit = SetUpPaveText(0.47, 0.58, 0.97, 0.84, 13, fitter->GetSignalMeanString());
      paveFit->AddText(fitter->GetSignalWidthString());
      paveFit->AddText(fitter->GetBkgPar1String());
      paveFit->AddText(fitter->GetTotalEntriesString());
      paveFit->Draw();

      if (fInvMassPlotNorm == kNormalizeBackground) {
        Double_t scaleFactor = 1. / fitter->GetBackground();
        histos[i]->Scale(scaleFactor, "width");
        hcopy->Scale(scaleFactor, "width");
        fitter->NormalizeBackground();
        yaxisTitle = "arb. units";
      }
    }

    if (fInvMassPlotNorm == kDivideByBinWidth) {
      hcopy->Scale(1., "width");
      histos[i]->Scale(1., "width");
      if (fitter) fitter->DivideByBinWidth();
      yaxisTitle = Form("counts / (%.2f MeV/#it{c}^{2})", hcopy->GetXaxis()->GetBinWidth(1)*1000);
    }

    blankHist->GetYaxis()->SetRangeUser(0, hcopy->GetMaximum()*1.8);
    blankHist->GetYaxis()->SetTitle(yaxisTitle);
      
    TPaveText* pave = SetUpPaveText(gPad->GetLeftMargin(), 0.81, 1-gPad->GetRightMargin(), 0.95, 15, histos[i]->GetTitle());
    pave->SetTextAlign(22);
    pave->SetTextFont(63);
    pave->Draw();

    if (pdgMass > 0 && !doFit) {
      TLine *line = new TLine(pdgMass, 0, pdgMass, hcopy->GetMaximum()*0.6);
      line->SetLineColor(kRed);
      line->SetLineWidth(1);
      line->Draw();
    }
    if (maxMassSel && minMassSel) {
      TLine *line1 = new TLine(minMassSel[i], 0, minMassSel[i], hcopy->GetMaximum());
      line1->SetLineColor(kGreen+2);
      line1->SetLineWidth(1);
      line1->SetLineStyle(2);
      line1->Draw();

      TLine *line2 = new TLine(maxMassSel[i], 0, maxMassSel[i], hcopy->GetMaximum());
      line2->SetLineColor(kGreen+2);
      line2->SetLineWidth(1);
      line2->SetLineStyle(2);
      line2->Draw();
    }
  }

  if (extraInfo) {
    for (Int_t i = 0; i < extraInfo->GetEntriesFast() && i < n; i++) {
      TObjString* info = static_cast<TObjString*>(extraInfo->At(i));
      canvas->cd(i+1);
      TPaveText* paveText = SetUpPaveText(0.18, 0.94, 1.0, 0.98, 12, info->GetString());
      if (i == 0) {
        paveText->SetTextFont(63);
        paveText->SetTextSize(15);
      }
      paveText->SetTextAlign(22);
      paveText->Draw();
    }
  }
  
  if (fSavePlots) SavePlot(canvas);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectQA()
{
  // Project QA histograms.

  TH1* evCount = dynamic_cast<TH1*>(fInputQAList->FindObject("fHistEventCount"));
  evCount->SetName("hEventCount");
  fOutputList->Add(evCount);

  TString labels[4] = {"Global", "Constr", "ConstrNoITS", "Hybrid"};

  TH1* hTrEtaSum = 0;
  TH1* hTrPhiSum = 0;
  TH1* hTrPtSum = 0;
  
  for (Int_t i = 0; i < 3; i++) {
    TH3* hTrPhiEtaPt = dynamic_cast<TH3*>(fInputQAList->FindObject(Form("fHistTrPhiEtaPt_0_%d",i)));
    if (!hTrPhiEtaPt) {
      Printf("Error-DJetCorrAnalysis::ProjectQA : Could not find track histogram for track type %d!", i);
      continue;
    }
    TH1* hTrEta = hTrPhiEtaPt->ProjectionX(Form("hTracks_%s_Eta", labels[i].Data()));
    TH1* hTrPhi = hTrPhiEtaPt->ProjectionY(Form("hTracks_%s_Phi", labels[i].Data()));
    TH1* hTrPt = hTrPhiEtaPt->ProjectionZ(Form("hTracks_%s_Pt", labels[i].Data()));

    if (i == 0) {
      hTrEtaSum = static_cast<TH1*>(hTrEta->Clone(Form("hTracks_%s_Eta", labels[3].Data())));
      hTrPhiSum = static_cast<TH1*>(hTrPhi->Clone(Form("hTracks_%s_Phi", labels[3].Data())));
      hTrPtSum = static_cast<TH1*>(hTrPt->Clone(Form("hTracks_%s_Pt", labels[3].Data())));

      fOutputList->Add(hTrEtaSum);
      fOutputList->Add(hTrPhiSum);
      fOutputList->Add(hTrPtSum);
    }
    else {
      hTrEtaSum->Add(hTrEta);
      hTrPhiSum->Add(hTrPhi);
      hTrPtSum->Add(hTrPt);
    }
    
    fOutputList->Add(hTrEta);
    fOutputList->Add(hTrPhi);
    fOutputList->Add(hTrPt);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectCorrD(DJetCorrAnalysisParams* params)
{
  // Project histograms related to the D meson correlated to a jet.

  TString cutsN;
  TString cutsD;
  
  TString prefix(params->GetName());

  Printf("Info-DJetCorrAnalysis::ProjectCorrD : Entering method.");

  ProjectDJetCorr(prefix, "AnyMatchingStatus", params, kAnyMatchingStatus, -1, -1, -1);
  ProjectDJetCorr(prefix, "NotMatched", params, kNotMatched, -1, -1, -1);
  ProjectDJetCorr(prefix, "Matched", params, kMatched, -1, -1, -1);

  TString dCuts(Form("DPt_%02.0f_%02.0f", params->GetMinDPt(), params->GetMaxDPt()));
  dCuts.ReplaceAll(".", "");

  TString jetCuts(Form("JetPt_%03.0f_%03.0f", params->GetMinJetPt(), params->GetMaxJetPt()));
  jetCuts.ReplaceAll(".", "");
  
  TString zCuts(Form("z_%.1f_%.1f", params->GetMinZ(), params->GetMaxZ()));
  zCuts.ReplaceAll(".", "");
  
  cutsN = Form("%s_%s_%s_Matched", dCuts.Data(), jetCuts.Data(), zCuts.Data());
  cutsD = Form("%s_AnyMatchingStatus", dCuts.Data());
  GenerateRatios(cutsN, cutsD);

  for (Int_t i = 0; i < params->GetNDPtBins(); i++) {
    ProjectDJetCorr(prefix, "AnyMatchingStatus", params, kAnyMatchingStatus, i, -1, -1, 0);
    ProjectDJetCorr(prefix, "NotMatched", params, kNotMatched, i, -1, -1, 0);
    ProjectDJetCorr(prefix, "Matched", params, kMatched, i, -1, -1, 0);

    dCuts = Form("DPt_%02.0f_%02.0f", params->GetDPtBin(i), params->GetDPtBin(i+1));
    dCuts.ReplaceAll(".", "");

    cutsN = Form("%s_%s_%s_Matched", dCuts.Data(), jetCuts.Data(), zCuts.Data());
    cutsD = Form("%s_AnyMatchingStatus", dCuts.Data());
    GenerateRatios(cutsN, cutsD);
  }

  for (Int_t i = 0; i < params->GetNJetPtBins(); i++) {
    ProjectDJetCorr(prefix, "Matched", params, kMatched, -1, i);
  }

  for (Int_t i = 0; i < params->GetNzBins(); i++) {
    ProjectDJetCorr(prefix, "Matched", params, kMatched, -1, -1, i);
    for (Int_t j = 0; j < params->GetNJetPtBins(); j++) {
      ProjectDJetCorr(prefix, "Matched", params, kMatched, -1, j, i);
    }
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectDJetCorr(TString prefix, TString suffix,
                                         DJetCorrAnalysisParams* params, EMatchingStatus st,
                                         Int_t dptBin, Int_t jetptBin, Int_t dzBin, Int_t minJetConst)
{
  // Project histograms related to the D meson with specified cuts.
  
  if (!LoadTHnSparse()) {
    return kFALSE;
  }

  Printf("Info-DJetCorrAnalysis::ProjectDJetCorr : Method called with prefix = '%s', suffix = '%s'", prefix.Data(), suffix.Data());

  TString cuts(params->GetCutString(st, dptBin, jetptBin, dzBin));

  Double_t minDPt = 0;
  Double_t maxDPt = 0;
  params->GetDPtBinRange(minDPt, maxDPt, dptBin);

  Double_t minz = 0;
  Double_t maxz = 0;
  params->GetzBinRange(minz, maxz, dzBin);

  Double_t minJetPt = 0;
  Double_t maxJetPt = 0;
  params->GetJetPtBinRange(minJetPt, maxJetPt, jetptBin);

  if (!suffix.IsNull()) suffix.Prepend("_");

  Int_t dPtAxis = params->GetAxisIndex(fDPtAxisTitle);
  if (dPtAxis < 0) return kFALSE;
  
  Int_t dEtaAxis = params->GetAxisIndex(fDEtaAxisTitle);
  if (dEtaAxis < 0) return kFALSE;
  
  Int_t dPhiAxis = params->GetAxisIndex(fDPhiAxisTitle);
  if (dPhiAxis < 0) return kFALSE;
  
  Int_t jetPtAxis = params->GetAxisIndex(fJetPtAxisTitle);
  if (jetPtAxis < 0) return kFALSE;
    
  Int_t jetEtaAxis = params->GetAxisIndex(fJetEtaAxisTitle);
  if (jetEtaAxis < 0) return kFALSE;
  
  Int_t jetPhiAxis = params->GetAxisIndex(fJetPhiAxisTitle);
  if (jetPhiAxis < 0) return kFALSE;

  Int_t dInvMassAxis = params->GetAxisIndex(fDInvMassAxisTitle);
  Int_t d2ProngInvMassAxis = params->GetAxisIndex(fD2prongInvMassAxisTitle);
  Int_t dDeltaInvMassAxis = params->GetAxisIndex(fDDeltaInvMassAxisTitle);
  Int_t dSoftPionPtAxis = params->GetAxisIndex(fDSoftPionPtAxisTitle);

  Int_t deltaRAxis = params->GetAxisIndex(fDeltaRAxisTitle);
  Int_t deltaEtaAxis = params->GetAxisIndex(fDeltaEtaAxisTitle);
  Int_t deltaPhiAxis = params->GetAxisIndex(fDeltaPhiAxisTitle);
  Int_t dzAxis = params->GetAxisIndex(fDzAxisTitle);
  Int_t jetConstAxis = params->GetAxisIndex(fJetConstAxisTitle);
  
  Int_t matchingStatusAxis = params->GetAxisIndex(fMatchingStatusAxisTitle);
  
  THnSparse* Dmesons = params->GetTHnSparse();
  
  if (matchingStatusAxis >= 0) {
    if (st == kMatched) {
      Dmesons->GetAxis(matchingStatusAxis)->SetRange(1,2);
    }
    else if (st == kNotMatched) {
      Dmesons->GetAxis(matchingStatusAxis)->SetRange(3,4);
    }
    else {
      Dmesons->GetAxis(matchingStatusAxis)->SetRange(1,4);
    }
  }

  Dmesons->GetAxis(dPtAxis)->SetRangeUser(minDPt, maxDPt);

  if (st == kMatched) { // apply cuts on z and jet pt only if matched requested
    Dmesons->GetAxis(jetPtAxis)->SetRangeUser(minJetPt, maxJetPt);

    if (dzAxis >= 0) Dmesons->GetAxis(dzAxis)->SetRangeUser(minz, maxz);
    
    if (jetConstAxis >= 0) {
      if (minJetConst > 0) {
        Dmesons->GetAxis(jetConstAxis)->SetRangeUser(params->GetMinJetConstituents(), Dmesons->GetAxis(jetConstAxis)->GetXmax()+fgkEpsilon);
      }
      else if (minJetConst < 0) {
        Dmesons->GetAxis(jetConstAxis)->SetRangeUser(0, params->GetMinJetConstituents());
      }
      else {
        Dmesons->GetAxis(jetConstAxis)->SetRange(0, Dmesons->GetAxis(jetConstAxis)->GetNbins()+1);
      }
    }
  }
  else {
    Dmesons->GetAxis(dEtaAxis)->SetRangeUser(params->GetMinDEta(), params->GetMaxDEta()); // this is for non-matched D mesons: look only for D mesons in a defined eta region (usually the eta region is defined by the jet axis)
    
    if (jetPtAxis >= 0) Dmesons->GetAxis(jetPtAxis)->SetRange(0, Dmesons->GetAxis(jetPtAxis)->GetNbins()+1);
    if (dzAxis >= 0) Dmesons->GetAxis(dzAxis)->SetRange(0, Dmesons->GetAxis(dzAxis)->GetNbins()+1);
    if (jetConstAxis >= 0) Dmesons->GetAxis(jetConstAxis)->SetRange(0, Dmesons->GetAxis(jetConstAxis)->GetNbins()+1);
  }

  if (dInvMassAxis >= 0) {
    TH1* hdinvmass = Dmesons->Projection(dInvMassAxis, "EA");
    hdinvmass->SetName(Form("h%s_InvMass_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
    Printf("Info-DJetCorrAnalysis::ProjectDJetCorr : Adding histogram '%s'", hdinvmass->GetName());
    if (params->GetInvMassRebinFactor() > 1 && params->IsD0()) {
      Printf("Rebinning %s by a factor %d", hdinvmass->GetName(), params->GetInvMassRebinFactor());
      hdinvmass->Rebin(params->GetInvMassRebinFactor());
    }
    fOutputList->Add(hdinvmass);
    
    Dmesons->GetAxis(dInvMassAxis)->SetRangeUser(params->GetInvMinMass(), params->GetInvMaxMass());
  }

  if (d2ProngInvMassAxis >= 0) {
    TH1* hd2pronginvmass = Dmesons->Projection(d2ProngInvMassAxis, "EA");
    hd2pronginvmass->SetName(Form("h%s_D0InvMass_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
    Printf("Info-DJetCorrAnalysis::ProjectDJetCorr : Adding histogram '%s'", hd2pronginvmass->GetName());
    fOutputList->Add(hd2pronginvmass);

    Double_t meanDPt = params->GetMeanDPt(dzBin, jetptBin);

    Dmesons->GetAxis(d2ProngInvMassAxis)->SetRangeUser(params->Get2ProngMinMass(meanDPt), params->Get2ProngMaxMass(meanDPt));
  }

  if (dDeltaInvMassAxis >= 0) {
    TH1* hddeltainvmass = Dmesons->Projection(dDeltaInvMassAxis, "EA");
    hddeltainvmass->SetName(Form("h%s_DeltaInvMass_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
    Printf("Info-DJetCorrAnalysis::ProjectDJetCorr : Adding histogram '%s'", hddeltainvmass->GetName());
    if (params->GetInvMassRebinFactor() > 1) {
      Printf("Rebinning %s by a factor %d", hddeltainvmass->GetName(), params->GetInvMassRebinFactor());
      hddeltainvmass->Rebin(params->GetInvMassRebinFactor());
    }
    fOutputList->Add(hddeltainvmass);

    if (dSoftPionPtAxis >= 0) {
      TH2* hdsoftpionptVsDeltaInvMass = Dmesons->Projection(dSoftPionPtAxis, dDeltaInvMassAxis, "EO");
      hdsoftpionptVsDeltaInvMass->SetName(Form("h%s_SoftPionVsDeltaInvMass_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hdsoftpionptVsDeltaInvMass);
    }
    
    Dmesons->GetAxis(dDeltaInvMassAxis)->SetRangeUser(params->GetDeltaInvMinMass(), params->GetDeltaInvMaxMass());
  }

  TH2* hDpos = Dmesons->Projection(dPhiAxis, dEtaAxis, "EO");
  hDpos->SetName(Form("h%s_MesonPhiVsEta_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
  fOutputList->Add(hDpos);

  if (dSoftPionPtAxis >= 0) {
    TH1* hdsoftpionpt = Dmesons->Projection(dSoftPionPtAxis, "EO");
    hdsoftpionpt->SetName(Form("h%s_SoftPion_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
    fOutputList->Add(hdsoftpionpt);
  }

  TString hname;
  
  if (dPtAxis >=0 && dptBin == -1) {
    hname = Form("h%s_MesonPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data());
    if (!fOutputList->Contains(hname)) {
      TH1* hDpt = Dmesons->Projection(dPtAxis, "EO");
      hDpt->SetName(hname);
      fOutputList->Add(hDpt);

      if (dSoftPionPtAxis >= 0) {
        TH2* hdsoftpionptVsDpt = Dmesons->Projection(dSoftPionPtAxis, dPtAxis, "EO");
        hdsoftpionptVsDpt->SetName(Form("h%s_SoftPionVsDPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hdsoftpionptVsDpt);
      }
    }
  }
  
  if (st == kMatched) {
    if (deltaRAxis >= 0) {
      TH1* hDeltaR = Dmesons->Projection(deltaRAxis, "EO");
      hDeltaR->SetName(Form("h%s_DeltaR_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDeltaR);

      if (dPtAxis >=0 && dptBin == -1) {
        TH2* hDeltaRVsDpt = Dmesons->Projection(deltaRAxis, dPtAxis, "EO");
        hDeltaRVsDpt->SetName(Form("h%s_DeltaRvsDPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hDeltaRVsDpt);
      }
    }

    if (jetConstAxis >= 0) {
      TH1* hJetConst = Dmesons->Projection(jetConstAxis, "EO");
      hJetConst->SetName(Form("h%s_JetConstituents_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hJetConst);

      if (dzAxis >= 0 && dzBin == -1) {
        TH2* hJetConstVsDz = Dmesons->Projection(jetConstAxis, dzAxis, "EO");
        hJetConstVsDz->SetName(Form("h%s_JetConstituentsVsDz_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hJetConstVsDz);
      }

      if (jetPtAxis >= 0 && jetptBin == -1) {
        TH2* hJetConstVsJetPt = Dmesons->Projection(jetConstAxis, jetPtAxis, "EO");
        hJetConstVsJetPt->SetName(Form("h%s_JetConstituentsVsJetPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hJetConstVsJetPt);
      }
    }

    if (deltaEtaAxis >= 0) {
      TH1* hDeltaEta = Dmesons->Projection(deltaEtaAxis, "EO");
      hDeltaEta->SetName(Form("h%s_DeltaEta_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDeltaEta);
    }

    if (deltaPhiAxis >= 0) {
      TH1* hDeltaPhi = Dmesons->Projection(deltaPhiAxis, "EO");
      hDeltaPhi->SetName(Form("h%s_DeltaPhi_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hDeltaPhi);
    }
    
    if (dzAxis >= 0 && dzBin == -1) {
      TH1* hdz = Dmesons->Projection(dzAxis, "EO");
      hdz->SetName(Form("h%s_MesonZ_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hdz);

      if (dPtAxis >=0 && dptBin == -1) {
        TH2* hdzVsDPt = Dmesons->Projection(dzAxis, dPtAxis, "EO");
        hdzVsDPt->SetName(Form("h%s_MesonZvsDPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hdzVsDPt);
      }
      
      if (jetPtAxis >= 0 && jetptBin == -1) {
        TH2* hdzVsJetPt = Dmesons->Projection(dzAxis, jetPtAxis, "EO");
        hdzVsJetPt->SetName(Form("h%s_MesonZvsJetPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hdzVsJetPt);
      }
      
      if (deltaRAxis >= 0) {
        TH2* hdzVsDeltaR = Dmesons->Projection(dzAxis, deltaRAxis, "EO");
        hdzVsDeltaR->SetName(Form("h%s_MesonZvsDeltaR_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hdzVsDeltaR);
      }
    }

    if (jetPtAxis >= 0 && jetptBin == -1) {
      TH1* hJetPt = Dmesons->Projection(jetPtAxis, "EO");
      hJetPt->SetName(Form("h%s_JetPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hJetPt);

      if (dPtAxis >=0 && dptBin == -1) {
        TH2* hDPtVsJetPt = Dmesons->Projection(dPtAxis, jetPtAxis, "EO");
        hDPtVsJetPt->SetName(Form("h%s_MesonPtvsJetPt_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
        fOutputList->Add(hDPtVsJetPt);
      }
    }

    if (jetEtaAxis >= 0 && jetPhiAxis >= 0) {
      TH2* hJetpos = Dmesons->Projection(jetPhiAxis, jetEtaAxis, "EO");
      hJetpos->SetName(Form("h%s_JetPhiVsEta_%s%s", prefix.Data(), cuts.Data(), suffix.Data()));
      fOutputList->Add(hJetpos);
    }
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::LoadQAList()
{
  // Load the input QA list.

  if (!OpenInputFile()) return kFALSE;
  
  if (!fInputQAList || fQAListName != fInputQAList->GetName()) {
    delete fInputQAList;
    Printf("Info-DJetCorrAnalysis::OpenInputFile : Getting list '%s' from file '%s'", fQAListName.Data(), fInputFile->GetName()); 
    fInputQAList = dynamic_cast<TList*>(fInputFile->Get(fQAListName));
  }
  
  if (!fInputQAList) {
    Printf("Error-DJetCorrAnalysis::OpenInputFile : Could not get list '%s' from file '%s'", fQAListName.Data(), fInputFile->GetName());
    return kFALSE;
  }
    
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::Regenerate()
{
  Bool_t result = kFALSE;

  result = GenerateDJetCorrHistograms();
  if (!result) return kFALSE;

  result = PlotDJetCorrHistograms(kTRUE);
  if (!result) return kFALSE;

  return kTRUE;
}

//____________________________________________________________________________________
TString DJetCorrAnalysis::GetTruthName(Int_t p)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";

  TString spectrum2Dcuts(params->GetCutString(kMatched, -1, -1, -1));
  TString hname = Form("h%s_MesonZvsJetPt_%s_Matched_Truth_Coarse", params->GetName(), spectrum2Dcuts.Data());
  
  return hname;
}

//____________________________________________________________________________________
TString DJetCorrAnalysis::GetMeasuredName(Int_t p)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";

  TString spectrum2Dcuts(params->GetCutString(kMatched, -1, -1, -1));
  TString hname(Form("h%s_DzSpectrum2D_%s_Matched", params->GetName(), spectrum2Dcuts.Data()));
  
  return hname;
}

//____________________________________________________________________________________
TString DJetCorrAnalysis::GetDzTruthName(Int_t /*p*/)
{
  //To do
  return "";
}

//____________________________________________________________________________________
TString DJetCorrAnalysis::GetDzMeasuredName(Int_t /*p*/)
{
  //To do
  return "";
}

//____________________________________________________________________________________
TString DJetCorrAnalysis::GetDPtTruthName(Int_t p, const char* matching)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";

  TString spectrumCuts;

  if (strcmp("Matched", matching)==0) {
    spectrumCuts = params->GetCutString(kMatched, -1, -1, -1);
  }
  else {
    spectrumCuts = params->GetCutString(kAnyMatchingStatus, -1, -1, -1);
  }

  TString hname = Form("h%s_MesonPt_%s_%s_Truth_Coarse", params->GetName(), spectrumCuts.Data(), matching);

  return hname;
}

//____________________________________________________________________________________
TString DJetCorrAnalysis::GetDPtMeasuredName(Int_t p, const char* matching)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";

  TString spectrumCuts(params->GetCutString(kAnyMatchingStatus, -1, -1, -1));
  TString hname(Form("h%s_DPtSpectrum_%s_%s", params->GetName(), spectrumCuts.Data(), matching));

  return hname;
}

//____________________________________________________________________________________
TString DJetCorrAnalysis::GetDEtaTruthName(Int_t p, const char* matching)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";

  TString spectrumCuts;

  if (strcmp("Matched", matching)==0) {
    spectrumCuts = params->GetCutString(kMatched, -1, -1, -1);
  }
  else {
    spectrumCuts = params->GetCutString(kAnyMatchingStatus, -1, -1, -1);
  }

  TString hname = Form("h%s_MesonEta_%s_%s_Truth", params->GetName(), spectrumCuts.Data(), matching);

  return hname;
}

//____________________________________________________________________________________
TString DJetCorrAnalysis::GetDEtaMeasuredName(Int_t /*p*/, const char* /*matching*/)
{
  return "";
}

//____________________________________________________________________________________
TString DJetCorrAnalysis::GetJetPtTruthName(Int_t p)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";

  TString spectrumCuts(params->GetCutString(kMatched, -1, -1, -1));
  TString hname = Form("h%s_JetPt_%s_Matched_Truth_Coarse", params->GetName(), spectrumCuts.Data());

  return hname;
}

//____________________________________________________________________________________
TString DJetCorrAnalysis::GetJetPtMeasuredName(Int_t p)
{
  DJetCorrAnalysisParams* params = static_cast<DJetCorrAnalysisParams*>(fAnalysisParams->At(p));
  if (!params) return "";

  TString hname;

  return hname;
}
