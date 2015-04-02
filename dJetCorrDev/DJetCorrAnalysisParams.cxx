// This class holds the analysis parameters of the Dmeson - jet correlation analysis
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TString.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

#include "DJetCorrAnalysisParams.h"

ClassImp(DJetCorrAnalysisParams);

//____________________________________________________________________________________
DJetCorrAnalysisParams::DJetCorrAnalysisParams() :
  TObject(),
  fName(),
  fNDPtBins(0),
  fDPtBins(0),
  fNJetPtBins(0),
  fJetPtBins(0),
  fNzBins(0),
  fzBins(0),
  fJetType(),
  fJetRadius(),
  fDmesonName(),
  fInputListName(),
  fInvMinMass(0),
  fInvMaxMass(0),
  f2ProngMinMass(0),
  f2ProngMaxMass(0), 
  fDeltaInvMinMass(0),
  fDeltaInvMaxMass(0)
{
}

//____________________________________________________________________________________
DJetCorrAnalysisParams::DJetCorrAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius) :
  TObject(),
  fName(),
  fNDPtBins(0),
  fDPtBins(0),
  fNJetPtBins(0),
  fJetPtBins(0),
  fNzBins(0),
  fzBins(0),
  fJetType(jetType),
  fJetRadius(jetRadius),
  fDmesonName(dmeson),
  fInputListName(),
  fInvMinMass(0),
  fInvMaxMass(0),
  f2ProngMinMass(0),
  f2ProngMaxMass(0), 
  fDeltaInvMinMass(0),
  fDeltaInvMaxMass(0)
{
  fInputListName = Form("AliAnalysisTaskDmesonJetCorrelations_%s_Jet_AKT%s%s_AODFilterTracks_pT0150_pt_scheme_TPC_histos", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data());
  fName = Form("%s_%s_%s", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data());

  if (fDmesonName == "DStar") {
    fNDPtBins = 11;
    fDPtBins = new Double_t[fNDPtBins+1];
    fDPtBins[ 0] =  2.5;
    fDPtBins[ 1] =  3.0;
    fDPtBins[ 2] =  4.0;
    fDPtBins[ 3] =  5.0;
    fDPtBins[ 4] =  6.0;
    fDPtBins[ 5] =  7.0;
    fDPtBins[ 6] =  8.0;
    fDPtBins[ 7] = 10.0;
    fDPtBins[ 8] = 15.0;
    fDPtBins[ 9] = 20.0;
    fDPtBins[10] = 30.0;
    fDPtBins[11] = 50.0;

    SetDeltaInvMassRange(413, 421, 0.08);
    SetInvMassRange(413, 0.15);
    Set2ProngMassRange(421, 0.30);
  }
  else {
    fNDPtBins = 9;
    fDPtBins = new Double_t[fNDPtBins+1];
    fDPtBins[ 0] =  0.5;
    fDPtBins[ 1] =  1.0;
    fDPtBins[ 2] =  2.0;
    fDPtBins[ 3] =  3.0;
    fDPtBins[ 4] =  4.0;
    fDPtBins[ 5] =  5.0;
    fDPtBins[ 6] =  6.0;
    fDPtBins[ 7] = 10.0;
    fDPtBins[ 8] = 20.0;
    fDPtBins[ 9] = 40.0;

    SetInvMassRange(421, 0.30);
  }

  fNJetPtBins = 9;
  fJetPtBins = new Double_t[fNJetPtBins+1];
  fJetPtBins[ 0] =   1.0;
  fJetPtBins[ 1] =   5.0;
  fJetPtBins[ 2] =  10.0;
  fJetPtBins[ 3] =  15.0;
  fJetPtBins[ 4] =  20.0;
  fJetPtBins[ 5] =  25.0;
  fJetPtBins[ 6] =  30.0;
  fJetPtBins[ 7] =  40.0;
  fJetPtBins[ 8] =  60.0;
  fJetPtBins[ 9] =  80.0;
  fJetPtBins[10] = 100.0;

  fNzBins = 6;
  fzBins = new Double_t[fNzBins+1];
  fzBins[ 0] =   0.0;
  fzBins[ 1] =   0.4;
  fzBins[ 2] =   0.6;
  fzBins[ 3] =   0.8;
  fzBins[ 4] =   1.0;
  fzBins[ 5] =   1.2;
  fzBins[ 6] =   2.0;
}

//____________________________________________________________________________________
DJetCorrAnalysisParams::DJetCorrAnalysisParams(const DJetCorrAnalysisParams& p) :
  TObject(p),
  fName(p.fName),
  fNDPtBins(p.fNDPtBins),
  fDPtBins(0),
  fNJetPtBins(p.fNJetPtBins),
  fJetPtBins(0),
  fNzBins(p.fNzBins),
  fzBins(0),
  fJetType(p.fJetType),
  fJetRadius(p.fJetRadius),
  fDmesonName(p.fDmesonName),
  fInputListName(p.fInputListName),
  fInvMinMass(p.fInvMinMass),
  fInvMaxMass(p.fInvMaxMass),
  f2ProngMinMass(p.f2ProngMinMass),
  f2ProngMaxMass(p.f2ProngMaxMass),
  fDeltaInvMinMass(p.fDeltaInvMinMass),
  fDeltaInvMaxMass(p.fDeltaInvMaxMass)
{
  fDPtBins = new Double_t[fNDPtBins+1];
  for (Int_t i = 0; i<= fNDPtBins; i++) fDPtBins[i] = p.fDPtBins[i];
  
  fJetPtBins = new Double_t[fNJetPtBins+1];
  for (Int_t i = 0; i<= fNJetPtBins; i++) fJetPtBins[i] = p.fJetPtBins[i];
  
  fzBins = new Double_t[fNzBins+1];
  for (Int_t i = 0; i<= fNzBins; i++) fzBins[i] = p.fzBins[i];
}

//____________________________________________________________________________________
void DJetCorrAnalysisParams::SetInvMassRange(Int_t pdg, Double_t range)
{
  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg));
  fInvMinMass = part->Mass() - range/2;
  fInvMaxMass = part->Mass() + range/2;
}

//____________________________________________________________________________________
void DJetCorrAnalysisParams::Set2ProngMassRange(Int_t pdg, Double_t range)
{
  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg));
  f2ProngMinMass = part->Mass() - range/2;
  f2ProngMaxMass = part->Mass() + range/2;
}

//____________________________________________________________________________________
void DJetCorrAnalysisParams::SetDeltaInvMassRange(Int_t pdg1, Int_t pdg2, Double_t range)
{
  TParticlePDG* part1 = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg1));
  TParticlePDG* part2 = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg2));
  
  fDeltaInvMinMass = part1->Mass() - part2->Mass() - range/2;
  fDeltaInvMaxMass = part1->Mass() - part2->Mass() + range/2;
}
