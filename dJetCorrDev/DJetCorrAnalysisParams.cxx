// This class holds the analysis parameters of the Dmeson - jet correlation analysis
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TString.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TF1.h>

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
  fTracksName(),
  fInputListName(),
  fInvMinMass(0),
  fInvMaxMass(0),
  f2ProngMinMass(0),
  f2ProngMaxMass(0), 
  fDeltaInvMinMass(0),
  fDeltaInvMaxMass(0),
  fMinDEta(-0.9),
  fMaxDEta(0.9),
  fMinJetConstituents(2),
  fMassFitTypeSig(MassFitter::kGaus),
  fMassFitTypeBkg(MassFitter::kExpo)
{
}

//____________________________________________________________________________________
DJetCorrAnalysisParams::DJetCorrAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius, const char* tracksName) :
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
  fTracksName(tracksName),
  fInputListName(),
  fInvMinMass(0),
  fInvMaxMass(0),
  f2ProngMinMass(0),
  f2ProngMaxMass(0), 
  fDeltaInvMinMass(0),
  fDeltaInvMaxMass(0),
  fMinDEta(-0.9),
  fMaxDEta(0.9),
  fMinJetConstituents(2),
  fMassFitTypeSig(MassFitter::kGaus),
  fMassFitTypeBkg(MassFitter::kExpo)
{
  fInputListName = Form("AliAnalysisTaskDmesonJetCorrelations_%s_rec_Jet_AKT%s%s_%s_pT0150_pt_scheme_TPC_histos", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data(), fTracksName.Data());
  fName = Form("%s_%s_%s", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data());

  if (fDmesonName == "DStar") {
    fNDPtBins = 8;
    fDPtBins = new Double_t[fNDPtBins+1];
    fDPtBins[ 0] =  1.0;
    fDPtBins[ 1] =  3.0;
    fDPtBins[ 2] =  5.0;
    fDPtBins[ 3] =  6.0;
    fDPtBins[ 4] =  7.0;
    fDPtBins[ 5] =  9.0;
    fDPtBins[ 6] = 12.0;
    fDPtBins[ 7] = 16.0;
    fDPtBins[ 8] = 30.0;

    SetDeltaInvMassRange(413, 421, 0.08);
    SetInvMassRange(413, 0.60);
    Set2ProngMassRange(421, 0.20);

    fMassFitTypeSig = MassFitter::kGaus;
    fMassFitTypeBkg = MassFitter::kExpoPower;
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

    fMassFitTypeSig = MassFitter::kGaus;
    fMassFitTypeBkg = MassFitter::kExpo;
  }

  fNJetPtBins = 3;
  fJetPtBins = new Double_t[fNJetPtBins+1];
  fJetPtBins[ 0] =   2.0;
  fJetPtBins[ 1] =   8.0;
  fJetPtBins[ 2] =  13.0;
  fJetPtBins[ 3] =  50.0;

  fNzBins = 5;
  fzBins = new Double_t[fNzBins+1];
  //fzBins[ 0] =   0.10;
  fzBins[ 0] =   0.30;
  fzBins[ 1] =   0.50;
  fzBins[ 2] =   0.80;
  fzBins[ 3] =   1.00;
  fzBins[ 4] =   1.20;
  fzBins[ 5] =   2.00;
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
  fDeltaInvMaxMass(p.fDeltaInvMaxMass),
  fMinDEta(p.fMinDEta),
  fMaxDEta(p.fMaxDEta),
  fMinJetConstituents(p.fMinJetConstituents),
  fMassFitTypeSig(p.fMassFitTypeSig),
  fMassFitTypeBkg(p.fMassFitTypeBkg)
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

//____________________________________________________________________________________
TString DJetCorrAnalysisParams::GetCutString(Int_t st, Int_t dptBin, Int_t jetptBin, Int_t dzBin)
{
  Double_t minDPt = GetMinDPt();
  Double_t maxDPt = GetMaxDPt();
  if (dptBin >= 0) {
    minDPt = GetDPtBin(dptBin);
    maxDPt = GetDPtBin(dptBin+1);
  }

  Double_t minJetPt = GetMinJetPt();
  Double_t maxJetPt = GetMaxJetPt();
  if (jetptBin >= 0) {
    minJetPt = GetJetPtBin(jetptBin);
    maxJetPt = GetJetPtBin(jetptBin+1);
  }

  Double_t minz = GetMinZ();
  Double_t maxz = GetMaxZ();
  if (dzBin >= 0) {
    minz = GetzBin(dzBin);
    maxz = GetzBin(dzBin+1);
  }

  TString dCuts(Form("DPt_%02.0f_%02.0f", minDPt, maxDPt));
  dCuts.ReplaceAll(".", "");

  TString jetCuts(Form("JetPt_%03.0f_%03.0f", minJetPt, maxJetPt));
  jetCuts.ReplaceAll(".", "");
  
  TString zCuts(Form("z_%.1f_%.1f", minz, maxz));
  zCuts.ReplaceAll(".", "");
    
  TString cuts(Form("%s_%s_%s", dCuts.Data(), jetCuts.Data(), zCuts.Data()));

  if (st != kMatched) {
    // Ignore cuts on z and jet pt
    cuts = dCuts;
  }

  return cuts;
}

//____________________________________________________________________________________
MassFitter* DJetCorrAnalysisParams::CreateMassFitter(const char* name) const
{
  Double_t minMass = fInvMinMass;
  Double_t maxMass = fInvMaxMass;
  Double_t startingSigma = 0.01;
  Double_t startingSigmaBkg = 0.01;
  if (fDmesonName == "DStar") {
    minMass = fDeltaInvMinMass;
    maxMass = fDeltaInvMaxMass;
    startingSigma = 7e-4;
    startingSigmaBkg = 5;
  }
  
  MassFitter* fitter = new MassFitter(name, fMassFitTypeSig, fMassFitTypeBkg, minMass, maxMass);
  fitter->GetFitFunction()->SetParameter(1, startingSigmaBkg);
  //fitter->GetFitFunction()->SetParLimits(1, startingSigmaBkg*0.1, startingSigmaBkg*10.);
  fitter->GetFitFunction()->SetParameter(4, startingSigma);
  fitter->GetFitFunction()->SetParLimits(4, startingSigma*0.5, startingSigma*1.5);
  
  return fitter;
}
