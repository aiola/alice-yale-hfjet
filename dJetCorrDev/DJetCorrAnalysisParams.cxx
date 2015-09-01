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
  fTruthInputListName(),
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
DJetCorrAnalysisParams::DJetCorrAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius,
                                               const char* tracksName, DJetCorrAnalysisType anaType, Bool_t isMC, Bool_t isBkgSub) :
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
  fTruthInputListName(),
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
  if (anaType == kInvMassAna) {
    if (isMC) {
      if (isBkgSub) {
        fInputListName = Form("AliAnalysisTaskDmesonJetCorrelations_%s_MCrec_Jet_AKT%s%s_%sMCrec_pT0150_pt_scheme_TPC_histos", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data(), fTracksName.Data());
      }
      else {
        fInputListName = Form("AliAnalysisTaskDmesonJetCorrelations_%s_Jet_AKT%s%s_%s_pT0150_pt_scheme_TPC_histos", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data(), fTracksName.Data());
      }
      fTruthInputListName = Form("AliAnalysisTaskDmesonJetCorrelations_%s_MC_Jet_AKT%s%s_mcparticles%s_pT0000_pt_scheme_TPC_histos", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data(), fDmesonName.Data());
    }
    else {
      fInputListName = Form("AliAnalysisTaskDmesonJetCorrelations_%s_rec_Jet_AKT%s%s_%s_pT0150_pt_scheme_TPC_histos", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data(), fTracksName.Data());
    }
  }
  else if (anaType == kResponseMatrixAna) {
    fInputListName = Form("AliJetResponseMaker_Jet_AKT%s%s_%s_pT0150_pt_scheme_Jet_AKT%s%s_mcparticles%s_pT0000_pt_scheme_Bias0_BiasType0_TPC_histos", fJetType.Data(), fJetRadius.Data(), fTracksName.Data(), fJetType.Data(), fJetRadius.Data(), fDmesonName.Data());
  }
  fName = Form("%s_%s_%s", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data());

  if (fDmesonName == "DStar") {
    fNDPtBins = 6;
    fDPtBins = new Double_t[fNDPtBins+1];
    fDPtBins[ 0] =  1.6;
    fDPtBins[ 1] =  4.0;
    fDPtBins[ 2] =  6.0;
    fDPtBins[ 3] =  8.0;
    fDPtBins[ 4] = 10.0;
    fDPtBins[ 5] = 15.0;
    fDPtBins[ 6] = 30.0;

    SetDeltaInvMassRange(413, 421, 0.08);
    SetInvMassRange(413, 0.60);
    Set2ProngMassRange(421, 0.20);

    fMassFitTypeSig = MassFitter::kGaus;
    fMassFitTypeBkg = MassFitter::kExpoPower;

    fNJetPtBins = 2;
    fJetPtBins = new Double_t[fNJetPtBins+1];
    fJetPtBins[ 0] =   8.0;
    fJetPtBins[ 1] =  13.0;
    fJetPtBins[ 2] =  30.0;
  }
  else {
    fNDPtBins = 6;
    fDPtBins = new Double_t[fNDPtBins+1];
    fDPtBins[ 0] =  1.0;
    fDPtBins[ 1] =  3.0;
    fDPtBins[ 2] =  5.0;
    fDPtBins[ 3] =  7.0;
    fDPtBins[ 4] = 10.0;
    fDPtBins[ 5] = 15.0;
    fDPtBins[ 6] = 30.0;

    SetInvMassRange(421, 0.30);

    fMassFitTypeSig = MassFitter::kGaus;
    fMassFitTypeBkg = MassFitter::kExpo;

    fNJetPtBins = 3;
    fJetPtBins = new Double_t[fNJetPtBins+1];
    fJetPtBins[ 0] =   5.0;
    fJetPtBins[ 1] =   8.0;
    fJetPtBins[ 2] =  13.0;
    fJetPtBins[ 3] =  30.0;
  }

  fNzBins = 4;
  fzBins = new Double_t[fNzBins+1];
  fzBins[ 0] =   0.20;
  fzBins[ 1] =   0.40;
  fzBins[ 2] =   0.60;
  fzBins[ 3] =   0.80;
  fzBins[ 4] =   1.00;
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
  fitter->GetFitFunction()->SetParameter(4, startingSigma);
  
  return fitter;
}

//____________________________________________________________________________________
void DJetCorrAnalysisParams::GetDPtBinRange(Double_t& minDPt, Double_t& maxDPt, Int_t dptBin) const
{
  if (dptBin >= 0 && dptBin < GetNDPtBins()) {
    minDPt = GetDPtBin(dptBin);
    maxDPt = GetDPtBin(dptBin+1);
  }
  else {
    minDPt = GetMinDPt();
    maxDPt = GetMaxDPt();
  }
}

//____________________________________________________________________________________
void DJetCorrAnalysisParams::GetzBinRange(Double_t& minZ, Double_t& maxZ, Int_t zBin) const
{
  if (zBin >= 0 && zBin < GetNzBins()) {
    minZ = GetzBin(zBin);
    maxZ = GetzBin(zBin+1);
  }
  else {
    minZ = GetMinZ();
    maxZ = GetMaxZ();
  }
}

//____________________________________________________________________________________
void DJetCorrAnalysisParams::GetJetPtBinRange(Double_t& minJetPt, Double_t& maxJetPt, Int_t jetPtBin) const
{
  if (jetPtBin >= 0 && jetPtBin < GetNJetPtBins()) {
    minJetPt = GetJetPtBin(jetPtBin);
    maxJetPt = GetJetPtBin(jetPtBin+1);
  }
  else {
    minJetPt = GetMinJetPt();
    maxJetPt = GetMaxJetPt();
  }
}
