#include <TLorentzVector.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TVector3.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>

#include "DecaySimulator.h"

ClassImp(DecaySimulator);

//____________________________________________________________________________________
DecaySimulator::DecaySimulator() :
  TNamed("DecaySimulator", "DecaySimulator"),
  fMass(1.86484),  // D0
  fDaughterMass1(0.139570),  // Pi+
  fDaughterMass2(0.493677),  // K-
  fNPtBins(0),
  fPtBins(0),
  fMinEta(-0.5),
  fMaxEta(0.5),
  fCMMomentum(0),
  fCMEnergy1(0),
  fCMEnergy2(0),
  fDecaySimulation(0),
  fMother(0),
  fDaughter1(0),
  fDaughter2(0)
{
  fNPtBins = 30;
  fPtBins = new Double_t[fNPtBins];
  for (Int_t i = 0; i < fNPtBins; i++) {
    fPtBins[i] = i;
  }
  fPtBins[0] = 0.5;
}

//____________________________________________________________________________________
DecaySimulator::DecaySimulator(Double_t mass, Double_t massDaughter1, Double_t massDaughter2) :
  TNamed("DecaySimulator", "DecaySimulator"),
  fMass(mass),
  fDaughterMass1(massDaughter1),
  fDaughterMass2(massDaughter2),
  fNPtBins(0),
  fPtBins(0),
  fMinEta(-0.5),
  fMaxEta(0.5),
  fCMMomentum(0),
  fCMEnergy1(0),
  fCMEnergy2(0),
  fDecaySimulation(0),
  fMother(0),
  fDaughter1(0),
  fDaughter2(0)
{
  fNPtBins = 30;
  fPtBins = new Double_t[fNPtBins];
  for (Int_t i = 0; i < fNPtBins; i++) {
    fPtBins[i] = i;
  }
  fPtBins[0] = 0.5;
}

//____________________________________________________________________________________
void DecaySimulator::SetPtBins(Int_t nbins, Double_t *bins)
{
  if (fPtBins) delete[] fPtBins;
  
  fNPtBins = nbins;
  fPtBins = new Double_t[fNPtBins];
  for (Int_t i = 0; i < fNPtBins; i++) {
    fPtBins[i] = bins[i];
  }
}

//____________________________________________________________________________________
void DecaySimulator::Init()
{
  if (gRandom) delete gRandom;
  gRandom = new TRandom3(0);

  fCMEnergy1 = (fMass*fMass + fDaughterMass1*fDaughterMass1 - fDaughterMass2*fDaughterMass2) / (2 * fMass);
  fCMEnergy2 = (fMass*fMass + fDaughterMass2*fDaughterMass2 - fDaughterMass1*fDaughterMass1) / (2 * fMass);

  fCMMomentum = TMath::Sqrt(fCMEnergy1*fCMEnergy1 - fDaughterMass1*fDaughterMass1);

  if (fDecaySimulation) delete fDecaySimulation;

  fDecaySimulation = new TTree("DecaySimulation", "DecaySimulation");
  fDecaySimulation->Branch("Mother", &fMother);
  fDecaySimulation->Branch("Daughter1", &fDaughter1);
  fDecaySimulation->Branch("Daughter2", &fDaughter2);
}

//____________________________________________________________________________________
void DecaySimulator::StartSimulation(Int_t n)
{
  Init();

  for (Int_t ipt = 0; ipt < fNPtBins; ipt++) {
    Printf("Pt = %.3f", fPtBins[ipt]);
    for (Int_t i = 0; i < n; i++) {

      // Generate a random eta-phi for the mother
      Double_t eta = gRandom->Rndm() * (fMaxEta - fMinEta);
      eta += fMinEta;
      Double_t phi = gRandom->Rndm() * TMath::TwoPi();

      fMother->SetPtEtaPhiM(fPtBins[ipt], eta, phi, fMass);

      // Generate a random decay direction
      Double_t phiDecay = gRandom->Rndm() * TMath::TwoPi();
      Double_t thetaDecay = gRandom->Rndm() * TMath::Pi();
      
      CalculateDaughters(thetaDecay, phiDecay);

      fDecaySimulation->Fill();
    }
  }
}

//____________________________________________________________________________________
void DecaySimulator::SaveTree(const char* fileName)
{
  TFile* file = TFile::Open(fileName, "recreate");
  if (!file || file->IsZombie()) {
    Printf("Could not open file '%s'.", fileName);
    return;
  }

  file->cd();
  fDecaySimulation->Write();

  file->Close();
  delete file;
}

//____________________________________________________________________________________
void DecaySimulator::CalculateDaughters(Double_t theta, Double_t phi)
{
  fDaughter1->SetPxPyPzE(fCMMomentum*TMath::Sin(theta)*TMath::Cos(phi),
                         fCMMomentum*TMath::Sin(theta)*TMath::Sin(phi),
                         fCMMomentum*TMath::Cos(theta),
                         fCMEnergy1);

  //fDaughter1->Print();


  fDaughter2->SetPxPyPzE(fCMMomentum*TMath::Sin(TMath::Pi() - theta)*TMath::Cos(TMath::Pi() + phi),
                         fCMMomentum*TMath::Sin(TMath::Pi() - theta)*TMath::Sin(TMath::Pi() + phi),
                         fCMMomentum*TMath::Cos(TMath::Pi() - theta),
                         fCMEnergy2);

  //fDaughter2->Print();
  
  if (fMother->P() > 0) {
    TVector3 boost = fMother->BoostVector();
    
    fDaughter1->Boost(boost);
    fDaughter2->Boost(boost);
  }
}
