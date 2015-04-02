#include <TLorentzVector.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <Riostream.h>
#include <TVector3.h>

Double_t CalculateDecayAngle(Double_t* x, Double_t* par)
{
  TLorentzVector mother(x[0]*TMath::Cos(x[1]), 0, x[0]*TMath::Sin(x[1]), TMath::Sqrt(x[0]*x[0] + par[0]*par[0]));
  
  Double_t p2 = (par[0]*par[0] - par[1]*par[1] - par[2]*par[2]) / 2; 
  Double_t p = TMath::Sqrt(p2);

  TLorentzVector decay1(0, 0, p, TMath::Sqrt(p2 + par[1]*par[1]));
  TLorentzVector decay2(0, 0, -p, TMath::Sqrt(p2 + par[2]*par[2]));

  TVector3 boost = mother.BoostVector();

  //boost.Print();
  
  decay1.Boost(boost);
  decay2.Boost(boost);

  //decay1.Print();
  //decay2.Print();

  Double_t angle = decay1.Angle(mother.Vect());

  return angle;
}

void DecayAngle()
{
  Double_t masses[] = {1.865, 0.139, 0.494};

  TF2* graph = new TF2("graph", CalculateDecayAngle, 0, 40, 0, TMath::TwoPi(), 3);
  graph->SetParameters(masses);
  graph->Draw("colz");

  //Printf("Decay angle = %.3f", graph->Eval(2, TMath::Pi()*3/2));
}
