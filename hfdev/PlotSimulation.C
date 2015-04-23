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
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLegendEntry.h>

TTree* CreateChain(const char* fname);
void PlotHist(TH1* hist, Style_t marker, Double_t markerSize, Color_t color, TLegend* leg = 0);

void PlotSimulation()
{  
  TTree* tree = CreateChain("simulation.root");
  if (!tree) {
    Printf("Unable to load tree.");
    return;
  }

  TLorentzVector* mother = 0;
  TLorentzVector* daughter1 = 0;
  TLorentzVector* daughter2 = 0;

  tree->SetBranchAddress("Mother", &mother);
  tree->SetBranchAddress("Daughter1", &daughter1);
  tree->SetBranchAddress("Daughter2", &daughter2);

  TH1* decayAngle1[30] = {0};
  TH1* decayAngle2[30] = {0};

  for (Int_t i = 0; i < 30; i++) {
    Printf("Creating histograms for pt %d", i);
    Double_t pt = i;
    if (i == 0) pt = 0.5;
    
    TString hname1(Form("decayAngle1_Pt%d", i));
    //TString htitle1(Form("K^{#pm} w/ #it{p}_{T,D^{0}} = %.1f GeV/#it{c}", pt));
    TString htitle1(Form("#it{p}_{T,D^{0}} = %.1f GeV/#it{c}", pt));
    decayAngle1[i] = new TH1F(hname1, htitle1, 50, 0, TMath::Pi());
    
    TString hname2(Form("decayAngle2_Pt%d", i));
    TString htitle2(Form("#pi^{#pm} w/ #it{p}_{T,D^{0}} = %.1f GeV/#it{c}", pt));
    decayAngle2[i] = new TH1F(hname2, htitle2, 50, 0, TMath::Pi());
  }

  Int_t nentries = tree->GetEntries();

  for (Int_t i = 0; i < nentries; i++) {
    if (i % 10000 == 0) Printf("Analyzing event %d", i);
    
    tree->GetEntry(i);

    if (!mother) {
      Printf("Could not retrive mother particle for event %d", i);
      continue;
    }

    if (!daughter1) {
      Printf("Could not retrive daughter 1 particle for event %d", i);
      continue;
    }

    if (!daughter2) {
      Printf("Could not retrive daughter 2 particle for event %d", i);
      continue;
    }

    if (daughter1->Pt() < 0.15) continue;
    if (daughter2->Pt() < 0.15) continue;

    Int_t ptbin = TMath::FloorNint(mother->Pt());
    if (ptbin >= 30) {
      Printf("Particle %d has pt = %.3f > 30 GeV/c.", i, mother->Pt());
      continue;
    }

    if (ptbin < 0) {
      Printf("Particle %d has negative pt = %.3f (?).", i, mother->Pt());
      continue;
    }

    Double_t a1 = 0;mother->Angle(daughter1->Vect());
    Double_t a2 = 0;mother->Angle(daughter2->Vect());

    if (daughter1->M() > 0.15) { // kaon
      a1 = mother->Angle(daughter1->Vect());
      a2 = mother->Angle(daughter2->Vect());
    }
    else {
      a1 = mother->Angle(daughter2->Vect());
      a2 = mother->Angle(daughter1->Vect());
    }
    
    decayAngle1[ptbin]->Fill(a1);
    decayAngle2[ptbin]->Fill(a2);
  }

  for (Int_t i = 0; i < 30; i++) {
    decayAngle1[i]->Sumw2();
    decayAngle2[i]->Sumw2();
    
    decayAngle1[i]->Scale(1. / decayAngle1[i]->Integral(), "width");
    decayAngle2[i]->Scale(1. / decayAngle2[i]->Integral(), "width");
  }

  TCanvas* canvas = new TCanvas("decayAngle", "decayAngle", 600, 600);
  canvas->cd();

  TH1* myBlankHist = new TH1F("myBlankHist", "myBlankHist", 100, 0, TMath::Pi());
  myBlankHist->GetXaxis()->SetTitle("angle_{lab} (daughter, mother)");
  myBlankHist->GetYaxis()->SetTitle("Probability density");
  myBlankHist->GetYaxis()->SetRangeUser(0,5);
  myBlankHist->Draw();

  TLegend* leg = new TLegend(0.28, 0.35, 0.68, 0.88);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(43);
  leg->SetTextSize(15);

  PlotHist(decayAngle1[1], kFullCircle, 0.8, kBlack, leg);
  PlotHist(decayAngle2[1], kOpenCircle, 0.8, kBlack, 0);
  
  PlotHist(decayAngle1[3], kFullSquare, 0.8, kRed+1, leg);
  PlotHist(decayAngle2[3], kOpenSquare, 0.8, kRed+1, 0);

  PlotHist(decayAngle1[5], kFullTriangleUp, 0.8, kBlue+1, leg);
  PlotHist(decayAngle2[5], kOpenTriangleUp, 0.8, kBlue+1, 0);

  PlotHist(decayAngle1[7], kFullTriangleDown, 0.8, kGreen+2, leg);
  PlotHist(decayAngle2[7], kOpenTriangleDown, 0.8, kGreen+2, 0);

  PlotHist(decayAngle1[12], kFullDiamond, 1.2, kMagenta+1, leg);
  PlotHist(decayAngle2[12], kOpenDiamond, 1.2, kMagenta+1, 0);

  leg->Draw();

  TLegend* leg2 = new TLegend(0.54, 0.22, 0.94, 0.30);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(15);
  leg2->AddEntry((TObject*)0, "Full symbols: K^{#pm}", "");
  leg2->AddEntry((TObject*)0, "Open symbols: #pi^{#pm}", "");
  leg2->Draw();
}

void PlotHist(TH1* hist, Style_t marker, Double_t markerSize, Color_t color, TLegend* leg)
{
  hist->SetMarkerStyle(marker);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize(markerSize);
  hist->SetLineColor(color);
  hist->Draw("same");
  if (leg) {
    TLegendEntry* entry = 0;

    entry = leg->AddEntry(hist, Form("%s, #mu = %.2f, #sigma = %.2f", hist->GetTitle(), hist->GetMean(), hist->GetRMS()), "pe");
    entry->SetTextFont(43);
    entry->SetTextSize(14);

    Double_t int06 = hist->Integral(0, hist->GetXaxis()->FindBin(0.6), "width") * 100;
    Double_t int04 = hist->Integral(0, hist->GetXaxis()->FindBin(0.4), "width") * 100;

    entry = leg->AddEntry((TObject*)0, Form("P(angle < 0.6) = %.1f %%, P(angle < 0.4) = %.1f %%", int06, int04), "");
    entry->SetTextFont(43);
    entry->SetTextSize(14);
  }
}
    

TTree* CreateChain(const char* fname)
{
  TFile* file = TFile::Open(fname);

  if (!file || file->IsZombie()) {
    Printf("Unable to open file %s", fname);
    return 0;
  }
  
  TTree* tree = dynamic_cast<TTree*>(file->Get("DecaySimulation"));
  return tree;
}
