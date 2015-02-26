#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TH1F.h>
#include <TF1.h>
#include <TList.h>
#include <Riostream.h>
#include <TString.h>
#include <TRandom3.h>
#include <TMath.h>

TH1* StatisticsSmearing(TH1* orig, const char* hname);
TH1* Fold(TH1* orig, TH1* res, const char* hname);
void PurgeErrors(TH1* hist, Double_t err = 0);
TH1* Fold2D(TH1* orig, TH2* res, const char* hname, TH1* prior);
void SaturationEffects(TH1* hist, Double_t avgNjets, Double_t R, Double_t acc);

void PrepareUnfoldingTest()
{
  Bool_t doSatEffects = kFALSE;

  Double_t bias = 5;
  
  Double_t R = 0.2;
  Double_t acc = (100*TMath::DegToRad() - 2*R) * (1.4 - 2*R);

  Double_t events = 1e7;

  TH1::AddDirectory(kFALSE);

  TString fileName("$JETRESULTS/LEGO_TRAIN_MC_360/MergedResults_Pt2-10.root");
  TString listName("AliJetResponseMaker_Jet_AKTFullR020_PicoTracks_pT0150_CaloClustersCorr_ET0300_pt_scheme_Jet_AKTFullR020_MCParticlesSelected_pT0000_pt_scheme_Bias0_BiasType0_EMCAL_histos");

  TFile *file = TFile::Open(fileName);
  TList *list = static_cast<TList*>(file->Get(listName));
  THnSparse *thn1 = static_cast<THnSparse*>(list->FindObject("fHistJets1"));
  THnSparse *thn2 = static_cast<THnSparse*>(list->FindObject("fHistJets2"));
  THnSparse *thnm = static_cast<THnSparse*>(list->FindObject("fHistMatching"));

  thn1->GetAxis(5)->SetRangeUser(bias,100);
  TH1* histDetectorMCReconstructed = thn1->Projection(2);
  histDetectorMCReconstructed->SetName("histDetectorMCReconstructed");

  thn2->GetAxis(5)->SetRangeUser(bias,1000);
  TH1* histDetectorMCTruth = thn2->Projection(2);
  histDetectorMCTruth->SetName("histDetectorMCTruth");
  
  thnm->GetAxis(7)->SetRangeUser(bias,100);
  TH2* histDetectorResponse = thnm->Projection(1,0);
  histDetectorResponse->SetName("histDetectorResponse");

  delete list;
  list = 0;
  
  file->Close();
  delete file;
  file = 0;

  TString fileNamedpT("$JETRESULTS/LEGO_TRAIN_374_JetType0.root");
  TString hnamedpT("DeltaPtRC_Cent0");
  
  file = TFile::Open(fileNamedpT);
  TH2* hdpt2d = static_cast<TH2*>(file->Get(hnamedpT));
  
  TH1* histBackgroundResponse = hdpt2d->ProjectionY("histBackgroundResponse");

  delete hdpt2d;
  hdpt2d = 0;

  file->Close();
  delete file;
  file = 0;

  TH1* histInclusiveTruth = static_cast<TH1*>(histDetectorMCTruth->Clone("histInclusiveTruth"));
  TH1* histInclusiveInterTruth = 0;

  if (doSatEffects) {
    Double_t avgNjets = 7.86;

    SaturationEffects(histInclusiveTruth, avgNjets, R, acc);
    histInclusiveInterTruth = Fold2D(histInclusiveTruth, histDetectorResponse, "histInclusiveInterTruth", histDetectorMCTruth);
  }
  else {
    histInclusiveInterTruth = static_cast<TH1*>(histDetectorMCReconstructed->Clone("histInclusiveInterTruth"));
  }

  TH1* histInclusiveMeasured_infinite = Fold(histInclusiveInterTruth, histBackgroundResponse, "histInclusiveMeasured_infinite");

  Double_t integ = histInclusiveMeasured_infinite->Integral(histInclusiveMeasured_infinite->GetXaxis()->FindBin(100.), histInclusiveMeasured_infinite->GetXaxis()->FindBin(105.));
  //Double_t scaling = 21. / 1e7 * events / integ;
  Double_t scaling = 1;

  PurgeErrors(histInclusiveTruth);
  PurgeErrors(histInclusiveInterTruth);
  PurgeErrors(histInclusiveMeasured_infinite, 1e-6);

  histInclusiveTruth->Scale(scaling);
  histInclusiveInterTruth->Scale(scaling);
  histInclusiveMeasured_infinite->Scale(scaling);

  TH1* histInclusiveMeasured = StatisticsSmearing(histInclusiveMeasured_infinite, "histInclusiveMeasured");

  TH1* histInclusivePriorMCTruth = histDetectorResponse->ProjectionY("histInclusivePriorMCTruth");
  TH1* histInclusivePriorMCInterTruth = histDetectorResponse->ProjectionX("histInclusivePriorMCInterTruth");
  TH1* histInclusivePriorMeasured = Fold(histInclusivePriorMCInterTruth, histBackgroundResponse, "histInclusivePriorMeasured");

  PurgeErrors(histInclusivePriorMCTruth);
  PurgeErrors(histInclusivePriorMCInterTruth);
  PurgeErrors(histInclusivePriorMeasured);

  histInclusivePriorMCTruth->Scale(scaling);
  histInclusivePriorMCInterTruth->Scale(scaling);
  histInclusivePriorMeasured->Scale(scaling);

  TH1* histEvents = new TH1F("histEvents", "histEvents", 1, 0, 1);
  histEvents->SetBinContent(1, events);
  histEvents->Sumw2();
  histEvents->SetBinError(1, TMath::Sqrt(histEvents->GetBinContent(1)));
  histEvents->SetEntries(histEvents->GetBinContent(1));

  TFile *output = TFile::Open("UnfoldingTestData_bias5_Pt2_10.root", "recreate");
  output->cd();
  histEvents->Write();
  histInclusiveTruth->Write();
  histInclusiveInterTruth->Write();
  
  histInclusiveMeasured_infinite->SetName("histInclusiveMeasured");
  histInclusiveMeasured_infinite->Write();
  //histInclusiveMeasured->Write();

  histInclusivePriorMCTruth->Write();
  histInclusivePriorMCInterTruth->Write();
  histInclusivePriorMeasured->Write();
  histBackgroundResponse->Write();
  histDetectorResponse->Write();
  histDetectorMCTruth->Write();
  histDetectorMCReconstructed->Write();
  output->Close();
  
  delete output;
  output = 0;
}

void PurgeErrors(TH1* hist, Double_t err)
{
  for (Int_t i = 0; i <= hist->GetNbinsX()+1; i++) hist->SetBinError(i, err*hist->GetBinContent(i));
}

TH1* StatisticsSmearing(TH1* orig, const char* hname)
{
  TRandom3 rndm(0);
  
  Int_t n = orig->GetNbinsX();
  Double_t minx = orig->GetXaxis()->GetXmin();
  Double_t maxx = orig->GetXaxis()->GetXmax();
  TH1* result = new TH1F(hname, hname, n, minx, maxx);
  result->Sumw2();

  Double_t entries = 0;

  for (Int_t i = 0; i <= n; i++) {
    Double_t err = TMath::Sqrt(orig->GetBinContent(i));
    Double_t content = rndm.Gaus(orig->GetBinContent(i), err) + 0.5;
    
    if (content < 0) content = 0;
    else content = TMath::Floor(content);

    result->SetBinContent(i, content);
    
    Double_t newerr = TMath::Sqrt(result->GetBinContent(i));
    // if (newerr / content < 1e-5) newerr = 1e-5 * content; // safe check, if errors are too small unfolding will fail!

    result->SetBinError(i, newerr);
    entries += result->GetBinContent(i);
  }

  result->SetEntries(entries);

  return result;
}

TH1* Fold(TH1* orig, TH1* res, const char* hname)
{  
  Double_t norm = res->Integral();
  Int_t n = orig->GetNbinsX();
  Double_t minx = orig->GetXaxis()->GetXmin();
  Double_t maxx = orig->GetXaxis()->GetXmax();
  TH1* result = new TH1F(hname, hname, n*2, -maxx, maxx);
  result->Sumw2();

  for (Int_t i = 1; i <= n; i++) {
    Double_t pt = orig->GetXaxis()->GetBinCenter(i);
    Double_t entries = orig->GetBinContent(i);
    for (Int_t j = 1; j <= res->GetNbinsX(); j++) {
      Double_t pt2 = res->GetXaxis()->GetBinCenter(j);
      Double_t entries2 = res->GetBinContent(j+1);
      result->Fill(pt+pt2, entries*entries2/norm);
    }
  }

  return result;
}

TH1* Fold2D(TH1* orig, TH2* res, const char* hname, TH1* prior)
{
  TH1* result = static_cast<TH1*>(orig->Clone(hname));

  for (Int_t i = 0; i <= result->GetNbinsX()+1; i++) {
    result->SetBinContent(i,0);
    result->SetBinError(i,0);
  }
  
  for (Int_t j = 0; j <= orig->GetNbinsX()+1; j++) {
    Double_t norm = res->Integral(0,res->GetNbinsX()+1,j,j);
    if (orig->GetBinContent(j) < 1e-20) continue;
    if (norm < 1e-20) continue;
    for (Int_t i = 0; i <= result->GetNbinsX()+1; i++) {
      Double_t eff = norm / prior->GetBinContent(j);
      Double_t w = res->GetBinContent(i,j) / norm * eff;
      result->AddBinContent(i, orig->GetBinContent(j) * w);
      result->SetBinError(i, TMath::Sqrt(result->GetBinError(i)*result->GetBinError(i) + orig->GetBinError(j)*orig->GetBinError(j)*w*w));
    }
  }

  result->SetEntries(orig->GetEntries());

  return result;
}

void SaturationEffects(TH1* hist, Double_t avgNjets, Double_t R, Double_t acc)
{
  Double_t Ntot = hist->Integral();
  Double_t events = Ntot /  avgNjets;
  Double_t occupancyRatio = TMath::Pi() * R * R / acc;
  Double_t fractionMergedJets = occupancyRatio * (avgNjets - 1) / 2;

  for (Int_t i = 1; i <= hist->GetNbinsX(); i++) {
    if (hist->GetBinContent(i) <= 0) continue;
    Double_t nup = hist->Integral(i, hist->GetNbinsX());
    Double_t removeThisBin = fractionMergedJets * hist->GetBinContent(i) / Ntot * nup / 2;
    Printf("Bin %d: removing %.2f %%", i, removeThisBin / hist->GetBinContent(i) * 100);
    hist->SetBinContent(i, hist->GetBinContent(i) - removeThisBin);
  }
}
