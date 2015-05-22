#include <TFile.h>
#include <TList.h>
#include <THnSparse.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>

void ProjectResponse(TH2*& resp, TH1*& eff, THnSparse* matching, THnSparse* truth, Double_t minJetPt, Double_t maxJetPt);

void GenerateResponse()
{
  gStyle->SetOptStat(0);
  
  const char* fname = "/Users/saiola/Documents/Work/ALICE/JetResults/Jets_EMC_pp_MC_420/AnalysisResults.root";

  TFile* file = TFile::Open(fname);
  //file->ls();
  //TList* list = (TList*)file->Get("AliJetResponseMaker_Jet_AKTChargedR060_DcandidatesAndTracksDStarMCrec_pT0150_pt_scheme_Jet_AKTChargedR060_mcparticlesDStar_pT0000_pt_scheme_Bias0_BiasType0_TPC_histos");
  TList* list = (TList*)file->Get("AliJetResponseMaker_Jet_AKTFullR060_DcandidatesAndTracksD0MCrec_pT0150_pt_scheme_Jet_AKTFullR060_mcparticlesD0_pT0000_pt_scheme_Bias0_BiasType0_TPC_histos");
  if (!list) {
    Printf("could not get list");
    return;
  }
  THnSparse* fHistMatching = (THnSparse*)list->FindObject("fHistMatching");
  THnSparse* fHistJets2 = (THnSparse*)list->FindObject("fHistJets2");
  
  TH2* resp[2] = {0};
  TH1* eff[2] = {0};

  ProjectResponse(resp[0], eff[0], fHistMatching, fHistJets2, 8, 13);
  ProjectResponse(resp[1], eff[1], fHistMatching, fHistJets2, 13, 50);
}

void ProjectResponse(TH2*& resp, TH1*& eff, THnSparse* matching, THnSparse* truth, Double_t minJetPt, Double_t maxJetPt)
{
  matching->GetAxis(0)->SetRangeUser(minJetPt, maxJetPt);
  resp = matching->Projection(10, 9);
  resp->SetName(Form("resp_%.0f_%.0f", minJetPt, maxJetPt));
  resp->SetTitle(Form("Response matrix %.1f < #it{p}_{T,jet}^{det} < %.1f GeV/#it{c}", minJetPt, maxJetPt));
  resp->GetXaxis()->SetTitle("#it{p}_{D}^{det} / #it{p}_{jet}^{det}");
  resp->GetYaxis()->SetTitle("#it{p}_{D}^{part} / #it{p}_{jet}^{part}");
  resp->Rebin2D(5,5);
  Printf("Resp x bins = %d, y bins = %d", resp->GetNbinsX(), resp->GetNbinsY());
  eff = resp->ProjectionY(Form("eff_%.0f_%.0f", minJetPt, maxJetPt));
  eff->SetTitle(Form("Efficiency %.1f < #it{p}_{T,jet}^{det} < %.1f GeV/#it{c}", minJetPt, maxJetPt));
  eff->GetXaxis()->SetTitle("#it{p}_{D}^{part} / #it{p}_{jet}^{part}");
  eff->GetYaxis()->SetTitle("Efficiency");
  Printf("Eff x bins = %d", eff->GetNbinsX());
  TH1* truth_proj = truth->Projection(4);
  truth_proj->Rebin(5);
  Printf("Truth x bins = %d", truth_proj->GetNbinsX());
  eff->Divide(truth_proj);

  resp->GetXaxis()->SetRangeUser(0.1, 1.01);
  resp->GetYaxis()->SetRangeUser(0.1, 1.01);
  eff->GetXaxis()->SetRangeUser(0.1, 1.01);
  
  TCanvas* c1 = new TCanvas(Form("resp_%.0f_%.0f", minJetPt, maxJetPt),Form("resp_%.0f_%.0f", minJetPt, maxJetPt));
  c1->SetLogz();
  resp->Draw("colz");
  c1->Update();
  c1->SaveAs(Form("%s.pdf",c1->GetName()));

  TCanvas* c2 = new TCanvas(Form("eff_%.0f_%.0f", minJetPt, maxJetPt),Form("resp_%.0f_%.0f", minJetPt, maxJetPt));
  eff->Draw();
  c2->Update();
  c2->SaveAs(Form("%s.pdf",c2->GetName()));
}


