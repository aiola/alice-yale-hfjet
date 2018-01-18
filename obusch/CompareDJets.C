#if !defined( __CINT__) || defined(__MAKECINT__)

#include "TSystem.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
#include "TLine.h"
#include "TF1.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TMath.h"
#include "TLatex.h"

#include <iostream>

#endif

Bool_t doWrite = kTRUE;

int   InclusiveStyle   = kFullSquare;
int   InclusiveColor   = kRed;
float InclusiveSize    = 1.2;

int   DJetStyle   = kFullCircle;
int   DJetColor   = kBlue;
float DJetSize    = 1.2;

const char* input_path = "/Volumes/DATA/ALICE/JetResults";

const Double_t lowLimPlot = 5; // lowest jet pt plotted - don't fill graphs below, to avoid plotting artifacts
const Double_t upLimPlot = 30; // highest jet pt plotted - don't fill graphs below, to avoid plotting artifacts

TList* getList(TString name);
Bool_t isConsistentDouble(Double_t x, Double_t y);
TGraphErrors* histoToGraph(TH1* hist);
TGraphErrors* histoToGraphSys(TH1* hist, TGraphErrors* grSysTot);
Bool_t GetDJetCrossSection(TH1*& stat, TGraphErrors*& syst);
Bool_t GetInclusiveJetCrossSection(TH1*& stat, TGraphErrors*& syst);
TGraphErrors* RebinGraph(TGraphErrors* input, TAxis* axis);
TH1* RebinHistogram(TH1* input, TAxis* axis);
TGraphErrors* DivideGraphs(TGraphErrors* num, TGraphErrors* den);

void CompareDJets()
{
  TH1::AddDirectory(kFALSE);

  TGraphErrors* InclSyst = nullptr;
  TGraphErrors* DJetSyst = nullptr;

  TH1* InclStat = nullptr;
  TH1* DJetStat = nullptr;

  Bool_t r = kFALSE;

  r = GetInclusiveJetCrossSection(InclStat, InclSyst);
  if (!r) exit(1);

  r = GetDJetCrossSection(DJetStat, DJetSyst);
  if (!r) exit(1);

  // normalize
  // spectra come already normalized per event and per bin

  Double_t sigmapp      = 73.2;        // published x-section in mbarn
  Double_t triggerEff   = 62.2 / 73.2; // eff of MB OR trigger to be applied to data
  Double_t bin0CorrData = 0.91;        // ratio event after vertexNcontributors cut to evts after phys selection

  InclStat->Scale(sigmapp*triggerEff);
  InclStat->Scale(bin0CorrData);

  // spectra histos to TGraph

  //TH1* hXSec = static_cast<TH1*>(InclStat->Clone("hXSec"));
  TH1* hXSec = RebinHistogram(InclStat, DJetStat->GetXaxis());

  TGraphErrors* grXSec = histoToGraph(hXSec);
  grXSec->SetName("grXSec");


  TGraphErrors* grXSecSys_orig = histoToGraphSys(InclStat, InclSyst);
  TGraphErrors* grXSecSys = RebinGraph(grXSecSys_orig, DJetStat->GetXaxis());
  grXSecSys->SetName("grXSecSys");

  // -------  
  // plot

  TCanvas *comparison = new TCanvas("comparison", "comparison",0,0,620,700);
  comparison->Range(0,0,1,1);
  comparison->SetFillColor(0);
  comparison->SetBorderMode(0);
  comparison->SetBorderSize(2);
  comparison->SetTickx(1);
  comparison->SetTicky(1);
  comparison->SetLeftMargin(0);
  comparison->SetRightMargin(0);
  comparison->SetTopMargin(0);
  comparison->SetBottomMargin(0);
  comparison->SetFrameBorderMode(0);

  TPad *pad01 = new TPad("pad01", "pad01",0,0.35,1,1);
  pad01->Draw();
  pad01->cd();
  pad01->SetFillColor(0);
  pad01->SetBorderMode(0);
  pad01->SetBorderSize(0);
  pad01->SetLogy();
  pad01->SetTickx(1);
  pad01->SetTicky(1);
  pad01->SetLeftMargin(0.17);
  pad01->SetRightMargin(0.05);
  pad01->SetTopMargin(0.05);
  pad01->SetBottomMargin(0); // set to 0 to hide axis
  pad01->SetFrameBorderMode(0);
  pad01->SetFrameBorderMode(0);

  grXSecSys->SetTitle("");

  grXSecSys->GetXaxis()->SetLabelFont(42);
  grXSecSys->GetXaxis()->SetLabelSize(0.03);
  grXSecSys->GetXaxis()->SetTitleSize(0.03);
  grXSecSys->GetXaxis()->SetTitleFont(42);
  grXSecSys->GetXaxis()->SetTitleOffset(1.1);
  grXSecSys->GetYaxis()->SetTitle("{p}_{T, ch jet} (GeV/#it{c})");

  grXSecSys->GetYaxis()->SetLabelFont(42);
  grXSecSys->GetYaxis()->SetTickLength(0.02);
  grXSecSys->GetYaxis()->SetLabelOffset(0);
  grXSecSys->GetYaxis()->SetLabelSize(0.06);
  grXSecSys->GetYaxis()->SetTitleSize(0.07);
  grXSecSys->GetYaxis()->SetTitleOffset(0.9);
  grXSecSys->GetYaxis()->SetTitleFont(42);
  grXSecSys->GetYaxis()->CenterTitle(false);
  grXSecSys->GetYaxis()->SetTitle("d^{2}#sigma^{ch jet}/d#it{p}_{T}d#eta (mb #it{c}/GeV)");

  grXSecSys->SetFillColor(kGray);
  grXSecSys->SetMinimum(2e-07);
  grXSecSys->SetMaximum(5e00);
  grXSecSys->GetXaxis()->SetRangeUser(lowLimPlot,30);

  grXSecSys->SetMarkerColor(InclusiveColor);
  grXSecSys->SetLineColor(InclusiveColor);
  grXSecSys->SetMarkerStyle(InclusiveStyle);
  grXSecSys->SetMarkerSize(InclusiveSize);
  grXSecSys->Draw("AP2");

  grXSec->SetMarkerColor(InclusiveColor);
  grXSec->SetLineColor(InclusiveColor);
  grXSec->SetMarkerStyle(InclusiveStyle);
  grXSec->SetMarkerSize(InclusiveSize);
  grXSec->Draw("P");

  DJetSyst->SetMarkerColor(DJetColor);
  DJetSyst->SetLineColor(DJetColor);
  DJetSyst->SetMarkerStyle(DJetStyle);
  DJetSyst->SetMarkerSize(DJetSize);
  DJetSyst->Draw("P2 same");

  DJetStat->SetMarkerColor(DJetColor);
  DJetStat->SetLineColor(DJetColor);
  DJetStat->SetMarkerStyle(DJetStyle);
  DJetStat->SetMarkerSize(DJetSize);
  DJetStat->Draw("P same");

  TLegend* leg = new TLegend(0.46,0.73,0.95,0.93,nullptr,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(18);
  leg->SetTextFont(43);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  TLegendEntry* entry = nullptr;

  entry=leg->AddEntry("Alice","Inclusive Jets","p");
  entry->SetLineColor(InclusiveColor);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(InclusiveColor);
  entry->SetLineColor(InclusiveColor);
  entry->SetMarkerStyle(InclusiveStyle);
  entry->SetMarkerSize(InclusiveSize);

  entry=leg->AddEntry("Alice","D^{0}-Tagged Jets","p");
  entry->SetLineColor(DJetColor);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(DJetColor);
  entry->SetLineColor(DJetColor);
  entry->SetMarkerStyle(DJetStyle);
  entry->SetMarkerSize(DJetSize);

  leg->Draw();

  TLatex* tex = new TLatex(10,6.0e-06,"pp  #sqrt{#it{s}} = 7 TeV");  
  tex->SetTextFont(42);
  tex->SetTextSize(0.06);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex* tex2 = new TLatex(10,1.4e-06,"anti-#it{k}_{T} #it{R} = 0.4");
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.06);
  tex2->SetLineWidth(2);
  tex2->Draw();

  // ratio

  comparison->cd();

  TPad* pad4 = new TPad("pad4", "pad4",0,0,1,0.36);
  pad4->Draw();
  pad4->cd();
  pad4->SetFillColor(0);
  pad4->SetBorderMode(0);
  pad4->SetBorderSize(2);
  pad4->SetTickx(1);
  pad4->SetTicky(1);  
  pad4->SetLeftMargin(0.17);
  pad4->SetRightMargin(0.05);
  pad4->SetTopMargin(0);
  pad4->SetBottomMargin(0.35);
  pad4->SetFrameBorderMode(0);
  pad4->SetFrameBorderMode(0);

  TGraphErrors* grRatioData = DivideGraphs(DJetSyst, grXSecSys);
  grRatioData->SetTitle("");

  grRatioData->GetXaxis()->SetTitle("#it{p}_{T, ch jet} (GeV/#it{c})");
  grRatioData->GetXaxis()->CenterTitle(false);
  grRatioData->GetXaxis()->SetLabelFont(42);
  grRatioData->GetXaxis()->SetLabelOffset(0);
  grRatioData->GetXaxis()->SetLabelSize(0.12);
  grRatioData->GetXaxis()->SetTitleSize(0.12);
  grRatioData->GetXaxis()->SetTitleFont(42);
  grRatioData->GetXaxis()->SetTitleOffset(1.1);

  grRatioData->GetYaxis()->SetTitle("Ratio");
  grRatioData->GetYaxis()->CenterTitle(true);
  grRatioData->GetYaxis()->SetTickLength(0.02);
  grRatioData->GetYaxis()->SetNdivisions(604);
  grRatioData->GetYaxis()->SetLabelFont(42);
  grRatioData->GetYaxis()->SetLabelSize(0.10);
  grRatioData->GetYaxis()->SetTitleSize(0.12);
  grRatioData->GetYaxis()->SetTitleOffset(0.6);
  grRatioData->GetYaxis()->SetTitleFont(42);

  grRatioData->SetFillColor(kGray);
  grRatioData->GetXaxis()->SetRangeUser(lowLimPlot,upLimPlot);
  grRatioData->GetYaxis()->SetRangeUser(0, 0.09);

  grRatioData->Draw("APM2");

  TH1* ratioStat = new TH1F("ratio", "ratio", DJetStat->GetNbinsX(), DJetStat->GetXaxis()->GetXbins()->GetArray());
  ratioStat->Divide(DJetStat, hXSec);
  ratioStat->SetMarkerColor(DJetColor);
  ratioStat->SetLineColor(DJetColor);
  ratioStat->SetMarkerStyle(DJetStyle);
  ratioStat->SetMarkerSize(DJetSize);
  ratioStat->Draw("same");

  TLine *line = new TLine(5,1,99,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  //line->Draw();

  gPad->RedrawAxis();

  if (doWrite) {
    TString strTit = "jetXsec_R04_compDJets.pdf";
    comparison->SaveAs(strTit);
  }  

}

TList* getList(TString name)
{
  TList* list = 0x0;

  TString dirName  = "PWGJE_FragmentationFunction_" + name;
  TString listName = "fracfunc_" + name;

  gDirectory->cd(dirName);

  list = (TList*) gDirectory->Get(listName);

  return list;
}

Bool_t isConsistentDouble(Double_t x, Double_t y)
{
  if (0.99999 * x < y && 1.00001 * x > y && 0.99999 * y < x && 1.00001 * y > x) return kTRUE;
  return kFALSE;
}

TGraphErrors* histoToGraph(TH1* hist)
{
  Int_t nBins = hist->GetNbinsX();
  TGraphErrors* gr = new TGraphErrors(nBins);

  for(Int_t bin=1; bin<=nBins; bin++){

    Double_t center = hist->GetBinCenter(bin);
    Double_t width  = hist->GetBinWidth(bin);
    Double_t cont   = hist->GetBinContent(bin);
    Double_t err    = hist->GetBinError(bin);

    if (center < lowLimPlot || center > upLimPlot) {
      cont = 0;
      err  = 0;
    }

    gr->SetPoint(bin-1,center,cont);
    gr->SetPointError(bin-1, 0.5*width,err);
  }

  return gr;
}

TGraphErrors* histoToGraphSys(TH1* hist, TGraphErrors* grSysTot)
{
  // convert histo to graph, with errors from grSysTot
  // doRatioGraph: set datapoint to 1, error = rel syst. error

  if (hist->GetNbinsX() != grSysTot->GetN()) {
    std::cout << " histoToGraphSys: inconsistent nBins histo / graph " << std::endl;
    exit(0);
  }

  Int_t nBins = hist->GetNbinsX();
  TGraphErrors* gr = new TGraphErrors(nBins);

  for(Int_t bin = 1; bin <= nBins; bin++) {
    Double_t center     = hist->GetBinCenter(bin);
    Double_t width      = hist->GetBinWidth(bin);
    Double_t cont       = hist->GetBinContent(bin);
    Double_t errStat    = hist->GetBinError(bin);
    Double_t relErrStat = cont ? errStat/cont : 0;

    Double_t sysCent, sysCont;
    grSysTot->GetPoint(bin-1,sysCent,sysCont);

    if(!isConsistentDouble(sysCent,center)){
      std::cout << " histoToGraphSys: inconsistent bin center histo / graph " << std::endl;
      exit(0);
    }

    Double_t relErrSys = grSysTot->GetErrorY(bin-1);
    Double_t errSys    = relErrSys * cont;

    if (center < lowLimPlot || center > upLimPlot) {
      cont     = 0;
      errSys   = 0;
    }

    gr->SetPoint(bin-1,center,cont);
    gr->SetPointError(bin-1, 0.5*width,errSys);
  }

  return gr;
}

Bool_t GetDJetCrossSection(TH1*& stat, TGraphErrors*& syst)
{
  TString fname = TString::Format("%s/JetPtSpectrum_DPt_30_Systematics.root", input_path);
  TFile file(fname, "read");
  if (file.IsZombie()) {
    std::cout << "Could not open file " << fname.Data() << std::endl;
    return kFALSE;
  }
  TList* list = static_cast<TList*>(file.Get("FinalSpectrum"));
  if (!list) {
    std::cout << "Could not find FinalSpectrum" << std::endl;
    return kFALSE;
  }
  stat = static_cast<TH1*>(list->FindObject("CentralPointsStatisticalUncertainty"));
  if (!stat) {
    std::cout << "Cannot get measured cross section with statistical uncertainty!" << std::endl;
    return kFALSE;
  }
  syst = static_cast<TGraphErrors*>(list->FindObject("CentralPointsSystematicUncertainty"));
  if (!syst) {
    std::cout << "Cannot get measured cross section with systematic uncertainty!" << std::endl;
    return kFALSE;
  }

  return kTRUE;
}

Bool_t GetInclusiveJetCrossSection(TH1*& stat, TGraphErrors*& syst)
{
  TString fname = "outData_spec_Bayes_combPtH.root";
  TString fnameSys = "outSys_tot_spec.root";

  // read data spectra

  TFile f(fname,"READ");
  stat =  static_cast<TH1*>(f.Get("hSpecComb"));
  if (!stat) stat = static_cast<TH1*>(f.Get("hSpecUnfolded_iter4"));
  f.Close();

  // read sys errors

  TFile g(fnameSys,"READ");
  syst = static_cast<TGraphErrors*>(g.Get("grSysErrTotNoUE"));
  g.Close();

  return kTRUE;
}

TH1* RebinHistogram(TH1* input, TAxis* axis)
{
  TString name = TString::Format("%s_rebinned", input->GetName());
  TH1* result = new TH1D(name, name, axis->GetNbins(), axis->GetXbins()->GetArray());
  result->Sumw2();

  for (Int_t ibin = 1; ibin <= input->GetNbinsX(); ibin++) {
    Double_t v = input->GetBinContent(ibin) * input->GetBinWidth(ibin);
    Double_t e = input->GetBinError(ibin) * input->GetBinWidth(ibin);
    Int_t ibin_r = result->GetXaxis()->FindBin(input->GetBinCenter(ibin));
    if (ibin_r < 1) continue;
    if (ibin_r > result->GetNbinsX()) break;
    result->SetBinContent(ibin_r, result->GetBinContent(ibin_r) + v);
    result->SetBinError(ibin_r, TMath::Sqrt(result->GetBinError(ibin_r)*result->GetBinError(ibin_r) + e*e));
  }

  result->Scale(1.0, "width");
  return result;
}

TGraphErrors* RebinGraph(TGraphErrors* input, TAxis* axis)
{
  TString name = TString::Format("%s_rebinned", input->GetName());
  TGraphErrors* result = new TGraphErrors(axis->GetNbins());
  result->SetName(name);

  Double_t* input_y_array = input->GetY();
  Double_t* input_y_err_array = input->GetEY();

  Double_t* result_y_array = result->GetY();
  Double_t* result_y_err_array = result->GetEY();

  for (Int_t ipoint = 0; ipoint < input->GetN(); ipoint++) {
    Double_t v = input_y_array[ipoint] * input->GetEX()[ipoint] * 2;
    Double_t e = input_y_err_array[ipoint] * input->GetEX()[ipoint] * 2;
    Int_t ipoint_r = axis->FindBin(input->GetX()[ipoint]) - 1;
    if (ipoint_r < 0) continue;
    if (ipoint_r >= axis->GetNbins()) break;
    result_y_array[ipoint_r] += v;
    Double_t upErr = result_y_err_array[ipoint_r] * result_y_err_array[ipoint_r] + e * e;
    result_y_err_array[ipoint_r] = TMath::Sqrt(upErr);
  }

  for (Int_t ipoint = 0; ipoint < result->GetN(); ipoint++) {
    result->GetX()[ipoint] = axis->GetBinCenter(ipoint+1);
    result->GetEX()[ipoint] = axis->GetBinWidth(ipoint+1) / 2;
    result->GetY()[ipoint] /= axis->GetBinWidth(ipoint+1);
    result->GetEY()[ipoint] /= axis->GetBinWidth(ipoint+1);
  }

  return result;
}

TGraphErrors* DivideGraphs(TGraphErrors* num, TGraphErrors* den)
{
  TGraphErrors* ratio = new TGraphErrors(num->GetN());
  ratio->SetName("ratio");

  for (Int_t ipoint = 0; ipoint < ratio->GetN(); ipoint++) {
    if (den->GetY()[ipoint] <= 0 || num->GetY()[ipoint] <= 0) continue;
    ratio->GetX()[ipoint] = num->GetX()[ipoint];
    ratio->GetEX()[ipoint] = num->GetErrorX(ipoint);
    ratio->GetY()[ipoint] = num->GetY()[ipoint] / den->GetY()[ipoint];
    ratio->GetEY()[ipoint] = TMath::Sqrt(num->GetErrorY(ipoint)*num->GetErrorY(ipoint) / num->GetY()[ipoint]/num->GetY()[ipoint] + den->GetErrorY(ipoint)*den->GetErrorY(ipoint) / den->GetY()[ipoint] /den->GetY()[ipoint]) * ratio->GetY()[ipoint];
  }
  return ratio;
}
