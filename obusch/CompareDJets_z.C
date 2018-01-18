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

int   Style[]   = {kFullSquare, kFullSquare};
int   Color[]   = {kRed, kRed};
float Size[]    = {1.2, 1.2};

const char* input_path = "/Volumes/DATA/ALICE/JetResults";

TList* getList(TString name);
Bool_t isConsistentDouble(Double_t x, Double_t y);
TGraphErrors* histoToGraph(TH1* hist);
TGraphErrors* histoToGraphSys(TH1* hist, TGraphErrors* grSysTot);
Bool_t GetDJetZDistribution(TString label, TH1*& stat, TGraphErrors*& syst);
Bool_t GetInclusiveJetCrossSection(TH1*& stat, TGraphErrors*& syst);
TGraphErrors* RebinGraph(TGraphErrors* input, TAxis* axis);
TH1* RebinHistogram(TH1* input, TAxis* axis);
TGraphErrors* DivideGraphs(TGraphErrors* num, TGraphErrors* den);
TCanvas* Plot(Int_t incl_bin, TGraphErrors* incl_stat, TGraphErrors* incl_syst, Double_t minZ, Double_t maxZ, Double_t minY, Double_t maxY, TH1* zStat, TGraphErrors* zSyst);

void CompareDJets_z()
{
  TH1::AddDirectory(kFALSE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TGraphErrors* InclSyst = nullptr;
  TH1* InclStat = nullptr;

  TGraphErrors* DJetLowPtSyst = nullptr;
  TH1* DJetLowPtStat = nullptr;

  TGraphErrors* DJetHighPtSyst = nullptr;
  TH1* DJetHighPtStat = nullptr;

  Bool_t r = kFALSE;

  r = GetInclusiveJetCrossSection(InclStat, InclSyst);
  if (!r) exit(1);

  r = GetDJetZDistribution("5_15", DJetLowPtStat, DJetLowPtSyst);
  if (!r) exit(1);

  r = GetDJetZDistribution("15_30", DJetHighPtStat, DJetHighPtSyst);
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
  Double_t jet_pt_lim[] = {5, 15, 30};
  TAxis* axis_jet_pt = new TAxis(2, jet_pt_lim);
  TH1* hXSec = RebinHistogram(InclStat, axis_jet_pt);

  TGraphErrors* grXSec = histoToGraph(hXSec);
  grXSec->SetName("grXSec");

  TGraphErrors* grXSecSys_orig = histoToGraphSys(InclStat, InclSyst);
  TGraphErrors* grXSecSys = RebinGraph(grXSecSys_orig, axis_jet_pt);
  grXSecSys->SetName("grXSecSys");

  // -------  
  // plot
  TCanvas* canvas_low_pt = Plot(0, grXSec, grXSecSys, 0.2, 1.0, 0, 0.15, DJetLowPtStat, DJetLowPtSyst);

  TCanvas* canvas_high_pt = Plot(1, grXSec, grXSecSys, 0.4, 1.0, 0, 0.15, DJetHighPtStat, DJetHighPtSyst);

  if (doWrite) {
    TString strTit_low = "R_momfrac_lowpt.pdf";
    canvas_low_pt->SaveAs(strTit_low);

    TString strTit_high = "R_momfrac_highpt.pdf";
    canvas_high_pt->SaveAs(strTit_high);
  }
}

TCanvas* Plot(Int_t incl_bin, TGraphErrors* incl_stat, TGraphErrors* incl_syst, Double_t minZ, Double_t maxZ, Double_t minY, Double_t maxY, TH1* zStat, TGraphErrors* zSyst)
{
  TString cname = TString::Format("JetPt_%d_comparison", int(incl_stat->GetX()[incl_bin]));
  TCanvas *canvas = new TCanvas(cname, cname,0,0,500,500);
  canvas->SetLeftMargin(0.15);
  canvas->SetBottomMargin(0.13);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  for (int i = 0; i < zStat->GetNbinsX(); i++) {
    Double_t v = zStat->GetBinContent(i + 1) / incl_stat->GetY()[incl_bin] / incl_stat->GetErrorX(incl_bin) / 2;
    Double_t e_stat = TMath::Sqrt(zStat->GetBinError(i + 1)*zStat->GetBinError(i + 1) / zStat->GetBinContent(i + 1) / zStat->GetBinContent(i + 1) +
        incl_stat->GetErrorY(incl_bin) * incl_stat->GetErrorY(incl_bin) / incl_stat->GetY()[incl_bin] / incl_stat->GetY()[incl_bin]) * v;
    Double_t e_syst = TMath::Sqrt(zSyst->GetErrorY(i)*zSyst->GetErrorY(i) / zStat->GetBinContent(i + 1) / zStat->GetBinContent(i + 1) +
        incl_syst->GetErrorY(incl_bin) * incl_syst->GetErrorY(incl_bin) / incl_stat->GetY()[incl_bin] / incl_stat->GetY()[incl_bin]) * v;
    zStat->SetBinContent(i+1, v);
    zStat->SetBinError(i+1, e_stat);
    zSyst->GetY()[i] = v;
    zSyst->GetEYhigh()[i] = e_syst;
    zSyst->GetEYlow()[i] = e_syst;
  }

  zSyst->GetXaxis()->SetNdivisions(int((maxZ-minZ)*20), 5, 0);

  zSyst->GetXaxis()->SetLabelFont(43);
  zSyst->GetXaxis()->SetLabelSize(20);
  zSyst->GetXaxis()->SetTitleFont(43);
  zSyst->GetXaxis()->SetTitleSize(22);
  zSyst->GetXaxis()->SetTitle("#it{z}_{||, ch jet}");

  zSyst->GetYaxis()->SetTickLength(0.02);
  zSyst->GetYaxis()->SetLabelOffset(0);
  zSyst->GetYaxis()->SetLabelFont(43);
  zSyst->GetYaxis()->SetLabelSize(20);
  zSyst->GetYaxis()->SetTitleOffset(1.5);
  zSyst->GetYaxis()->SetTitleFont(43);
  zSyst->GetYaxis()->SetTitleSize(22);
  zSyst->GetYaxis()->CenterTitle(false);
  zSyst->GetYaxis()->SetTitle("#it{R}(#it{p}_{T}, #it{z}_{||}) / #Delta #it{z}_{||}");

  zSyst->SetFillColor(kGray);
  zSyst->GetXaxis()->SetRangeUser(minZ, maxZ);
  zSyst->GetYaxis()->SetRangeUser(minY, maxY);

  zSyst->SetMarkerColor(Color[incl_bin]);
  zSyst->SetLineColor(Color[incl_bin]);
  zSyst->SetMarkerStyle(Style[incl_bin]);
  zSyst->SetMarkerSize(Size[incl_bin]);
  zSyst->Draw("AP2");

  zStat->SetMarkerColor(Color[incl_bin]);
  zStat->SetLineColor(Color[incl_bin]);
  zStat->SetMarkerStyle(Style[incl_bin]);
  zStat->SetMarkerSize(Size[incl_bin]);
  zStat->Draw("same");

  /*
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
  entry->SetLineColor(Color[incl_bin]);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(Color[incl_bin]);
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
*/
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

  return canvas;
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

    gr->SetPoint(bin-1,center,cont);
    gr->SetPointError(bin-1, 0.5*width,errSys);
  }

  return gr;
}

Bool_t GetDJetZDistribution(TString label, TH1*& stat, TGraphErrors*& syst)
{
  TString fname = TString::Format("%s/JetZCrossSection_JetPt_%s_Systematics.root", input_path, label.Data());
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
