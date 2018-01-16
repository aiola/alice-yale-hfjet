#if !defined( __CINT__) || defined(__MAKECINT__)

#include "TSystem.h"
#include "TColor.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TVector2.h"
#include "TVector3.h"
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
#include "TRandom.h"
#include "TROOT.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TH3F.h"
#include "TMultiGraph.h"
#include "TLatex.h"

#include <iostream>

#endif

using namespace std;

Bool_t doWrite = kTRUE;

int gray  = TColor::GetColor("#cccccc");
int red   = TColor::GetColor("#FF0000");
int blue  = TColor::GetColor("#0000FF");
int pink  = TColor::GetColor("#FF00FF");
int black = TColor::GetColor("#000000");
int green = TColor::GetColor("#006600");
int other = TColor::GetColor("#776600");

int solidSquare      = 21;
int openSquare       = 25;
int openCircle       = 24;
int openTriangleUp   = 26;
int openTriangleDown = 32;
int openDiamond      = 27;
int openStar         = 30;

int   AliceStyle   = solidSquare;
int   AliceColor   = red;
float AliceSize    = 1.5;

int   P11Style     = openCircle;
int   P11Color     = blue;
float P11Size      = 1.4;

int   MonashStyle  = openDiamond;
int   MonashColor  = green;
float MonashSize    = 1.4;

const Double_t lowLimPlot = 5; // lowest jet pt plotted - don't fill graphs below, to avoid plotting artifacts

TString fnameMC = "/Users/obusch1/work/FFprep/miniJets/PythiaFastJet/tunes/PythiaFastJet_tunes_spec_MinB.root";

// ---------------------------------------------------

TList* getList(TString name = ""){

  TList* list = 0x0; 

    
  TString dirName  = "PWGJE_FragmentationFunction_" + name;
  TString listName = "fracfunc_" + name;    

  gDirectory->cd(dirName);

  list = (TList*) gDirectory->Get(listName);

  return list;
}

// ------------------------------------------------------------------------------ 

Bool_t isConsistentDouble(Double_t x, Double_t y){

  if(0.99999*x < y && 1.00001*x>y && 0.99999*y<x && 1.00001*y>x) return kTRUE;
  return kFALSE;
}

// --------------------------------------------------------

TGraphErrors* histoToGraph(TH1* hist){

  Int_t nBins = hist->GetNbinsX();
  TGraphErrors* gr = new TGraphErrors(nBins);

  for(Int_t bin=1; bin<=nBins; bin++){ 
    
    Double_t center = hist->GetBinCenter(bin);
    Double_t width  = hist->GetBinWidth(bin);
    Double_t cont   = hist->GetBinContent(bin);
    Double_t err    = hist->GetBinError(bin);

    if(center < lowLimPlot){
      cont = 0;
      err  = 0;
    }
    
    gr->SetPoint(bin-1,center,cont);
    gr->SetPointError(bin-1, 0.5*width,err);
  }

  return gr;
}

// --------------------------------------------------------

TGraphErrors* histoToGraphSys(TH1* hist, TGraphErrors* grSysTot, Bool_t doRatioGraph = kFALSE){

  // convert histo to graph, with errors from grSysTot
  // doRatioGraph: set datapoint to 1, error = rel syst. error
  
  if(hist->GetNbinsX() != grSysTot->GetN()){
    cout<<" histoToGraphSys: inconsistent nBins histo / graph "<<endl;
    exit(0);
  }

  Int_t nBins = hist->GetNbinsX();
  TGraphErrors* gr = new TGraphErrors(nBins);
  
  for(Int_t bin=1; bin<=nBins; bin++){ 

    Double_t center     = hist->GetBinCenter(bin);
    Double_t width      = hist->GetBinWidth(bin);
    Double_t cont       = hist->GetBinContent(bin);
    Double_t errStat    = hist->GetBinError(bin);
    Double_t relErrStat = cont ? errStat/cont : 0; 
    
    Double_t sysCent, sysCont;
    grSysTot->GetPoint(bin-1,sysCent,sysCont);

    if(!isConsistentDouble(sysCent,center)){
      cout<<" histoToGraphSys: inconsistent bin center histo / graph "<<endl;
      exit(0);
    }

    Double_t relErrSys = grSysTot->GetErrorY(bin-1);
    Double_t errSys    = relErrSys * cont;

    //cout<<" bin "<<bin<<" center "<<center<<" relErrSys "<<relErrSys<<" cont "<<cont<<" errSys "<<errSys<<endl;

    if(center < lowLimPlot){
      cont     = 0;
      errSys   = 0;
    }

    if(doRatioGraph){
      //Double_t ratioErr = TMath::Sqrt(relErrSys*relErrSys + relErrStat*relErrStat);
      Double_t ratioErr = relErrSys;
      if(center < lowLimPlot) ratioErr = 0;

      gr->SetPoint(bin-1,center,1);
      gr->SetPointError(bin-1,0.5*width,ratioErr);
    }
    else{
      gr->SetPoint(bin-1,center,cont);
      gr->SetPointError(bin-1, 0.5*width,errSys);
    }
  }

  return gr;
}

// --------------------------------------------------------

TGraphErrors* readMC(TString strGener, Bool_t doUESub = kFALSE){

  TFile f(fnameMC,"READ");

  gDirectory->cd(strGener);

  TH1D* hSpec;
  if(doUESub) hSpec = (TH1D*) gDirectory->Get("hSpecUESub");
  else        hSpec = (TH1D*) gDirectory->Get("hSpec");

  hSpec->SetDirectory(0);

  f.Close();
  
  TGraphErrors* gr = histoToGraph(hSpec);
  gr->SetName(Form("grSpec%s",strGener.Data()));
  
  delete hSpec;

  return gr;
}

// ---------------------------------------------------------

TGraphErrors* ratioDataMC(TGraphErrors* grData,TGraphErrors* grMC){

  // ratio MC / data
  // old: use only stat. errors of MC (not data)
  // new: use stat. error of data + MC
  
  Int_t nPointsData = grData->GetN();
  Int_t nPointsMC   = grMC->GetN();
  Int_t nPoints     = TMath::Min(nPointsData,nPointsMC); // MC histos go up to 150 GeV, data to 200

  TGraphErrors* grRatio = new TGraphErrors(nPoints);
  grRatio->SetName(Form("grRatio%s",grMC->GetName()));
  
  for(Int_t p=0; p<nPoints; p++){ 

    Double_t dataCent, dataCont, dataErrX, dataErrY;
    grData->GetPoint(p,dataCent,dataCont);
    dataErrX = grData->GetErrorX(p);
    dataErrY = grData->GetErrorY(p);

    Double_t MCCent, MCCont, MCErrX, MCErrY;
    grMC->GetPoint(p,MCCent,MCCont);
    MCErrX = grMC->GetErrorX(p);
    MCErrY = grMC->GetErrorY(p);

    if(!isConsistentDouble(dataCent,MCCent)){
      cout<<" ratioDataMC: inconsistent bin center histo / graph "<<endl;
      cout<<" dataCent "<<dataCent<<" MCCent "<<MCCent<<endl;
      exit(0);
    }
    if(!isConsistentDouble(dataErrX,MCErrX)){
      cout<<" ratioDataMC: inconsistent bin width histo / graph "<<endl;
      exit(0);
    }

    Double_t ratio = 0;
    if(dataCont) ratio = MCCont / dataCont;
    
    Double_t errRatio = 0;
    if(ratio){
      Double_t errRatioMC   = MCErrY / MCCont * ratio;
      Double_t errRatioData = dataErrY / dataCont * ratio;
      errRatio = TMath::Sqrt(errRatioMC*errRatioMC + errRatioData*errRatioData);
    }
      
    if(dataCent < lowLimPlot){
      ratio   = 0;
      errRatio = 0;
    }

    grRatio->SetPoint(p,dataCent,ratio);
    grRatio->SetPointError(p,MCErrX,errRatio);
  }

  return grRatio;
  
}

// --------------------------------------------------------

void compSpecMC(Bool_t doUESub = kFALSE){
   
  // read unfolded spectra

  //TString fname = "../unfoldSpec_2016_02_27/outData_spec_Bayes.root";
  //if(doUESub) fname = "../unfoldSpec_2016_02_27/outData_spec_Bayes_UESub.root";

  //TString fname = "../unfoldSpec_lego141/outData_spec_Bayes.root";
  //if(doUESub) fname = "../unfoldSpec_lego141/outData_spec_Bayes_UESub.root";

  //TString fname = "../unfoldSpec_lego141_MinBPtHMC/outData_spec_Bayes_combPtH.root";
  //if(doUESub) fname = "../unfoldSpec_lego141_MinBPtHMC/outData_spec_Bayes_UESub_combPtH.root";

  TString fname = "outData_spec_Bayes_combPtH.root";
  if(doUESub) fname = "../unfoldSpec_lego182_MinBPtHMC/outData_spec_Bayes_UESub_combPtH.root";

  
  TString fnameSys = "outSys_tot_spec.root";
  
  // read data spectra

  TFile f(fname,"READ");
  
  TH1D* hSpec =  (TH1D*) f.Get("hSpecComb");
  if(!hSpec) hSpec = (TH1D*) f.Get("hSpecUnfolded_iter4");

  hSpec->SetDirectory(0);
  
  f.Close();

  // read sys errors
  
  TFile g(fnameSys,"READ");
  TGraphErrors *grSysTot;
  if(doUESub) grSysTot = (TGraphErrors*) g.Get("grSysErrTot");
  else        grSysTot = (TGraphErrors*) g.Get("grSysErrTotNoUE");
  g.Close();

  //cout<<" gr "<<grSysTot<<endl;
  // 
  
  // normalize
  // spectra come already normalized per event and per bin

  TH1D* hXSec = (TH1D*) hSpec->Clone("hXSec");

  Double_t sigmapp      = 73.2;        // published x-section in mbarn
  Double_t triggerEff   = 62.2 / 73.2; // eff of MB OR trigger to be applied to data
  Double_t bin0CorrData = 0.91;        // ratio event after vertexNcontributors cut to evts after phys selection

  hXSec->Scale(sigmapp*triggerEff);
  hXSec->Scale(bin0CorrData);

  // spectra histos to TGraph
  
  TGraphErrors* grXSec = histoToGraph(hXSec);
  grXSec->SetName("grXSec");

  TGraphErrors* grXSecSys = histoToGraphSys(hXSec,grSysTot);
  grXSecSys->SetName("grXSecSys");
  
  TGraphErrors* grRatioData = histoToGraphSys(hXSec,grSysTot,kTRUE);
  grXSecSys->SetName("grRatioData");

  // read MC 
  //TGraphErrors* grXSecP2011  = readMC("Perugia2011",doUESub);
  //TGraphErrors* grXSecMonash = readMC("Monash",doUESub);

  //TGraphErrors* grRatioP2011  = ratioDataMC(grXSec,grXSecP2011);
  //TGraphErrors* grRatioMonash = ratioDataMC(grXSec,grXSecMonash);
  
  // -------  
  // plot

  TCanvas *mcComparison = new TCanvas("mcComparison", "mcComparison",0,0,620,700);
  mcComparison->Range(0,0,1,1);
  mcComparison->SetFillColor(0);
  mcComparison->SetBorderMode(0);
  mcComparison->SetBorderSize(2);
  mcComparison->SetTickx(1);
  mcComparison->SetTicky(1);
  mcComparison->SetLeftMargin(0);
  mcComparison->SetRightMargin(0);
  mcComparison->SetTopMargin(0);
  mcComparison->SetBottomMargin(0);
  mcComparison->SetFrameBorderMode(0);

  TPad *pad01 = new TPad("pad01", "pad01",0,0.35,1,1);
  pad01->Draw();
  pad01->cd();
  //pad01->Range(-2.024097,-6.522879,104,-1.44508);
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

  grXSecSys->SetMarkerColor(AliceColor);
  grXSecSys->SetLineColor(AliceColor);
  grXSecSys->SetMarkerStyle(AliceStyle);
  grXSecSys->SetMarkerSize( AliceSize);

  grXSecSys->GetXaxis()->SetLabelFont(42);
  grXSecSys->GetXaxis()->SetLabelSize(0.03);
  grXSecSys->GetXaxis()->SetTitleSize(0.03);
  grXSecSys->GetXaxis()->SetTitleFont(42);
  grXSecSys->GetXaxis()->SetTitleOffset(1.1);
  grXSecSys->GetYaxis()->SetTitle("{p}_{T}^{ch jet}");
  
  grXSecSys->GetYaxis()->SetLabelFont(42);
  grXSecSys->GetYaxis()->SetTickLength(0.02);
  grXSecSys->GetYaxis()->SetLabelOffset(0);
  grXSecSys->GetYaxis()->SetLabelSize(0.06);
  grXSecSys->GetYaxis()->SetTitleSize(0.07);
  grXSecSys->GetYaxis()->SetTitleOffset(0.9);
  grXSecSys->GetYaxis()->SetTitleFont(42);
  grXSecSys->GetYaxis()->CenterTitle(false);
  grXSecSys->GetYaxis()->SetTitle("d^{2}#sigma^{ch jet}/d#it{p}_{T}d#eta (mb #it{c}/GeV)");

  grXSecSys->SetFillColor(gray);
  grXSecSys->SetMinimum(2e-07);
  grXSecSys->SetMaximum(5e00);
  grXSecSys->GetXaxis()->SetRangeUser(lowLimPlot,99);

  grXSecSys->Draw("AP2");

  grXSec->SetMarkerColor(AliceColor);
  grXSec->SetLineColor(AliceColor);
  grXSec->SetMarkerStyle(AliceStyle);
  grXSec->SetMarkerSize( AliceSize);
  grXSec->Draw("P");

  // check
  if(1){
    double checkCent, checkCont,checkErrStat,checkErrSys;
    int checkPoint = 25;

    grXSec->GetPoint(checkPoint,checkCent,checkCont);
    checkErrStat = grXSec->GetErrorY(checkPoint);
    checkErrSys  = grXSecSys->GetErrorY(checkPoint);
    Double_t checkRelErrSys  = checkErrSys/checkCont;
    Double_t checkRelErrStat = checkErrStat/checkCont;
    cout<<" pt "<<checkCent<<" sigma "<<checkCont<<" err stat "<<checkErrStat<<" sys "<<checkErrSys
  <<" relErrStat "<<checkRelErrStat<<" sys "<<checkRelErrSys<<endl;
  }
  
  
  /*
  grXSecP2011->SetMarkerColor(P11Color);
  grXSecP2011->SetLineColor(P11Color);
  grXSecP2011->SetMarkerStyle(P11Style);
  grXSecP2011->SetMarkerSize(P11Size);
  grXSecP2011->Draw("P");
  

  grXSecMonash->SetMarkerColor(MonashColor);
  grXSecMonash->SetLineColor(MonashColor);
  grXSecMonash->SetMarkerStyle(MonashStyle);
  grXSecMonash->SetMarkerSize(MonashSize);
  grXSecMonash->Draw("P");
*/
  
  TLegend *leg = new TLegend(0.46,0.62,0.95,0.9,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.055);
  leg->SetTextFont(42);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
    
  TLegendEntry *entry;

  entry=leg->AddEntry("Alice","ALICE","p");
  entry->SetLineColor(AliceColor);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(AliceColor);
  entry->SetLineColor(AliceColor);
  entry->SetMarkerStyle(AliceStyle);
  entry->SetMarkerSize(AliceSize);
    
  
  entry=leg->AddEntry("P11","PYTHIA Perugia 2011","p");
  entry->SetLineColor(P11Color);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(P11Color);
  entry->SetLineColor(P11Color);
  entry->SetMarkerStyle(P11Style);
  entry->SetMarkerSize(P11Size);
  

  entry=leg->AddEntry("Monash","PYTHIA8 Monash","p");
  entry->SetLineColor(MonashColor);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(MonashColor);
  entry->SetLineColor(MonashColor);
  entry->SetMarkerStyle(MonashStyle);
  entry->SetMarkerSize(1.4*MonashSize); // draw diamond a bit larger

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

  // TLatex* tex = new TLatex(10,1.4e-06,"|#it{#eta}^{jet}| < 0.5");
  // tex->SetTextFont(42);
  // tex->SetTextSize(0.06);
  // tex->SetLineWidth(2);
  // tex->Draw();

  
  //pad01->Modified();


  // ratio

  mcComparison->cd();
  
  TPad* pad4 = new TPad("pad4", "pad4",0,0,1,0.36);
  pad4->Draw();
  pad4->cd();
  //pad4->Range(-2.024097,-0.9153846,104,2.7);
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


  grRatioData->SetTitle("");

  grRatioData->GetXaxis()->SetTitle("#it{p}_{T}^{ch jet} (GeV/#it{c})");
  grRatioData->GetXaxis()->CenterTitle(false);
  grRatioData->GetXaxis()->SetLabelFont(42);
  grRatioData->GetXaxis()->SetLabelOffset(0);
  grRatioData->GetXaxis()->SetLabelSize(0.12);
  grRatioData->GetXaxis()->SetTitleSize(0.12);
  grRatioData->GetXaxis()->SetTitleFont(42);
  grRatioData->GetXaxis()->SetTitleOffset(1.1);
  
  grRatioData->GetYaxis()->SetTitle("MC/data ");
  grRatioData->GetYaxis()->CenterTitle(true);
  grRatioData->GetYaxis()->SetTickLength(0.02);
  grRatioData->GetYaxis()->SetNdivisions(604);
  grRatioData->GetYaxis()->SetLabelFont(42);
  grRatioData->GetYaxis()->SetLabelSize(0.10);
  grRatioData->GetYaxis()->SetTitleSize(0.12);
  grRatioData->GetYaxis()->SetTitleOffset(0.6);
  grRatioData->GetYaxis()->SetTitleFont(42);

  grRatioData->SetFillColor(gray);
  grRatioData->SetMinimum(0.4);
  grRatioData->SetMaximum(2.2);
  grRatioData->GetXaxis()->SetRangeUser(lowLimPlot,99);
  grRatioData->GetYaxis()->SetRangeUser(0.4,1.6);

  grRatioData->Draw("APM2");
  
  TLine *line = new TLine(5,1,99,1);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();
/*
  grRatioP2011->SetMarkerColor(P11Color);
  grRatioP2011->SetLineColor(P11Color);
  grRatioP2011->SetMarkerStyle(P11Style);
  grRatioP2011->SetMarkerSize(P11Size);
  grRatioP2011->Draw("P");


  grRatioMonash->SetMarkerColor(MonashColor);
  grRatioMonash->SetLineColor(MonashColor);
  grRatioMonash->SetMarkerStyle(MonashStyle);
  grRatioMonash->SetMarkerSize(MonashSize);
  grRatioMonash->Draw("P");
*/
  gPad->RedrawAxis();
  
  
  if(doWrite){
    TString strTit     = "jetXsec_R04_compMC.pdf";
    if(doUESub) strTit = "jetXsec_R04_compMC_UEsub.pdf";
    mcComparison->SaveAs(strTit);
  }  

}
