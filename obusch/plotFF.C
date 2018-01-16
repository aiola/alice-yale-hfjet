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

int   MonashStyle  = openStar;
int   MonashColor  = green;
float MonashSize    = 1.4;

// int   AliceStyle   = solidSquare;
// int   AliceColor   = red;
// float AliceSize    = 1.5;

int   P0Style      = openSquare;
int   P0Color      = kCyan; //other;
float P0Size       = 1.4;


static  int myDarkRed     = TColor::GetColor(128,0,0);
static  int myDarkGreen   = TColor::GetColor(0,128,0);
static  int myDarkBlue    = TColor::GetColor(0,0,128);

int FFcolor[] = {myDarkRed, myDarkGreen, myDarkBlue}; 

Double_t jetPtLim[] = {5,10,15,20}; 
const Int_t nBinsJetPt = 3;


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

  for(Int_t bin=1; bin<=hist->GetNbinsX(); bin++){ 
    
    Double_t center = hist->GetBinCenter(bin);
    Double_t width  = hist->GetBinWidth(bin);
    Double_t cont   = hist->GetBinContent(bin);
    Double_t err    = hist->GetBinError(bin);

    gr->SetPoint(bin-1,center,cont);
    gr->SetPointError(bin-1, 0.5*width,err);
  }

  return gr;
}

// --------------------------------------------------------

TGraphErrors* histoToGraphSys(TH1* hist, TGraphErrors* grSysTot, Bool_t doRatioGraph = kFALSE){

  // convert histo to graph, bin cont from hist, errors from grSysTot
  // doRatioGraph flag: set bin cont to 1, error = rel syst. error
  
  if(hist->GetNbinsX() != grSysTot->GetN()){
    cout<<" histoToGraphSys: inconsistent nBins histo / graph "<<endl;
    exit(0);
  }

  Int_t nBins = hist->GetNbinsX();
  TGraphErrors* gr = new TGraphErrors(nBins);
  
  for(Int_t bin=1; bin<=hist->GetNbinsX(); bin++){ 
	
    Double_t center     = hist->GetBinCenter(bin);
    Double_t width      = hist->GetBinWidth(bin);
    Double_t cont       = hist->GetBinContent(bin);
    Double_t errStat    = hist->GetBinError(bin);
    Double_t relErrStat = cont ? errStat/cont : 0; 
    
    Double_t sysCent, sysCont;
    grSysTot->GetPoint(bin-1,sysCent,sysCont);

    if(!isConsistentDouble(sysCent,center)){
      cout<<" histoToGraphSys: inconsistent bin center histo / graph : sysCent "<<sysCent<<" center "<<center<<endl;
      exit(0);
    }

    Double_t relErrSys = grSysTot->GetErrorY(bin-1);
    Double_t errSys    = relErrSys * cont;

    //cout<<" bin "<<bin<<" center "<<center<<" relErrSys "<<relErrSys<<" cont "<<cont<<" errSys "<<errSys<<endl;

    if(doRatioGraph){
      //Double_t ratioErr = TMath::Sqrt(relErrSys*relErrSys + relErrStat*relErrStat);
      Double_t ratioErr = relErrSys;

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

TGraphErrors* readMC(TString strGener, TString fnameMC, Double_t jetPtLo, Double_t jetPtHi){

  TFile f(fnameMC,"READ");

  gDirectory->cd(strGener);

  TH1D* hFF = (TH1D*) gDirectory->Get(Form("fh1FFZGen_%02d_%02d",(int)jetPtLo,(int)jetPtHi));

  hFF->SetDirectory(0);

  f.Close();
  
  TGraphErrors* gr = histoToGraph(hFF);
  gr->SetName(Form("grFF%s_%02d_%02d",strGener.Data(),(int)jetPtLo,(int)jetPtHi));
  
  delete hFF;

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
  
    grRatio->SetPoint(p,dataCent,ratio);
    grRatio->SetPointError(p,MCErrX,errRatio);
  }

  return grRatio;
  
}

// --------------------------------------------------------

void plotFF(){
   
  // read unfolded spectra

  TString fname       = "outData_FFZ.root";
  TString fnameSys = "outSys_tot_FF.root";
  //TString fnameMC     = "/Users/obusch1/work/FFprep/miniJets/PythiaFastJet/tunes/PythiaFastJet_tunes_minB_20161206.root";

  
  TGraphErrors* grFF[nBinsJetPt];
  TGraphErrors* grFFSys[nBinsJetPt];
  TGraphErrors* grRatioData[nBinsJetPt];

  TGraphErrors* grSysTot[nBinsJetPt];
/*
  TGraphErrors* grFFP0[nBinsJetPt]; 
  TGraphErrors* grFFP2011[nBinsJetPt]; 
  TGraphErrors* grFFPNoCR[nBinsJetPt]; 
  TGraphErrors* grFFMonash[nBinsJetPt]; 
  
  TGraphErrors* grRatioP0[nBinsJetPt]; 
  TGraphErrors* grRatioP2011[nBinsJetPt]; 
  TGraphErrors* grRatioPNoCR[nBinsJetPt]; 
  TGraphErrors* grRatioMonash[nBinsJetPt]; 
  */
  // read sys errors
  
  TFile g(fnameSys,"READ");

  for(int i=0; i<nBinsJetPt; i++){

    TString strName(Form("grSysErrTot_%02d_%02d",(Int_t)jetPtLim[i],(Int_t)jetPtLim[i+1]));
    grSysTot[i] = (TGraphErrors*) g.Get(strName);
  }
  
  g.Close();
  
  // read data spectra

  TFile f(fname,"READ");

  for(int i=0; i<nBinsJetPt; i++){

    TString strName(Form("dNdzUnfolded_iter4_%02d_%02d",(Int_t)jetPtLim[i],(Int_t)jetPtLim[i+1]));

    TH1D* hFF = (TH1D*) f.Get(strName);

    grFF[i] = histoToGraph(hFF);

    grFFSys[i] = histoToGraphSys(hFF,grSysTot[i]);
    grFFSys[i]->SetName(Form("grFFSys%d",i));
    
    grRatioData[i] = histoToGraphSys(hFF,grSysTot[i],kTRUE);
    grRatioData[i]->SetName(Form("grRatioData%d",i));
  }
    
  f.Close();

  // ---
  // read MC 

  for(int i=0; i<nBinsJetPt; i++){

    TString strName(Form("dNdzUnfolded_iter4_%02d_%02d",(Int_t)jetPtLim[i],(Int_t)jetPtLim[i+1]));
/*
    grFFP0[i]    = readMC("Perugia0Gen",fnameMC,jetPtLim[i],(Int_t)jetPtLim[i+1]);
    grFFP2011[i] = readMC("Perugia2011Gen",fnameMC,jetPtLim[i],(Int_t)jetPtLim[i+1]);
    grFFPNoCR[i] = readMC("PerugiaNoCRGen",fnameMC,jetPtLim[i],(Int_t)jetPtLim[i+1]);
    grFFMonash[i] = readMC("tuneMonashGen",fnameMC,jetPtLim[i],(Int_t)jetPtLim[i+1]);
    
    grRatioP0[i]     = ratioDataMC(grFF[i],grFFP0[i]);
    grRatioP2011[i]  = ratioDataMC(grFF[i],grFFP2011[i]);
    grRatioPNoCR[i]  = ratioDataMC(grFF[i],grFFPNoCR[i]);
    grRatioMonash[i] = ratioDataMC(grFF[i],grFFMonash[i]);
    */
  }

  
  // -------  
  // plot

  TCanvas *cFF = new TCanvas("cFF", "cFF",0,0,620,700);
  cFF->Range(0,0,1,1);
  cFF->SetFillColor(0);
  cFF->SetBorderMode(0);
  cFF->SetBorderSize(2);
  cFF->SetTickx(1);
  cFF->SetTicky(1);
  cFF->SetLeftMargin(0);
  cFF->SetRightMargin(0);
  cFF->SetTopMargin(0);
  cFF->SetBottomMargin(0);
  cFF->SetFrameBorderMode(0);

  TPad *pad01 = new TPad("pad01", "pad01",0,0,1,1);
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
  pad01->SetBottomMargin(0.1); // set to 0 to hide axis
  pad01->SetFrameBorderMode(0);
  pad01->SetFrameBorderMode(0);


  for(Int_t bin=0; bin<nBinsJetPt; bin++){

    grFFSys[bin]->SetTitle("");
    
    grFFSys[bin]->SetMarkerColor(FFcolor[bin]);
    grFFSys[bin]->SetLineColor(FFcolor[bin]);
    grFFSys[bin]->SetMarkerStyle(AliceStyle);
    grFFSys[bin]->SetMarkerSize( AliceSize);
    
    grFFSys[bin]->GetXaxis()->SetLabelFont(42);
    grFFSys[bin]->GetXaxis()->SetLabelSize(0.05);
    grFFSys[bin]->GetXaxis()->SetTitleSize(0.05);
    grFFSys[bin]->GetXaxis()->SetTitleOffset(0.9);
    grFFSys[bin]->GetXaxis()->SetTitleFont(42);
    grFFSys[bin]->GetXaxis()->CenterTitle(false);
    grFFSys[bin]->GetXaxis()->SetTitle("z^{ch}");

    grFFSys[bin]->GetYaxis()->SetLabelFont(42);
    grFFSys[bin]->GetYaxis()->SetLabelOffset(0);
    grFFSys[bin]->GetYaxis()->SetLabelSize(0.05);
    grFFSys[bin]->GetYaxis()->SetTitleSize(0.05);
    grFFSys[bin]->GetYaxis()->SetTitleOffset(1.1);
    grFFSys[bin]->GetYaxis()->SetTitleFont(42);
    grFFSys[bin]->GetYaxis()->CenterTitle(false);
    grFFSys[bin]->GetYaxis()->SetTitle("1/N_{jets} dN/dz^{ch}");
    
    grFFSys[bin]->SetFillColor(gray);
    //grFFSys[bin]->SetMinimum(2e-07);
    //grFFSys[bin]->SetMaximum(5e00);
    grFFSys[bin]->GetXaxis()->SetRangeUser(0,1.06);
    grFFSys[bin]->GetYaxis()->SetRangeUser(9e-02,4.8e01);
    
    if(bin ==0) grFFSys[bin]->Draw("AP2");
    else        grFFSys[bin]->Draw("P2");
    
    grFF[bin]->SetMarkerColor(FFcolor[bin]);
    grFF[bin]->SetLineColor(FFcolor[bin]);
    grFF[bin]->SetMarkerStyle(AliceStyle);
    grFF[bin]->SetMarkerSize( AliceSize);
    grFF[bin]->Draw("P");
  }

  TLegend *leg = new TLegend(0.51,0.74,0.70,0.92,NULL,"brNDC");

  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->SetTextFont(42);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
    
  leg->AddEntry(grFF[0],"#it{p}_{T}^{ch jet}  5 - 10 GeV/#it{c}","P");
  leg->AddEntry(grFF[1],"#it{p}_{T}^{ch jet} 10 - 15 GeV/#it{c}","P");
  leg->AddEntry(grFF[2],"#it{p}_{T}^{ch jet} 15 - 20 GeV/#it{c}","P");
  leg->Draw();


  TLatex* tex = new TLatex(0.08,0.35,"ALICE");
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  tex->Draw();

  TLatex* tex2 = new TLatex(0.08,0.24,"pp  #sqrt{#it{s}} = 7 TeV");
  tex2->SetTextFont(42);
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  tex2->Draw();

  TLatex* tex3 = new TLatex(0.08,0.16,"anti-#it{k}_{T} #it{R} = 0.4");
  tex3->SetTextFont(42);
  tex3->SetTextSize(0.04);
  tex3->SetLineWidth(2);
  tex3->Draw();

  
  // ------------------------------
  // ratios

  TCanvas *cRatio = new TCanvas("cRatio", "cRatio",0,0,620,700);
  cRatio->Range(0,0,1,1);
  cRatio->SetFillColor(0);
  cRatio->SetBorderMode(0);
  cRatio->SetBorderSize(2);
  cRatio->SetTickx(1);
  cRatio->SetTicky(1);
  cRatio->SetLeftMargin(0);
  cRatio->SetRightMargin(0);
  cRatio->SetTopMargin(0);
  cRatio->SetBottomMargin(0);
  cRatio->SetFrameBorderMode(0);

  TPad* padR[nBinsJetPt];
  TLatex *TLJetBinR[nBinsJetPt];

  Double_t padCoord[] = {0,0.40,0.7,1};
  
  for(Int_t bin=0; bin<nBinsJetPt; bin++){

    cRatio->cd();
    
    //padR[bin] = new TPad(Form("padR%d",bin),"",0,bin*0.33,1,(bin+1)*0.33);
    if(bin ==1) padR[bin] = new TPad(Form("padR%d",bin),"",0,padCoord[bin]-0.01,1,padCoord[bin+1]);

    padR[bin] = new TPad(Form("padR%d",bin),"",0,padCoord[bin],1,padCoord[bin+1]);
    padR[bin]->Draw();
    padR[bin]->cd();
    padR[bin]->SetFillColor(0);
    padR[bin]->SetBorderMode(0);
    padR[bin]->SetBorderSize(0);
    padR[bin]->SetTickx(1);
    padR[bin]->SetTicky(1);
    padR[bin]->SetLeftMargin(0.17);
    padR[bin]->SetRightMargin(0.05);

    padR[bin]->SetTopMargin(0);
    padR[bin]->SetBottomMargin(0); // set to 0 to hide axis
    if(bin == 2) padR[bin]->SetTopMargin(0.05);
    if(bin == 0) padR[bin]->SetBottomMargin(0.25); 

    cout<<" bin "<<bin<<" xlowndc "<<padR[bin]->GetAbsXlowNDC()<<endl;
    padR[bin]->SetFrameBorderMode(0);

    grRatioData[bin]->SetTitle("");

    if(bin == 0) grRatioData[bin]->GetXaxis()->SetTitle("z^{ch}");
    grRatioData[bin]->GetXaxis()->CenterTitle(false);
    grRatioData[bin]->GetXaxis()->SetLabelFont(42);
    grRatioData[bin]->GetXaxis()->SetLabelOffset(0);
    grRatioData[bin]->GetXaxis()->SetLabelSize(0.10);
    grRatioData[bin]->GetXaxis()->SetTitleSize(0.12);
    grRatioData[bin]->GetXaxis()->SetTitleFont(42);
    grRatioData[bin]->GetXaxis()->SetRangeUser(0,1.0);
    grRatioData[bin]->GetXaxis()->SetTitleOffset(0.8);

    
    grRatioData[bin]->GetYaxis()->SetTitle("");
    if(bin == 2) grRatioData[bin]->GetYaxis()->SetTitle("MC/data");
    grRatioData[bin]->GetYaxis()->CenterTitle(false);
    grRatioData[bin]->GetYaxis()->SetLabelFont(43); // 42->43: size will be in pixels
    grRatioData[bin]->GetYaxis()->SetLabelSize(25);
    grRatioData[bin]->GetXaxis()->SetTitleFont(42); 
    grRatioData[bin]->GetYaxis()->SetTitleSize(0.14);
    grRatioData[bin]->GetYaxis()->SetTickLength(0.01);
      
    grRatioData[bin]->GetYaxis()->SetTitleOffset(0.4);
    grRatioData[bin]->SetFillColor(gray);
    grRatioData[bin]->GetYaxis()->SetNdivisions(406);
    //grRatioData[bin]->GetYaxis()->SetRangeUser(0.61,1.39);
     grRatioData[bin]->GetYaxis()->SetRangeUser(0.51,1.29);

    grRatioData[bin]->Draw("APM2");

    // Draw horizontal line at unity
    TLine *line = new TLine(0,1,1.0,1);
    line->SetLineStyle(2);
    line->Draw();
/*
    grRatioPNoCR[bin]->SetMarkerColor(P0Color);
    grRatioPNoCR[bin]->SetLineColor(P0Color);
    grRatioPNoCR[bin]->SetMarkerStyle(P0Style);
    grRatioPNoCR[bin]->SetMarkerSize(P0Size);
    grRatioPNoCR[bin]->Draw("P");

    grRatioP2011[bin]->SetMarkerColor(P11Color);
    grRatioP2011[bin]->SetLineColor(P11Color);
    grRatioP2011[bin]->SetMarkerStyle(P11Style);
    grRatioP2011[bin]->SetMarkerSize(P11Size);
    grRatioP2011[bin]->Draw("P");
  
    grRatioMonash[bin]->SetMarkerColor(MonashColor);
    grRatioMonash[bin]->SetLineColor(MonashColor);
    grRatioMonash[bin]->SetMarkerStyle(MonashStyle);
    grRatioMonash[bin]->SetMarkerSize(MonashSize);
    grRatioMonash[bin]->Draw("P");
*/

    // jet pt bins
    TString strTitleTLR(Form("#it{p}_{T}^{ch jet} %d - %d GeV/#it{c}",(Int_t)jetPtLim[bin],(Int_t)jetPtLim[bin+1]));

    if(bin==0)      TLJetBinR[bin] = new TLatex(0.64,0.88,strTitleTLR);
    else if(bin==1) TLJetBinR[bin] = new TLatex(0.64,0.83,strTitleTLR);
    else            TLJetBinR[bin] = new TLatex(0.64,0.78,strTitleTLR);

    TLJetBinR[bin]->SetNDC();   
    TLJetBinR[bin]->SetTextFont(43);
    TLJetBinR[bin]->SetTextSize(22);
    TLJetBinR[bin]->Draw();
  }

 // lower pad 

  padR[0]->cd();

  TLatex* texR2 = new TLatex(0.20,0.30,"ALICE   pp #sqrt{#it{s}} = 7 TeV   anti-#it{k}_{T} #it{R} = 0.4");
  texR2->SetNDC();
  texR2->SetTextFont(42);
  texR2->SetTextSize(0.09);
  texR2->SetLineWidth(2);
  texR2->Draw();
  
  // legend: upper pad 

  padR[2]->cd();

 padR[2]->cd();
  TLegend *leg2 = new TLegend(0.20,0.05,0.38,0.35,NULL,"brNDC");
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.08);
  leg2->SetTextFont(42);
  leg2->SetLineColor(1);
  leg2->SetLineStyle(1);
  leg2->SetLineWidth(1);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
    /*
  leg2->AddEntry(grRatioP2011[0],"PYTHIA6 Perugia 2011","P");
  leg2->AddEntry(grRatioPNoCR[0],"PYTHIA6 Perugia NoCR","P");
  leg2->AddEntry(grRatioMonash[0],"PYTHIA8 Monash ","P");
*/
  leg2->Draw();

  
  
  if(doWrite){
    cFF->SaveAs("FF_R04_dNdz.pdf");
    cRatio->SaveAs("FF_R04_ratio.pdf");
  }
}
