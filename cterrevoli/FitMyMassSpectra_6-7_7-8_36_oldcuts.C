#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TInterpreter.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TDatabasePDG.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TPaveText.h>

#include "AliAODRecoDecayHF.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliHFMassFitterVAR.h"
#endif


// MACRO to perform fits to D meson invariant mass spectra
// and store raw yields and cut object into a root output file
// Origin: F. Prino (prino@to.infn.it)
// D0: C. Bianchin (cbianchi@pd.infn.it)



//
enum {kD0toKpi, kDplusKpipi};
enum {kBoth, kParticleOnly, kAntiParticleOnly};
enum {kExpo=0, kLinear, kPol2};
enum {kGaus=0, kDoubleGaus, kReflTempl};


// Common variables: to be configured by the user

const Int_t nPtBins=12;
Double_t ptlims[nPtBins+1]={1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,16.,24.,36.};
Int_t rebin[nPtBins+1]=    {4,10,10,10,8,6,12,10,10,14,24,30,38};

Float_t fixSoverRefAt[nPtBins]={-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1}; //for data 0-36
Float_t minRangeFit[nPtBins]={1.72,1.72,1.72,1.72,1.72,1.72,1.72,1.72,1.72,1.72,1.65,1.68};
Float_t maxRangeFit[nPtBins]={2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.05,2.15,2.2};



Int_t typeb=kExpo;
Int_t types=kGaus;
Int_t optPartAntiPart=kBoth;//kAntiParticleOnly;//Both;
Int_t factor4refl=0;
Float_t massRangeForCounting[nPtBins]={0.03,0.03,0.03,0.03,0.035,0.05,0.05,0.055,0.055,0.06,0.08,0.08}; // GeV 0-24

TH2F* hPtMass=0x0;

//for D0only
Bool_t cutsappliedondistr=kFALSE;//kTRUE;
const Int_t nsamples=1;//3;
//Int_t nevents[nsamples]={8.5635859e+07 /*LHC10b+cpart*/,8.9700624e+07/*LHC10dpart*/};
//Int_t nevents[nsamples]={9.0374946e+07 /*LHC10b+c*/,9.7593785e+07/*LHC10d*/};
//Int_t nevents[nsamples]={1.1777812e+08 /*LHC10dnewTPCpid*/,9.7593785e+07/*LHC10d*/};
//Int_t nevents[nsamples]={1.1777812e+08 /*LHC10dnewTPCpid*/,9.7593785e+07/*LHC10d*/,9.0374946e+07 /*LHC10b+c*/};
Int_t nevents[nsamples]={372e+06}; /*LHC10dnewTPCpid*/
//,9.0374946e+07 /*LHC10b+c*/};
//Int_t nevents[nsamples]={1. /*LHC10dnewTPCpid*/,1 /*LHC10dnewTPCpidrescale*/};
// Functions

Bool_t LoadDplusHistos(TObjArray* listFiles, TH1F** hMass);
Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH1F** hMass, const char *CutsType,bool isTHnSparse=false, bool isARtask=true);




void FitMassSpectra(Int_t analysisType=kD0toKpi,
//AnalysisResults_1304b.root  AnalysisResults_1305c.root  AnalysisResults_1306d.root  AnalysisResults_1307e.root 
		    TString fileNameb="../AnalysisResults1406_b.root", 
		    TString fileNamec="../AnalysisResults1407_c.root", 
		    TString fileNamed="../AnalysisResults1408_d.root",
		    TString fileNamee="../AnalysisResults1409_e.root",
		    const char *CutsType="_pass4_central_kFirst3", Int_t useTempl=1,
		    bool isARtask=false){
  //

  gROOT->SetStyle("Plain");	
  //gInterpreter->ExecuteMacro("$ALICE_ROOT/PWGHF/vertexingHF/macros/LoadLibraries.C");
  gStyle->SetOptTitle(0);

  if(useTempl==1) types=kReflTempl;

  TObjArray* listFiles=new TObjArray();
  if(fileNameb!="") listFiles->AddLast(new TObjString(fileNameb.Data()));
  if(fileNamec!="") listFiles->AddLast(new TObjString(fileNamec.Data()));
  if(fileNamed!="") listFiles->AddLast(new TObjString(fileNamed.Data()));
  if(fileNamee!="") listFiles->AddLast(new TObjString(fileNamee.Data()));
  if(listFiles->GetEntries()==0){
    printf("Missing file names in input\n");
    return;
  }

  TH1F** hmass=new TH1F*[nPtBins];
  TH1F** hTemplRefl=new TH1F*[nPtBins];
  TH1F** hSignMC=new TH1F*[nPtBins];
  for(Int_t i=0;i<nPtBins;i++) {
    hmass[i]=0x0;
    hTemplRefl[i]=0x0;
    hSignMC[i]=0x0;
  }

  TString filenameT="";
  TString filenameS="";

  Float_t massD;
  Bool_t retCode;
  if(analysisType==kD0toKpi){

    filenameT=Form("../reflections_fittedDoubleGaus.root");//CHANGE
    filenameS=Form("../ToyInvariantMassDistributions_pass4_central_kFirst3_StepSmearing0.root");//CHANGE
    TFile *ftmp=TFile::Open(filenameT.Data());
    TFile *ftmp2=TFile::Open(filenameS.Data());
    for(Int_t i=0;i<nPtBins;i++){
  cout<< " i  " << i << endl;
      if(useTempl>0){
	hTemplRefl[i]=(TH1F*)ftmp->Get(Form("hRflMCScaledFittedDoubleGaus_ptBin%d_stepSmearing%d",i,0));// THIS IS TO USE ALWAYS THE SAME TEMPLATE
	if(fixSoverRefAt[i]<0){
	  hSignMC[i]=(TH1F*)ftmp2->Get(Form("hSigMCScaled_ptBin%d_stepSmearing%d",i,0));
	  fixSoverRefAt[i]=hTemplRefl[i]->Integral(hTemplRefl[i]->FindBin(minRangeFit[i]*1.0001),hTemplRefl[i]->FindBin(maxRangeFit[i]*0.999))/hSignMC[i]->Integral(hSignMC[i]->FindBin(minRangeFit[i]*1.0001),hSignMC[i]->FindBin(maxRangeFit[i]*0.999));
	  Printf("Sign over Refl bin %d = %f",i,fixSoverRefAt[i]);
	}
      }
    }
    Bool_t isTHnSparse=kTRUE;
    retCode=LoadD0toKpiHistos(listFiles,hmass,CutsType,isTHnSparse,isARtask);
    massD=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  }
  else if(analysisType==kDplusKpipi){
    retCode=LoadDplusHistos(listFiles,hmass);
    massD=TDatabasePDG::Instance()->GetParticle(411)->Mass();
  }
  else{
    printf("Wronganalysistype parameter\n");
    return;
  }
  if(!retCode){
    printf("ERROR in reading input files\n");
    return;
  } 
  
  Printf("massD------->%f",massD);

  TH1D* hCntSig1=new TH1D("hCntSig1","hCntSig1",nPtBins,ptlims);
  TH1D* hCntSig2=new TH1D("hCntSig2","hCntSig2",nPtBins,ptlims);
  TH1D* hNDiffCntSig1=new TH1D("hNDiffCntSig1","hNDiffCntSig1",nPtBins,ptlims);
  TH1D* hNDiffCntSig2=new TH1D("hNDiffCntSig2","hNDiffCntSig2",nPtBins,ptlims);
  TH1D* hSignal=new TH1D("hSignal","hSignal",nPtBins,ptlims);
  TH1D* hRelErrSig=new TH1D("hRelErrSig","hRelErrSig",nPtBins,ptlims);
  TH1D* hInvSignif=new TH1D("hInvSignif","hInvSignif",nPtBins,ptlims);
  TH1D* hBackground=new TH1D("hBackground","hBackground",nPtBins,ptlims);
  TH1D* hBackgroundNormSigma=new TH1D("hBackgroundNormSigma","hBackgroundNormSigma",nPtBins,ptlims);
  TH1D* hSignificance=new TH1D("hSignificance","hSignificance",nPtBins,ptlims);
  TH1D* hMass=new TH1D("hMass","hMass",nPtBins,ptlims);
  TH1D* hSigma=new TH1D("hSigma","hSigma",nPtBins,ptlims);


  Int_t nMassBins=hmass[1]->GetNbinsX();
  Double_t hmin=hmass[1]->GetBinLowEdge(3);
  Double_t hmax=hmass[1]->GetBinLowEdge(nMassBins-2)+hmass[1]->GetBinWidth(nMassBins-2);
  Float_t minBinSum1=hmass[1]->FindBin(massD-massRangeForCounting[1]);
  Float_t maxBinSum1=hmass[1]->FindBin(massD+massRangeForCounting[1]);
  Int_t iPad=1;
  // //RefPreliminary
  // //Set7
   Float_t sigmaRef[nPtBins]={0.008,0.010,0.011,0.013,0.015,0.014,0.015,0.015,0.017,0.018,0.026,0.030};
   //Float_t meanRef[nPtBins]={1.865,1.865,1.865,1.865,1.866,1.866,1.866,1.866,1.867,1.867}; 

  TF1** funBckStore1=new TF1*[nPtBins];
  TF1** funBckStore2=new TF1*[nPtBins];
  TF1** funBckStore3=new TF1*[nPtBins];
  

  AliHFMassFitterVAR** fitter=new AliHFMassFitterVAR*[nPtBins];
  Double_t arrchisquare[nPtBins];
  TCanvas* c1= new TCanvas("c1","MassSpectra",1000,800);
  Int_t nx=4, ny=3;
  if(nPtBins==10){
    nx=4;
    ny=3;
  }
  if(nPtBins==11){
    nx=4;
    ny=3;
  }
  
  if(nPtBins==12){
    nx=4;
    ny=3;
  }
  
  c1->Divide(nx,ny);

  TCanvas *myCanvas[nPtBins];
  TPaveText* ptBin[nPtBins];

  Double_t sig,errsig,s,errs,b,errb;
  for(Int_t iBin=0; iBin<nPtBins; iBin++){

    minBinSum1=hmass[iBin]->FindBin(massD-massRangeForCounting[iBin]);
    maxBinSum1=hmass[iBin]->FindBin(massD+massRangeForCounting[iBin]);

    c1->cd(iPad++);
    gPad->SetTicks();
    
    
    Int_t origNbins=hmass[iBin]->GetNbinsX();
    hmass[iBin]->GetXaxis()->SetTitle("Invariant Mass (K#pi) (GeV/c^{2})");
    hmass[iBin]->GetYaxis()->SetTitle(Form("Entries / %1.0f MeV/c^{2}", (1000*(hmass[iBin]->GetXaxis()->GetBinWidth(1))*rebin[iBin])));
    fitter[iBin]=new AliHFMassFitterVAR(hmass[iBin],hmin, hmax,rebin[iBin],typeb,types);
   
    if(useTempl>=1){
      fitter[iBin]->SetTemplateReflections(hTemplRefl[iBin]);
      fitter[iBin]->SetFixReflOverS(fixSoverRefAt[iBin],kTRUE);//fixSoverRefAt[iBin]);
    }

    fitter[iBin]->SetRangeFit(minRangeFit[iBin],maxRangeFit[iBin]);
    
   
    rebin[iBin]=origNbins/fitter[iBin]->GetBinN();
    Printf("*******************");
    Printf("REBIN = %d", rebin[iBin]);
    Printf("*******************");
    fitter[iBin]->SetReflectionSigmaFactor(factor4refl);
    fitter[iBin]->SetInitialGaussianMean(massD);
//  fitter[iBin]->SetFixGaussianMean(meanRef[iBin]);
    if(iBin==11)fitter[iBin]->SetFixGaussianSigma(sigmaRef[iBin]);
    Bool_t out=fitter[iBin]->MassFitter(0);
    if(!out) {
      fitter[iBin]->GetHistoClone()->Draw();
      continue;
    }
    Double_t mass=fitter[iBin]->GetMean();
    Double_t massErr=fitter[iBin]->GetMeanUncertainty();
    Double_t sigma=fitter[iBin]->GetSigma();
    Double_t sigmaErr=fitter[iBin]->GetSigmaUncertainty();
 
    Float_t halfRange=3*sigma;
    Float_t minBinSum=hmass[iBin]->FindBin(massD-halfRange);
    Float_t maxBinSum=hmass[iBin]->FindBin(massD+halfRange);
    Printf("minBinSum (5sigma)=%f\t minBinSum (100MeV)=%f", minBinSum, minBinSum1);
    Printf("maxBinSum (5sigma)=%f\t maxBinSum (100MeV)=%f", maxBinSum, maxBinSum1);

    arrchisquare[iBin]=fitter[iBin]->GetReducedChiSquare();
    TF1* fB1=fitter[iBin]->GetBackgroundFullRangeFunc();
    TF1* fB2=fitter[iBin]->GetBackgroundRecalcFunc();
    TF1* fM=fitter[iBin]->GetMassFunc();
    /*if(iBin==0 && fB1) funBckStore1=new TF1(*fB1);
    if(iBin==0 && fB2) funBckStore2=new TF1(*fB2);
    if(iBin==0 && fM) funBckStore3=new TF1(*fM);*/
    funBckStore1[iBin]=(TF1*)fB1->Clone(Form("BkgFunction1_ptBin%d",iBin));
    funBckStore2[iBin]=(TF1*)fB2->Clone(Form("BkgFunction2_ptBin%d",iBin));
    funBckStore3[iBin]=(TF1*)fM->Clone(Form("MassFunction_ptBin%d",iBin));

    fitter[iBin]->DrawHere(gPad);
    fitter[iBin]->Signal(3,s,errs);
    fitter[iBin]->Background(3,b,errb);
    fitter[iBin]->Significance(3,sig,errsig);
    Double_t ry=fitter[iBin]->GetRawYield();
    Double_t ery=fitter[iBin]->GetRawYieldError();

    ptBin[iBin] = new TPaveText(0.6,0.55,0.8,0.6,"NDC");
    ptBin[iBin]->SetFillStyle(0);
    ptBin[iBin]->SetBorderSize(0);
    ptBin[iBin]->AddText(0.,0.,Form("%1.1f<p_{T}<%1.1f",ptlims[iBin],ptlims[iBin+1]));
    ptBin[iBin]->SetTextFont(42);
    ptBin[iBin]->SetTextSize(0.08);
    gPad->cd();
    ptBin[iBin]->Draw();
    gPad->Update();

    myCanvas[iBin] = new TCanvas(Form("myCanvas_%d",iBin),Form("Invariant mass pt bin %d",iBin));
    fitter[iBin]->DrawHere(gPad);
    gPad->SetTicks();
    gPad->cd();
    ptBin[iBin]->Draw();
    gPad->Update();

    Float_t cntSig1=0.;
    Float_t cntSig2=0.;
    Float_t cntErr=0.;
    for(Int_t iMB=minBinSum1; iMB<=maxBinSum1; iMB++){
      Float_t bkg1=fB1 ? fB1->Eval(hmass[iBin]->GetBinCenter(iMB))/rebin[iBin] : 0;
      Float_t bkg2=fB2 ? fB2->Eval(hmass[iBin]->GetBinCenter(iMB))/rebin[iBin] : 0;
      cntSig1+=(hmass[iBin]->GetBinContent(iMB)-bkg1);
      cntSig2+=(hmass[iBin]->GetBinContent(iMB)-bkg2);
      cntErr+=(hmass[iBin]->GetBinContent(iMB));
    }
    hCntSig1->SetBinContent(iBin+1,cntSig1);
    hCntSig1->SetBinError(iBin+1,TMath::Sqrt(cntErr));
    // hNDiffCntSig1->SetBinContent(iBin+1,(s-cntSig1)/s);
    // hNDiffCntSig1->SetBinError(iBin+1,TMath::Sqrt(cntErr)/s);
    hNDiffCntSig1->SetBinContent(iBin+1,(ry-cntSig1)/ry);
    hNDiffCntSig1->SetBinError(iBin+1,TMath::Sqrt(cntErr)/ry);
    hCntSig2->SetBinContent(iBin+1,cntSig2);
    // hNDiffCntSig2->SetBinContent(iBin+1,(s-cntSig2)/s);
    // hNDiffCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr)/s);
    hNDiffCntSig2->SetBinContent(iBin+1,(ry-cntSig2)/ry);
    hNDiffCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr)/ry);
    hCntSig2->SetBinError(iBin+1,TMath::Sqrt(cntErr));
    // hSignal->SetBinContent(iBin+1,s);
    // hSignal->SetBinError(iBin+1,errs);
    hSignal->SetBinContent(iBin+1,ry);
    hSignal->SetBinError(iBin+1,ery);
    hRelErrSig->SetBinContent(iBin+1,errs/s);
    hInvSignif->SetBinContent(iBin+1,1/sig);
    hInvSignif->SetBinError(iBin+1,errsig/(sig*sig));
    hBackground->SetBinContent(iBin+1,b); //consider sigma
    hBackground->SetBinError(iBin+1,errb);
    hBackgroundNormSigma->SetBinContent(iBin+1,b/(3*fitter[iBin]->GetSigma())*(3*0.012)); //consider sigma
    hBackgroundNormSigma->SetBinError(iBin+1,errb);
    hSignificance->SetBinContent(iBin+1,sig);
    hSignificance->SetBinError(iBin+1,errsig);
    hMass->SetBinContent(iBin+1,mass);
    hMass->SetBinError(iBin+1,massErr);
    hSigma->SetBinContent(iBin+1,sigma);
    hSigma->SetBinError(iBin+1,sigmaErr);
    
  }
  
  Float_t min=0.,max=0.,interval=0.;
  for(Int_t iBin=0; iBin<nPtBins; iBin++){
    //sig=hSigma->GetBinContent(iBin+1);
    interval=maxBinSum1-minBinSum1; //5*sig;
    min=massD-interval;
    max=massD+interval;
    Printf("pTbin=%d\t min=%f\t max=%f\n",iBin+1,min,max);
  }

  /*
  c1->cd(1); // is some cases the fitting function of 1st bin get lost
  funBckStore1->Draw("same");
  funBckStore2->Draw("same");
  funBckStore3->Draw("same");
  */

  TCanvas *cpar=new TCanvas("cpar","Fit params",1200,600);
  cpar->Divide(2,1);
  cpar->cd(1);
  hMass->SetMarkerStyle(20);
  hMass->Draw("PE");
  hMass->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hMass->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
  cpar->cd(2);
  hSigma->SetMarkerStyle(20);
  hSigma->Draw("PE");
  hSigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hSigma->GetYaxis()->SetTitle("Sigma (GeV/c^{2})");

  TCanvas* csig=new TCanvas("csig","Results",1200,600);
  csig->Divide(3,1);
  csig->cd(1);
  hSignal->SetMarkerStyle(20);
  hSignal->SetMarkerColor(4);
  hSignal->SetLineColor(4);
  hSignal->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hSignal->GetYaxis()->SetTitle("Signal");
  hSignal->Draw("P");
  hCntSig1->SetMarkerStyle(26);
  hCntSig1->SetMarkerColor(2);
  hCntSig1->SetLineColor(2);
  hCntSig1->Draw("PSAME");
  hCntSig2->SetMarkerStyle(29);
  hCntSig2->SetMarkerColor(kGray+1);
  hCntSig2->SetLineColor(kGray+1);
  hCntSig2->Draw("PSAME");
  TLegend* leg=new TLegend(0.4,0.7,0.89,0.89);
  leg->SetFillColor(0);
  TLegendEntry* ent=leg->AddEntry(hSignal,"From Fit","PL");
  ent->SetTextColor(hSignal->GetMarkerColor());
  ent=leg->AddEntry(hCntSig1,"Bin Count. only bkg","PL");
  ent->SetTextColor(hCntSig1->GetMarkerColor());
  ent=leg->AddEntry(hCntSig2,"Bon Count. sign+bkg","PL");
  ent->SetTextColor(hCntSig2->GetMarkerColor());
  leg->Draw();
  csig->cd(2);
  hBackground->SetMarkerStyle(20);
  hBackground->Draw("P");
  hBackground->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hBackground->GetYaxis()->SetTitle("Background");
  csig->cd(3);
  hSignificance->SetMarkerStyle(20);
  hSignificance->Draw("P");
  hSignificance->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hSignificance->GetYaxis()->SetTitle("Significance");

  TCanvas* cDiffS=new TCanvas("cDiffS","Difference",1200,600);
  cDiffS->Divide(2,1);
  cDiffS->cd(1);
  hRelErrSig->SetMarkerStyle(20); //fullcircle
  hRelErrSig->SetTitleOffset(1.2);  
  hRelErrSig->SetTitle("Rel Error from Fit;p_{T} (GeV/c);Signal Relative Error");
  hRelErrSig->Draw("P");
  hInvSignif->SetMarkerStyle(21); //fullsquare
  hInvSignif->SetMarkerColor(kMagenta+1);
  hInvSignif->SetLineColor(kMagenta+1);
  hInvSignif->Draw("PSAMES");
  TLegend* leg2=new TLegend(0.4,0.7,0.89,0.89);
  leg2->SetFillColor(0);
  TLegendEntry* ent2=leg2->AddEntry(hRelErrSig,"From Fit","P");
  ent2->SetTextColor(hRelErrSig->GetMarkerColor());
  ent2=leg2->AddEntry(hInvSignif,"1/Significance","PL");
  ent2->SetTextColor(hInvSignif->GetMarkerColor());
  leg2->Draw();

  cDiffS->cd(2);
  hNDiffCntSig1->SetMarkerStyle(26);
  hNDiffCntSig1->SetMarkerColor(2);
  hNDiffCntSig1->SetLineColor(2);
  hNDiffCntSig1->SetTitle("Cmp Fit-Count;p_{T} (GeV/c);(S_{fit}-S_{count})/S_{fit}");
  hNDiffCntSig1->Draw("P");
  hNDiffCntSig2->SetMarkerStyle(29);
  hNDiffCntSig2->SetMarkerColor(kGray+1);
  hNDiffCntSig2->SetLineColor(kGray+1);
  hNDiffCntSig2->Draw("PSAME");
  TLegend* leg1=new TLegend(0.4,0.7,0.89,0.89);
  leg1->SetFillColor(0);
  TLegendEntry* ent1=leg1->AddEntry(hNDiffCntSig1,"Bin Count. only bkg","PL");
  ent1->SetTextColor(hNDiffCntSig1->GetMarkerColor());
  ent1=leg1->AddEntry(hNDiffCntSig2,"Bin Count. sign+bkg","PL");
  ent1->SetTextColor(hNDiffCntSig2->GetMarkerColor());
  leg1->Draw();

  TCanvas *cDiffonly = new TCanvas("SystErrBinCount","SystErrBinCount");
  cDiffonly->cd();
  hNDiffCntSig1->SetMarkerStyle(22);
  hNDiffCntSig1->SetMarkerColor(2);
  hNDiffCntSig1->SetMarkerSize(1.3);
  hNDiffCntSig1->SetLineColor(2);
  hNDiffCntSig1->SetLineWidth(2);
  hNDiffCntSig1->SetTitle("Cmp Fit-Count;p_{T} (GeV/c);(S_{fit}-S_{count})/S_{fit}");
  hNDiffCntSig1->Draw("P");
  hNDiffCntSig2->SetMarkerStyle(20);
  hNDiffCntSig2->SetMarkerColor(kBlue);
  hNDiffCntSig2->SetLineColor(kBlue);
  hNDiffCntSig2->SetLineWidth(2);
  hNDiffCntSig2->SetMarkerSize(1.3);
  hNDiffCntSig2->Draw("PSAME");
  TLegend* legOnly=new TLegend(0.4,0.7,0.89,0.89);
  legOnly->SetFillColor(0);
  TLegendEntry* entOnly=legOnly->AddEntry(hNDiffCntSig1,"Bin Count. only bkg","PL");
  entOnly->SetTextColor(hNDiffCntSig1->GetMarkerColor());
  entOnly=legOnly->AddEntry(hNDiffCntSig2,"Bin Count. sign+bkg","PL");
  entOnly->SetTextColor(hNDiffCntSig2->GetMarkerColor());
  TPaveText* tsyst=new TPaveText(0.2,0.90,0.8,0.95,"NDC");
  tsyst->SetFillStyle(0);
  tsyst->SetBorderSize(0);
  tsyst->AddText(0.,0.,"D^{0} #rightarrow K#pi p-Pb, #sqrt{s_{NN}} = 2.76 TeV, Min. bias");
  tsyst->SetTextFont(42);
  tsyst->Draw();

  TPaveText* tsystlabel=new TPaveText(0.8,0.90,1.,0.95,"NDC");
  tsystlabel->SetFillStyle(0);
  tsystlabel->SetBorderSize(0);
  tsystlabel->AddText(0.,0.,"ALICE Performance");
  tsystlabel->SetTextFont(42);
  tsystlabel->Draw();
 
  legOnly->Draw();

  

  TGraph* grReducedChiSquare=new TGraph(nPtBins,ptlims,arrchisquare);
  grReducedChiSquare->SetName("grReducedChiSquare");
  grReducedChiSquare->SetTitle("Reduced Chi2;p_{T} (GeV/c);#tilde{#chi}^{2}");
  TCanvas *cChi2=new TCanvas("cChi2","reduced chi square",600,600);
  cChi2->cd();
  grReducedChiSquare->SetMarkerStyle(21);
  grReducedChiSquare->Draw("AP");

  TCanvas* cbkgNormSigma=new TCanvas("cbkgNormSigma","Background normalized to sigma",400,600);
  cbkgNormSigma->cd();
  hBackgroundNormSigma->SetMarkerStyle(20);
  hBackgroundNormSigma->Draw("P");
  hBackgroundNormSigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hBackgroundNormSigma->GetYaxis()->SetTitle("Background #times 3 #times 0.012/ (3 #times #sigma)");
  hBackgroundNormSigma->Draw();
 


  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) {
    if(analysisType==kD0toKpi) partname="D0";
    if(analysisType==kDplusKpipi) partname="Dplus";
  }
  if(optPartAntiPart==kAntiParticleOnly) {
    if(analysisType==kD0toKpi) partname="D0bar";
    if(analysisType==kDplusKpipi) partname="Dminus";
  }

  TString templ="_NoTemplate";
  if(useTempl>=1) templ="_Template";

  //  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"update");
  TFile* outf=new TFile(Form("RawYield%s_%s%s_topofinal_ndl1Last.root",partname.Data(),CutsType,templ.Data()),"recreate");
  outf->cd();
  hMass->Write();
  hSigma->Write();
  hCntSig1->Write();
  hCntSig2->Write();
  hNDiffCntSig1->Write();
  hNDiffCntSig2->Write();
  hSignal->Write();
  hRelErrSig->Write();
  hInvSignif->Write();
  hBackground->Write();
  hBackgroundNormSigma->Write();
  hSignificance->Write();
  grReducedChiSquare->Write();
  //  hPtMass->Write();

//cri
// TFile *miofile =new TFile("massfile_central.root","recreate");
//  for(Int_t i=0;i<nPtBins;i++){
//  hmass[i]->Write();
//  }
// miofile->Close();
//cri


  for(Int_t iBin=0; iBin<nPtBins; iBin++){
    outf->cd();
    funBckStore1[iBin]->Write();
    funBckStore2[iBin]->Write();
    funBckStore3[iBin]->Write();
  }

  outf->Close();
}


Bool_t LoadDplusHistos(TObjArray* listFiles, TH1F** hMass){

  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  AliRDHFCutsDplustoKpipi** dcuts=new AliRDHFCutsDplustoKpipi*[nFiles];

  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){
    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    
    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName()); 
    TDirectory *dir = (TDirectory*)f->Get("PWG3_D2H_InvMassDplus");
    if(!dir){
      printf("ERROR: directory PWG3_D2H_InvMassDplus not found in %s\n",fName.Data());
      continue;
    }
    hlist[nReadFiles]=(TList*)dir->Get("coutputDplus");
    TList *listcut = (TList*)dir->Get("coutputDplusCuts");
    dcuts[nReadFiles]=(AliRDHFCutsDplustoKpipi*)listcut->At(1);
    if(nReadFiles>0){
      Bool_t sameCuts=dcuts[nReadFiles]->CompareCuts(dcuts[0]);
      if(!sameCuts){
	printf("ERROR: Cut objects do not match\n");
	return kFALSE;
      }
    }
    nReadFiles++;
  }
  if(nReadFiles<nFiles){
    printf("WARNING: not all requested files have been found\n");
  }

  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1];//.;

  Int_t iFinBin=0;
  for(Int_t i=0;i<nPtBinsCuts;i++){
    if(ptlimsCuts[i]>=ptlims[iFinBin+1]) iFinBin+=1; 
    if(iFinBin>nPtBins) break;
    if(ptlimsCuts[i]>=ptlims[iFinBin] && 
       ptlimsCuts[i+1]<=ptlims[iFinBin+1]){
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	TString histoName;
	if(optPartAntiPart==kBoth) histoName.Form("hMassPt%dTC",i);
	else if(optPartAntiPart==kParticleOnly) histoName.Form("hMassPt%dTCPlus",i);
	else if(optPartAntiPart==kAntiParticleOnly) histoName.Form("hMassPt%dTCMinus",i);
	TH1F* htemp=(TH1F*)hlist[iFile]->FindObject(histoName.Data());
	if(!htemp){
	  printf("ERROR: Histogram %s not found\n",histoName.Data());
	  return kFALSE;
	}
	if(!hMass[iFinBin]){
	  hMass[iFinBin]=new TH1F(*htemp);
	}else{
	  hMass[iFinBin]->Add(htemp);
	}
      }
    }
  }
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="Dplus";
  if(optPartAntiPart==kAntiParticleOnly) partname="Dminus";

  for(Int_t iFile=0; iFile<nReadFiles; iFile++){
    TH2F* htemp2=(TH2F*)hlist[iFile]->FindObject("hPtVsMassTC");
    if(iFile==0){
      hPtMass=new TH2F(*htemp2);
    }else{
      hPtMass->Add(htemp2);
    }
  }

  TFile* outf=new TFile(Form("RawYield%s.root",partname.Data()),"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();

  return kTRUE;

}

Bool_t LoadD0toKpiHistos(TObjArray* listFiles, TH1F** hMass, const char *CutsType, bool isTHnSparse, bool isARtask){
  //
  Int_t nFiles=listFiles->GetEntries();
  TList **hlist=new TList*[nFiles];
  AliRDHFCutsD0toKpi** dcuts=new AliRDHFCutsD0toKpi*[nFiles];

  Int_t nReadFiles=0;
  for(Int_t iFile=0; iFile<nFiles; iFile++){


    TString fName=((TObjString*)listFiles->At(iFile))->GetString();    

    TFile *f=TFile::Open(fName.Data());
    if(!f){
      printf("ERROR: file %s does not exist\n",fName.Data());
      continue;
    }
    printf("Open File %s\n",f->GetName());

    TString dirname="PWG3_D2H_D0InvMass";
    if(optPartAntiPart==kParticleOnly) dirname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) dirname+="D0bar";
    if(cutsappliedondistr) dirname+="C";
    if(isARtask) dirname="PWG3_D2H_d0D0";
    else {
      dirname += CutsType;
    }
    TDirectory *dir = (TDirectory*)f->Get(dirname);
    if(!dir){
      printf("ERROR: directory %s not found in %s\n",dirname.Data(),fName.Data());
      continue;
    }
    TString listmassname="coutputmassD0Mass";
    if(optPartAntiPart==kParticleOnly) listmassname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) listmassname+="D0bar";
    if(cutsappliedondistr) listmassname+="C";
    if(isARtask) {
      if( strcmp(CutsType,"tight")==0){
	listmassname="clistTGHCsign_d0D0";
      }
      else if( strcmp(CutsType,"loose")==0){
	listmassname="clistLSCsign_d0D0";
      }
      else if( strcmp(CutsType,"noPID")==0){
	listmassname="clistNCsign_d0D0";
      }
    }
    else {
      listmassname += CutsType;
      listmassname += "0100";
    }
    printf("List name %s\n",listmassname.Data());
    hlist[nReadFiles]=(TList*)dir->Get(listmassname);

    TString cutsobjname="cutsD0";
    if(optPartAntiPart==kParticleOnly) cutsobjname+="D0";
    else if(optPartAntiPart==kAntiParticleOnly) cutsobjname+="D0bar";
    if(cutsappliedondistr) cutsobjname+="C";
    if(isARtask) {
      if( strcmp(CutsType,"tight")==0){
	cutsobjname="ccutsObjectTight_d0D0";
      }
      else if( strcmp(CutsType,"loose")==0 || strcmp(CutsType,"noPID")==0 ){
	cutsobjname="ccutsObjectLoose_d0D0";
      }
    }else{
      cutsobjname += CutsType;
      cutsobjname += "0100";
    }
    printf("Cuts name %s\n",cutsobjname.Data());

    dcuts[nReadFiles]=(AliRDHFCutsD0toKpi*)dir->Get(cutsobjname);
    if(!dcuts[nReadFiles]) {
      printf("ERROR: Cut objects do not match\n");
      return kFALSE;
    }
    if(nReadFiles>0){
      Bool_t sameCuts=dcuts[nReadFiles]->CompareCuts(dcuts[0]);
      if(!sameCuts){
	printf("ERROR: Cut objects do not match\n");
	return kFALSE;
      }
    }
    nReadFiles++;
  }
  if(nReadFiles<nFiles){
    printf("WARNING: not all requested files have been found\n");
    if (nReadFiles==0) {
      printf("ERROR: Any file/dir found\n");
      return kFALSE;
    }
  }
  //  printf("Cuts type %s, particle/antipart %d\n",CutsType,optPartAntiPart);

  Int_t nPtBinsCuts=dcuts[0]->GetNPtBins();
  printf("Number of pt bins for cut object = %d\n",nPtBins);
  Float_t *ptlimsCuts=dcuts[0]->GetPtBinLimits();
  //ptlimsCuts[nPtBinsCuts]=ptlimsCuts[nPtBinsCuts-1]+4.;

  printf("Loop over the pt bins\n");
  Float_t ptlimsCutsT[13]={1.,2.,3.,4.,5.,6.,7.,8.,10.,12.,16.,24.,36.};
  for(Int_t i=0;i<12;i++){
    cout<< " loop su nPtBinsCuts " << i << endl;
      for(Int_t iFile=0; iFile<nReadFiles; iFile++){
	//	TH1F* htemp=(TH1F*)hlist[iFile]->FindObject(Form("histMass_%d",i));
	const char *histoname;
        if(!isARtask) histoname = Form("histMass_%d",i);
	TH1F * htemp;
	if(!isARtask) { 
          if(isTHnSparse){
          //qui qui qui
            cout<< " sono nel caso di thnsparse " << endl;
            THnSparseF *hspSign=(THnSparseF*)hlist[iFile]->FindObject("fhStudyImpParSingleTrackCand");
            cout<< " hspSign " << hspSign << " file " << iFile<< endl;
            TAxis *axpt=hspSign->GetAxis(0);//pt
            axpt->SetRangeUser(1.001*ptlimsCutsT[i],0.999*ptlimsCutsT[i+1]);
            TAxis *axInvMassD0=hspSign->GetAxis(5);//invMAssD0
            TAxis *axNormResImpPartTrk1=hspSign->GetAxis(1);//normIP1
            TAxis *axNormResImpPartTrk2=hspSign->GetAxis(2);//normIP2
            TAxis *axLxy=hspSign->GetAxis(3);//lxy
            TAxis *axNormLxy=hspSign->GetAxis(4);//lyx norm
            TAxis *axCut=hspSign->GetAxis(6);//cut 
            TAxis *axPID=hspSign->GetAxis(7);//pid 
            TAxis *axD0=hspSign->GetAxis(8);//d0d0bar             
//            if(ptlimsCutsT[i]>=12. && ptlimsCutsT[i]<24.)  axNormLxy->SetRangeUser(3.,15.);
            if(ptlimsCutsT[i]>=24.)  axNormLxy->SetRangeUser(1.,15.);
//            axNormLxy->SetRangeUser(1.001*4.,0.999*15.);//ndlcut


            if(ptlimsCutsT[i]>=1 && ptlimsCutsT[i]<=5){
              axNormResImpPartTrk1->SetRangeUser(-2.,2.);//2
              axNormResImpPartTrk2->SetRangeUser(-2.,2.);//2
              }
            if(ptlimsCutsT[i]>=5 && ptlimsCutsT[i]<12){
              axNormResImpPartTrk1->SetRangeUser(-3.,3.);//3
              axNormResImpPartTrk2->SetRangeUser(-3.,3.);//3
              }
            if(ptlimsCutsT[i]>=12 && ptlimsCutsT[i]<=16){
              axNormResImpPartTrk1->SetRangeUser(-4.,4.);//4
              axNormResImpPartTrk2->SetRangeUser(-4.,4.);//4
              }
            if(ptlimsCutsT[i]>=16 && ptlimsCutsT[i]<=24){
              axNormResImpPartTrk1->SetRangeUser(-3.,3.);//3
              axNormResImpPartTrk2->SetRangeUser(-3.,3.);//3
              }
            if(ptlimsCutsT[i]>=24){
              axNormResImpPartTrk1->SetRangeUser(-2.5*0.999,2.5*0.999);//2
              axNormResImpPartTrk2->SetRangeUser(-2.5*0.999,2.5*0.999);//2
              }

            axCut->SetRange(1,1);//D0
            axPID->SetRange(2,2);//D0
            axD0->SetRange(1,1);//D0
            htemp=(TH1F*)hspSign->Projection(5);

            axCut->SetRange(2,2);//D0bar
            axPID->SetRange(3,3);//D0bar
            axD0->SetRange(2,2);//D0bar
            htemp->Add(hspSign->Projection(5));

             axCut->SetRange(1,1); //D0
            axPID->SetRange(4,4);//both
            axD0->SetRange(1,1);//D0
            htemp->Add(hspSign->Projection(5));

            axCut->SetRange(2,2); //D0bar
            axPID->SetRange(4,4);//both
            axD0->SetRange(2,2);//D0bar
            htemp->Add(hspSign->Projection(5));

            axCut->SetRange(3,3); //D0both
            axPID->SetRange(4,4);//both
            axD0->SetRangeUser(0,2);//D0both
            htemp->Add(hspSign->Projection(5));

            axCut->SetRange(3,3); //D0both
            axPID->SetRange(2,2);//D0
            axD0->SetRange(1,1);//D0
            htemp->Add(hspSign->Projection(5));

            axCut->SetRange(3,3); //D0 both
            axPID->SetRange(3,3);//D0 bar
            axD0->SetRange(2,2);//D0bar
            htemp->Add(hspSign->Projection(5));

	}else {
//        htemp=(TH1F*)hlist[iFile]->FindObject(Form("%s%d",histoname,i));
          htemp=(TH1F*)hlist[iFile]->FindObject(Form("histMass_%d",i));//CRIIII;
        }
        }
        else {
          printf(" Looking for histo %s in pt bin %d",histoname,i);
          htemp=(TH1F*)hlist[iFile]->FindObject(Form("%s%d",histoname,i));
        }

        if(!hMass[i]){
          hMass[i]=new TH1F(*htemp);
        }else{
          cout<< " lo aggiungo " << endl;
          hMass[i]->Add(htemp);
        }
 
                if(optPartAntiPart==kBoth && isARtask){
          if( strcmp(CutsType,"tight")==0 ) histoname="hInvMassD0barTGHCsign_pt";
          if( strcmp(CutsType,"loose")==0 ) histoname="hInvMassD0barLSCsign_pt";
          printf(" Looking for histo %s in pt bin %d",histoname,i);
          TH1F* htemp2=(TH1F*)hlist[iFile]->FindObject(Form("%s%d",histoname,i));
          if(!hMass[i]){
            hMass[i]=new TH1F(*htemp2);
          }else{
            hMass[i]->Add(htemp2);
          }
        }
      }
}
  TString partname="Both";
  if(optPartAntiPart==kParticleOnly) partname="D0";
  if(optPartAntiPart==kAntiParticleOnly) partname="D0bar";

//   for(Int_t iFile=0; iFile<nReadFiles; iFile++){
//     TH2F* htemp2=(TH2F*)hlist[iFile]->FindObject("hPtVsMassTC");
//     if(iFile==0){
//       hPtMass=new TH2F(*htemp2);
//     }else{
//       hPtMass->Add(htemp2);
//     }
//   }
  
//cri
 TFile *miofile =new TFile("massfile_central.root","recreate");
  for(Int_t i=0;i<nPtBins;i++){
  hMass[i]->Write();
  }
 miofile->Close();
//cri

  TFile* outf=new TFile(Form("RawYield%s_%s_topofinal_ndl1Last.root",partname.Data(),CutsType),"recreate");
  outf->cd();
  dcuts[0]->Write();
  outf->Close();
  return kTRUE;
}

void CompareFitTypes(TString* paths, TString* legtext,Int_t ncmp=3,TString* filenameYield=0x0){
  //read ncmp RawYield.roots and draw them together
  //set the global variable nevents before running
  //arguments:
  // - paths= vector of ncmp dimension with the paths of the file RawYield.root
  // - legtext= vector of ncmp dimension with the label for the legend
  // - ncmp= number of files to compare (default is 3)
  // - filenameYield= optional ncmp-dimensional array with the difference between the names of the files to be compared (useful if the 2 files are in the same directory but have different names)

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFrameBorderMode(0);

  if(!filenameYield) filenameYield=new TString[ncmp];

  for(Int_t k=0;k<ncmp;k++){
    if(!filenameYield) filenameYield[k]="RawYield.root";
    filenameYield[k].Prepend(paths[k]);
  }
  
  TCanvas* cSig=new TCanvas("cSig","Results",1200,600);
  cSig->Divide(3,1);
  TCanvas* cBkgN=new TCanvas("cBkgN","Background normalized to sigma",400,600);
  TCanvas* cDiffS=new TCanvas("cDiffS","Difference",1200,600);
  cDiffS->Divide(2,1);
  TCanvas *cChi2=new TCanvas("cChi2","reduced chi square",600,600);

  TLegend* leg1=new TLegend(0.4,0.7,0.89,0.89);
  leg1->SetFillColor(0);
  TLegend* leg2=(TLegend*)leg1->Clone();
  TLegend* leg3=(TLegend*)leg1->Clone();
  TLegend* leg4=new TLegend(0.4,0.6,0.8,0.89);
  leg4->SetFillColor(0);

  for(Int_t iTypes=0;iTypes<ncmp;iTypes++){
    TFile* fin=new TFile(filenameYield[iTypes]);
    if(!fin){
      printf("WARNING: %s not found",filenameYield[iTypes].Data());
      continue;
    }

    TH1F* hSignal=(TH1F*)fin->Get("hSignal");
    TH1F* hBackground=(TH1F*)fin->Get("hBackground");
    TH1F* hBackgroundNormSigma=(TH1F*)fin->Get("hBackgroundNormSigma");
    TH1F* hSignificance=(TH1F*)fin->Get("hSignificance");
    hSignal->SetName(Form("%s%d",hSignal->GetName(),iTypes));
    hBackground->SetName(Form("%s%d",hBackground->GetName(),iTypes));
    hBackgroundNormSigma->SetName(Form("%s%d",hBackgroundNormSigma->GetName(),iTypes));
    hSignificance->SetName(Form("%s%d",hSignificance->GetName(),iTypes));

    hSignal->SetMarkerColor(iTypes+2);
    hSignal->SetLineColor(iTypes+2);
    hBackground->SetMarkerColor(iTypes+2);
    hBackground->SetLineColor(iTypes+2);
    hBackgroundNormSigma->SetMarkerColor(iTypes+2);
    hBackgroundNormSigma->SetLineColor(iTypes+2);
    hSignificance->SetMarkerColor(iTypes+2);
    hSignificance->SetLineColor(iTypes+2);

    TLegendEntry* ent4=leg4->AddEntry(hSignal,Form("%s",legtext[iTypes].Data()),"PL");
    ent4->SetTextColor(hSignal->GetMarkerColor());

    if(ncmp==nsamples){
      printf("Info: Normalizing signal, background and significance to the number of events\n");
      hSignal->Scale(1./nevents[iTypes]);
      hBackground->Scale(1./nevents[iTypes]);
      hBackgroundNormSigma->Scale(1./nevents[iTypes]);
      hSignificance->Scale(1./TMath::Sqrt(nevents[iTypes]));
    }

    if(iTypes==0){
      cSig->cd(1);
      hSignal->DrawClone("P");
      cSig->cd(2);
      hBackground->DrawClone("P");
      cSig->cd(3);
      hSignificance->DrawClone("P");
      cBkgN->cd();
      hBackgroundNormSigma->DrawClone("P");
    } else{
      cSig->cd(1);
      hSignal->DrawClone("Psames");
      cSig->cd(2);
      hBackground->DrawClone("Psames");
      cSig->cd(3);
      hSignificance->DrawClone("Psames");
      cBkgN->cd();
      hBackgroundNormSigma->DrawClone("Psames");
    }

    TH1F* hRelErrSig=(TH1F*)fin->Get("hRelErrSig");
    TH1F* hInvSignif=(TH1F*)fin->Get("hInvSignif");
    hRelErrSig->SetName(Form("%s%d",hRelErrSig->GetName(),iTypes));
    hInvSignif->SetName(Form("%s%d",hInvSignif->GetName(),iTypes));

    hRelErrSig->SetMarkerColor(iTypes+2);
    hRelErrSig->SetLineColor(iTypes+2);
    hInvSignif->SetMarkerColor(iTypes+2);
    hInvSignif->SetLineColor(iTypes+2);

    TLegendEntry* ent1=leg1->AddEntry(hRelErrSig,Form("From Fit (%s)",legtext[iTypes].Data()),"P");
    ent1->SetTextColor(hRelErrSig->GetMarkerColor());
    ent1=leg1->AddEntry(hInvSignif,Form("1/Significance (%s)",legtext[iTypes].Data()),"PL");
    ent1->SetTextColor(hInvSignif->GetMarkerColor());

    cDiffS->cd(1);
    if(iTypes==0){
      hRelErrSig->DrawClone("P");
      hInvSignif->DrawClone();
    } else{
      hRelErrSig->DrawClone("Psames");
      hInvSignif->DrawClone("sames");
    }

    TH1F* hNDiffCntSig1=(TH1F*)fin->Get("hNDiffCntSig1");
    TH1F* hNDiffCntSig2=(TH1F*)fin->Get("hNDiffCntSig2");
    hNDiffCntSig1->SetName(Form("%s%d",hNDiffCntSig1->GetName(),iTypes));
    hNDiffCntSig2->SetName(Form("%s%d",hNDiffCntSig2->GetName(),iTypes));

    hNDiffCntSig1->SetMarkerColor(iTypes+2);
    hNDiffCntSig1->SetLineColor(iTypes+2);
    hNDiffCntSig2->SetMarkerColor(iTypes+2);
    hNDiffCntSig2->SetLineColor(iTypes+2);
    TLegendEntry* ent2=leg2->AddEntry(hNDiffCntSig1,Form("Bin Count. only bkg (%s)",legtext[iTypes].Data()),"PL");
    ent2->SetTextColor(hNDiffCntSig1->GetMarkerColor());
    ent2=leg2->AddEntry(hNDiffCntSig2,Form("Bin Count. signal+bkg (%s)",legtext[iTypes].Data()),"PL");
    ent2->SetTextColor(hNDiffCntSig2->GetMarkerColor());

    cDiffS->cd(2);
    if(iTypes==0){
      hNDiffCntSig1->DrawClone();
      hNDiffCntSig2->DrawClone();
    }else{
     hNDiffCntSig1->DrawClone("sames");
     hNDiffCntSig2->DrawClone("sames");
    }

    TGraph* grReducedChiSquare=(TGraph*)fin->Get("grReducedChiSquare");
    grReducedChiSquare->SetName(Form("%s%d",grReducedChiSquare->GetName(),iTypes));

    grReducedChiSquare->SetMarkerColor(iTypes+2);
    grReducedChiSquare->SetLineColor(iTypes+2);
    TLegendEntry* ent3=leg3->AddEntry(grReducedChiSquare,Form("%s",legtext[iTypes].Data()),"PL");
    ent3->SetTextColor(grReducedChiSquare->GetMarkerColor());

    cChi2->cd();
    if(iTypes==0) grReducedChiSquare->DrawClone("AP");
    else grReducedChiSquare->DrawClone("P");
  }

  cSig->cd(1);
  leg4->Draw();

  cDiffS->cd(1);
  leg1->Draw();

  cDiffS->cd(2);
  leg2->Draw();

  cChi2->cd();
  leg3->Draw();

  TFile* fout=new TFile("ComparisonRawYield.root","RECREATE");
  fout->cd();
  cDiffS->Write();
  cChi2->Write();
  fout->Close();
}
