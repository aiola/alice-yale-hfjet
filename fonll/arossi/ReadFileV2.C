#include <fstream>
void ReadFileAndGetNtupleAndHisto(TString strfile="Predictionbquark14TeVMidrapidity.dat",Int_t conversion=2/*1->no scaling, 2: pb->mub */){
  TNtupleD *mynt=new TNtupleD("ntData","ntData","pt:centr:min:max:minsc:maxsc:minmass:maxmass:min_pdf:max_pdf:fr0505:fr22:fr21:fr12:fr105:fr051");
  //Double_t pt,centr,min,max,minsc,maxsc,minmass,maxmass,fr0505,fr22,fr21,fr12,fr105,fr051;


  FILE *f=fopen(strfile.Data(),"r");

  std::ifstream infile(strfile.Data());
  mynt->ReadStream(infile);
  mynt->Print();
  
  TString strfileout=strfile;
  strfileout.ReplaceAll(".dat",".root");
  strfileout.ReplaceAll(".txt",".root");
  Double_t scalingFactor=1.;
  if(conversion==2){
    scalingFactor=1.e-6;
    strfileout.ReplaceAll(".root","microbarn.root");
  }

  Int_t nvars=mynt->GetNvar();
  Double_t *var=new Double_t[nvars];// going this way was just an exercise

  TGraphErrors **gr=new TGraphErrors*[nvars-1];// one variable is pt
  TH1D **hist=new TH1D*[nvars-1];
  TObjArray *obja=mynt->GetListOfBranches();
//   mynt->SetBranchAddress("pt",&pt);
//   mynt->SetBranchAddress("");
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();
//   mynt->SetBranchAddress();

  for(Int_t j=0;j<nvars;j++){// assume first variable is pt
    mynt->SetBranchAddress(obja->At(j)->GetName(),&var[j]);
  }

  // Now get histo granularity
  Int_t nent=mynt->GetEntries();
  mynt->GetEntry(0);
  Double_t ptA=var[0];
  mynt->GetEntry(1);
  Double_t ptB=var[0];
  Double_t dpt=ptB-ptA;
  mynt->GetEntry(nent-1);
  Double_t ptMax=var[0];
  Double_t minptHist=0.;
  if(ptA/dpt-(Int_t)(ptA/dpt)<1.e-6)minptHist=-dpt/2.;
  else if((ptA-dpt/2.)/dpt-(Int_t)((ptA-dpt/2.)/dpt)<1.e-6)minptHist=0;
  else {
    printf("pt binning issues \n");
    return;
  }
  Double_t *maxPred=new Double_t[nent];
  Double_t *minPred=new Double_t[nent];

  if(TMath::Abs(ptMax-((Double_t)nent-1.)*dpt-ptA)>0.001){
    printf("Not fixed pt bin width: ptmin=%f, ptmax=%f, nbins=%d, dpt=%f \n",ptA,ptMax,nent,dpt);
    return;
  }

  // Create histograms
  TH1D *hMaxPred=new TH1D("histMaxPred","histMaxPred",(Int_t)((ptMax-minptHist)/dpt)+1,minptHist,ptMax+dpt/2.);
  TH1D *hMinPred=new TH1D("histMinPred","histMinPred",(Int_t)((ptMax-minptHist)/dpt)+1,minptHist,ptMax+dpt/2.);  
  TGraphAsymmErrors *grCentMinMax=new TGraphAsymmErrors();
  grCentMinMax->SetName("grCentMinMaxPred");
  if(conversion==2){
    hMaxPred->SetYTitle("#frac{d#sigma}{dp_{T}} (#mub/(GeV/c))");
    hMinPred->SetYTitle("#frac{d#sigma}{dp_{T}} (#mub/(GeV/c))");
    grCentMinMax->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}} (#mub/(GeV/c))");
    
  }
  else {
    hMaxPred->SetYTitle("#frac{d#sigma}{dp_{T}} (pb/(GeV/c))");
    hMinPred->SetYTitle("#frac{d#sigma}{dp_{T}} (pb/(GeV/c))");
    grCentMinMax->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}} (pb/(GeV/c))");
  }
  hMaxPred->SetXTitle("p_{T} (GeV/c)");
  hMinPred->SetXTitle("p_{T} (GeV/c)");
  grCentMinMax->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  

  for(Int_t j=1;j<nvars;j++){// assume first variable is pt

    gr[j]=new TGraphErrors();
    gr[j]->SetName(Form("gr%s",obja->At(j)->GetName()));
    hist[j]=new TH1D(Form("hist%s",obja->At(j)->GetName()),Form("hist%s",obja->At(j)->GetName()),(Int_t)((ptMax-minptHist)/dpt)+1,minptHist,ptMax+dpt/2.);//always start the hist range from 0
    if(conversion==2){
      hist[j]->SetYTitle("#frac{d#sigma}{dp_{T}} (#mub/(GeV/c))");
      gr[j]->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}} (#mub/(GeV/c))");
	  
    }
    else {
      hist[j]->SetYTitle("#frac{d#sigma}{dp_{T}} (pb/(GeV/c))");
      gr[j]->GetYaxis()->SetTitle("#frac{d#sigma}{dp_{T}} (pb/(GeV/c))");
    }
    hist[j]->SetXTitle("p_{T} (GeV/c)");
    gr[j]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  }
  
  for(Int_t j=0;j<nent;j++){
    mynt->GetEntry(j);
    maxPred[j]=var[1];// initialize to central prediction
    minPred[j]=var[1];
    Int_t binh=hMaxPred->FindBin(var[0]);
    for(Int_t i=1;i<nvars;i++){// assume first variable is pt
      gr[i]->SetPoint(j,var[0],var[i]*scalingFactor);
      gr[i]->SetPointError(j,dpt/2.,var[i]*scalingFactor*1.e-6);// the error is a fake number
      hist[i]->SetBinContent(binh,var[i]*scalingFactor);
      hist[i]->SetBinError(binh,var[i]*scalingFactor*1.e-6);
      if(var[i]>maxPred[j])maxPred[j]=var[i];
      if(var[i]<minPred[j])minPred[j]=var[i];
    }
    
    maxPred[j]*=scalingFactor;
    minPred[j]*=scalingFactor;
    hMaxPred->SetBinContent(binh,maxPred[j]);
    hMaxPred->SetBinError(binh,maxPred[j]*1.e-6);
    hMinPred->SetBinContent(binh,minPred[j]);
    hMinPred->SetBinError(binh,minPred[j]*1.e-6);
    grCentMinMax->SetPoint(j,var[0],var[1]*scalingFactor);
    grCentMinMax->SetPointError(j,dpt/2.,dpt/2.,var[1]*scalingFactor-minPred[j],maxPred[j]-var[1]*scalingFactor);// the error is a fake number
 

  }
  
  

  TFile *fout=new TFile(strfileout.Data(),"RECREATE");
  mynt->Write();
  hMaxPred->Write();
  hMinPred->Write();
  grCentMinMax->Write();
  for(Int_t j=1;j<nvars;j++){

    gr[j]->Write();
    hist[j]->Write();
  }

  fout->Close();
}
void ReadFile(TString strfile="Predictionbquark14TeVMidrapidity.dat",const int npoints=400){

  Float_t pt[npoints],centr[npoints],min[npoints],max[npoints],minsc[npoints],maxsc[npoints],minmass[npoints],maxmass[npoints],fr0505[npoints],fr22[npoints],fr21[npoints],fr12[npoints],fr105[npoints],fr051[npoints];

  FILE *f=fopen(strfile.Data(),"r");

  std::ifstream infile(strfile.Data());
  Int_t last=0;
  std::string line;
  while (std::getline(infile, line))
    {
      std::size_t found = line.find_first_of("#");
      if(found==std::string::npos){
	sscanf(line.data(),"%f %f %f %f %f %f %f %f",&pt[last],&centr[last],&min[last],&max[last],&minsc[last],&maxsc[last],&minmass[last],&maxmass[last]);// need to go for a loop over the string to get also the others (sscanf has a maximum of 12 arguments)
	printf("pt=%f, centr =%f, fr051=%f\n",pt[last],centr[last],maxmass[last]);
	last++;
      }
      else printf("Reading Introductory line \n");
   }
  

//   char *strInt[255];
//   TString str;
//   while (1){
//     fscanf(f,"%s",strInt);
//     if(strInt[0]=="#")printf("Reading Introductory line \n");
//     //.Contains("#"))printf("Reading Introductory line \n");
//     else break;
//   }
//   printf("Introdcution finished \n");



}


void ReadFile2(TString strfile="Predictionbquark14TeVMidrapidity.dat"){// does not work...

  Float_t pt,centr,min,max,min_sc,max_sc,min_mass,max_mass,fr0505,fr22,fr21,fr12,fr105,fr051;
  FILE *f=fopen(strfile.Data(),"r");
  char *strInt[255];
  TString str;
  while (1){
    fscanf(f,"%s",strInt);
    if(strInt[0]=="#")printf("Reading Introductory line \n");
    //.Contains("#"))printf("Reading Introductory line \n");
    else break;
  }
  printf("Introdcution finished \n");

  while (    fscanf(f,"%f %f %f %f",&pt,&centr,&min,&max)==4){
    printf("pt=%f, centr =%f, fr051=%f",pt,centr,max);
    
  }
}
