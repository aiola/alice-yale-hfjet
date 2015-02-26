class JetUnfoldPbPb;
class JetUnfold;
class JetUnfoldParams;
class JetUnfoldPtRanges;

const char *dir = "";

JetUnfoldPbPb *gUnfold = 0;
JetUnfold *gUnfoldBkg = 0;

void AddUnfoldingParams(JetUnfold* unfold, UInt_t meth, Int_t minReg = 1, Int_t maxReg = 8, 
			const char* prior = "Measured", Int_t steep1=0, Int_t steep2=-999, Int_t steep3=-999,
			const char* prior2 = 0);

void LoadTestUnfold()
{
  // RooUnfold
  gSystem->Load("libRooUnfold");

  // AliUnfolding
  gSystem->Load("libVMC");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  // My unfolding classes
  gROOT->LoadMacro("JetUnfold.cxx+g");
  gROOT->LoadMacro("JetUnfoldPbPb.cxx+g");

  gUnfold = new JetUnfoldPbPb();
}

void AddUnfoldingParams(JetUnfold* unfold, UInt_t meth, Int_t minReg, Int_t maxReg, 
			const char* prior, Int_t steep1, Int_t steep2, Int_t steep3,
			const char* prior2)
{
  JetUnfoldParams params;
  params.SetRegBase(1e-6);
 
  params.SetUnfoldingMethod(meth);
  params.SetPriorName(prior);
  params.SetPriorSteepness(steep1);

  for (Int_t i = minReg; i <=maxReg; i++) {
    params.SetRegParam(i);
    unfold->AddUnfoldParams(params);
  }

  if (steep2 != -999) {
    params.SetUnfoldingMethod(meth);
    params.SetPriorName(prior);
    params.SetPriorSteepness(steep2);
    
    for (Int_t i = minReg; i <=maxReg; i++) {
      params.SetRegParam(i);
      unfold->AddUnfoldParams(params);
    }
  }

  if (steep3 != -999) {
    params.SetUnfoldingMethod(meth);
    params.SetPriorName(prior);
    params.SetPriorSteepness(steep3);
    
    for (Int_t i = minReg; i <=maxReg; i++) {
      params.SetRegParam(i);
      unfold->AddUnfoldParams(params);
    }
  }

  if (prior2 != 0) {
    params.SetUnfoldingMethod(meth);
    params.SetPriorName(prior2);
    params.SetPriorSteepness(steep1);
    
    for (Int_t i = minReg; i <=maxReg; i++) {
      params.SetRegParam(i);
      unfold->AddUnfoldParams(params);
    }
  }
}

void SetUnfold()
{
  const char* suffix = "_bias0_stat1";

  gUnfold->SetDebugLevel(JetUnfold::kToyDebug1);

  gUnfold->SetTwoStepMethod(JetUnfoldPbPb::kNoTwoStep);
  //gUnfold->SetTwoStepMethod(JetUnfoldPbPb::kNormal);

  gUnfold->SetShowTruth(kTRUE);
  gUnfold->SetEMCalAcceptance(0.2);

  gUnfold->SetPlotPath(Form("TestPlots%s_ptrange", suffix));
  gUnfold->SetFileFormat("pdf");

  TString workPath(Form("./%s",dir));
  gUnfold->SetWorkingPath(workPath);

  gUnfold->SetInputFileName(Form("UnfoldingTestData%s.root", suffix));
  gUnfold->SetOutputFileName(Form("UnfoldingTestOutput%s_ptrange.root", suffix));
 
  const Int_t    nBinsTrue                = 15;
  Double_t       binsTrue[nBinsTrue+1]    = {0,10,20,30,40,50,60,70,80,90,100,110,120,150,170,200};

  const Int_t    nBinsInter               = 31;
  Double_t       binsInter[nBinsInter+1]  = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,150,170,200};
      
  const Int_t    nBinsMeas                = 30;
  Double_t       binsMeas[nBinsMeas+1]    = {-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120};

  gUnfold->SetTrueBins(nBinsTrue,binsTrue);
  gUnfold->SetInterBins(nBinsInter,binsInter);
  gUnfold->SetMeasBins(nBinsMeas-12,binsMeas+12);

  gUnfold->SetFinalLimits(0,119.9);
  gUnfold->SetDrawLimits(-30,119.9);

  gUnfoldBkg = new JetUnfold(*gUnfold);
  gUnfoldBkg->SetInputResponseName("Background");
  gUnfoldBkg->SetTruthName("InterTruth");
  gUnfoldBkg->SetOutputFileName(Form("UnfoldingTestOutput%s_bkg.root", suffix));

  AddUnfoldingParams(gUnfold,JetUnfold::kBayesian,1,8,"MCTruth",-1);
  AddUnfoldingParams(gUnfold,JetUnfold::kSVD,1,6,"MCTruth",-1,1,0);
  AddUnfoldingParams(gUnfold,JetUnfold::kChi2,6,15,"MCTruth",-1);
 
  gUnfold->SetDefaultParams(4, JetUnfold::kBayesian); // reg=5
  gUnfold->SetDefaultParams(10);//, JetUnfold::kSVD); // reg=5
  gUnfold->SetDefaultParams(27,JetUnfold::kChi2);     // reg=2e-2

  AddUnfoldingParams(gUnfoldBkg,JetUnfold::kBayesian,1,8,"MCInterTruth",-1);
  AddUnfoldingParams(gUnfoldBkg,JetUnfold::kSVD,1,6,"MCInterTruth",-1,1,0);
  AddUnfoldingParams(gUnfoldBkg,JetUnfold::kChi2,6,15,"MCInterTruth",-1);
 
  gUnfoldBkg->SetDefaultParams(4, JetUnfold::kBayesian); // reg=5
  gUnfoldBkg->SetDefaultParams(9);//, JetUnfold::kSVD); // reg=5
  gUnfoldBkg->SetDefaultParams(27,JetUnfold::kChi2);     // reg=2e-2

}
