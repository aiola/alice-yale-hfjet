//
//  Execute with:
//  gSystem->SetIncludePath("-I. -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -I$ALICE_ROOT/ANALYSIS/macros -I$ROOTSYS/include");
//  gROOT->LoadMacro("AliDJetRawYieldUncertainty.cxx++")
//  .L ExtractDJetRawYieldUncertainty.C
//  EvaluateBinPerBinUncertainty(...) //to be done for each pT bin in which you have a mass spectrum
//  ExtractDJetRawYieldUncertainty(...) //to build the uncertainty for the various bins of the jet pT spectrum
// 

void SetInputParametersDzero(AliDJetRawYieldUncertainty *interface);
void SetInputParametersDstar(AliDJetRawYieldUncertainty *interface);

void EvaluateBinPerBinUncertainty(
   Int_t specie=AliDJetRawYieldUncertainty::kD0toKpi,  //D-meson decay channel
   Int_t method=AliDJetRawYieldUncertainty::kEffScale,  //yield extraction method
   Double_t ptmin=0.,  //lower pT edge of mass plot
   Double_t ptmax=99., //upper pT edge of mass plot
   Double_t zmin=0.,   //lower z edge
   Double_t zmax=2.    //upper z edge
   )
{

  AliDJetRawYieldUncertainty *interface = new AliDJetRawYieldUncertainty();
  Bool_t flagSpecie = interface->SetDmesonSpecie((AliDJetRawYieldUncertainty::DMesonSpecies)specie);
  if (!flagSpecie) return;
  interface->SetYieldMethod((AliDJetRawYieldUncertainty::YieldMethod)method);
  interface->SetPtBinEdgesForMassPlot(ptmin,ptmax);
  interface->SetZedges(zmin,zmax);

  if (specie==0) {
    SetInputParametersDzero(interface);  // check the names and the values in the method!!
  }
  else if(specie==1) {
    SetInputParametersDstar(interface);  // check the names and the values in the method!!
  }
  else {
    printf("Error in setting the D-meson specie! Exiting...\n");
    return;
  }

  interface->SetDebugLevel(2); //0 = just do the job; 1 = additional printout; 2 = print individual fits

  Bool_t extract = interface->ExtractInputMassPlot();
  if(!extract) {
    printf("Error in extracting the mass plot! Exiting...\n");
    return;
  }
  
  Bool_t multitrial = interface->RunMultiTrial();
  if(!multitrial) {
    printf("Error in running the MultiTrial code! Exiting...\n");
    return;
  }

  interface->ClearObjects();

  return;
}

//________________________________________
void ExtractDJetRawYieldUncertainty(
   Int_t specie=AliDJetRawYieldUncertainty::kD0toKpi,  //D-meson decay channel
   Int_t method=AliDJetRawYieldUncertainty::kEffScale,  //yield extraction method
   Int_t nTrials=10,  	     //only for SB method: number of random trials for each pT(D) bin to build pT(jet) spectrum variations
   Bool_t allowRepet=kFALSE  //only for SB method: allow repetitions in the extraction of trials in a give pT(D) bin
   )
{

  AliDJetRawYieldUncertainty *interface = new AliDJetRawYieldUncertainty();
  Bool_t flagSpecie = interface->SetDmesonSpecie((AliDJetRawYieldUncertainty::DMesonSpecies)specie);
  if(!flagSpecie) return;
  interface->SetYieldMethod((AliDJetRawYieldUncertainty::YieldMethod)(method));
  interface->SetMaxNTrialsForSidebandMethod(nTrials); //only for SB method: number of random trials for each pT(D) bin to build pT(jet) spectrum variations
  interface->SetAllowRepetitionOfTrialExtraction(allowRepet);
  
  if(specie==0) SetInputParametersDzero(interface);  // here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
  else if(specie==1) SetInputParametersDstar(interface);  // here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
  else {printf("Error in setting the D-meson specie! Exiting...\n"); return;}

  interface->SetDebugLevel(2); //0 = just do the job; 1 = additional printout; 2 = print individual fits

  Bool_t evalunc = interface->EvaluateUncertainty();
  if(!evalunc) {
    printf("Error in evaluating the yield uncertainty! Exiting...\n");
    return;
  }

  interface->ClearObjects();

  return;
}

//________________________________________
void ExtractDJetRawYieldUncertainty_FromSB_CoherentTrialChoice(
   Int_t specie=AliDJetRawYieldUncertainty::kD0toKpi,  //D-meson decay channel
   Int_t nTrials=10
   ) //number of variations is fixed (all the variations in the pT(D) bins, which should match among the various pT(D) bins!)
{

  AliDJetRawYieldUncertainty *interface = new AliDJetRawYieldUncertainty();
  Bool_t flagSpecie = interface->SetDmesonSpecie((AliDJetRawYieldUncertainty::DMesonSpecies)specie);
  if(!flagSpecie) return;
  interface->SetYieldMethod(AliDJetRawYieldUncertainty::kSideband);
  interface->SetMaxNTrialsForSidebandMethod(nTrials); 
  
  if(specie==0) SetInputParametersDzero(interface);  // here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
  else if(specie==1) SetInputParametersDstar(interface);  // here most of the configuration is dummy (not used in the evaluation), you need just the files and some bin ranges
  else {printf("Error in setting the D-meson specie! Exiting...\n"); return;}

  interface->SetDebugLevel(2); //0 = just do the job; 1 = additional printout; 2 = print individual fits

  Bool_t evalunc = interface->EvaluateUncertainty_CoherentTrialChoice();
  if(!evalunc) {
    printf("Error in evaluating the yield uncertainty! Exiting...\n");
    return;
  }

  interface->ClearObjects();

  return;
}

//________________________________________
void SetInputParametersDzero(AliDJetRawYieldUncertainty *interface){

  //Dzero cfg
  const Int_t nDbins = 8;
  Double_t ptDbins[nDbins+1] = {3, 4, 5, 6, 7, 8, 10, 16, 30};
  const Int_t nJetbins = 6;
  Double_t ptJetbins[nJetbins+1] = {5, 6, 8, 10, 14, 20, 30}; //used for eff.scale approach, but also in sideband approach to define the bins of the output jet spectrum
  Double_t DMesonEff[nDbins] = {/*0.0118323, 0.02011807,  0.03644752, */0.05664352 ,0.07682878 ,0.08783701, 0.09420746, 0.1047988, 0.1338670, 0.2143196, 0.2574591}; //chopping 0-1, 1-2

  Double_t sigmafixed=0.014;
  Double_t chi2cut=3;
  Bool_t meansigmaVar[6] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE}; //set mean/sigma variations: fixedS, fixedS+15%, fixedS+15%, freeS&M, freeS/fixedM, fixedS&M
  Bool_t bkgVar[8] = {kTRUE,kFALSE,kTRUE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE}; //set bgk variations: exp, lin, pol2, pol3, pol4, pol5, PowLaw, PowLaw*Exp
  Int_t nRebinSteps=1;
  Int_t rebinStep[1]={1};
  Int_t nMinMassSteps=2;
  Double_t minMassStep[2]={1.72,1.74};
  Int_t nMaxMassSteps=2;
  Double_t maxMassStep[2]={2.00,2.03};
  Int_t nStepsBC=2;
  Double_t nSigmasBC[2]={3.5,4.0};
  //WARNING! set nmask value to active mean/sigma*active bkg variations! 
  //And adjust consequently the following matrix (put an entry for each variation, with value: 0=don't consider it, 1=consider it in the final syst eval) 
  Int_t nmask = 12;
  Bool_t mask[12] =    {1,1,   // fixed sigma (Expo, Lin, Pol2, Pol3, Pol4, Pol5, PowLaw, PowLaw*Exp)
			1,1,   // fixed sigma+15%
			1,1,   // fixed sigma-15%
			1,1,   // free sigma, free mean
			1,1,   // free sigma, fixed mean
			1,1};  // fixed mean, fixed sigma

  interface->AddInputFileName("/Volumes/DATA/ALICE/JetResults/Jets_EMC_pp_823_824_825_826/LHC10b/merge/AnalysisResults.root");
  interface->AddInputFileName("/Volumes/DATA/ALICE/JetResults/Jets_EMC_pp_823_824_825_826/LHC10c/merge/AnalysisResults.root");
  interface->AddInputFileName("/Volumes/DATA/ALICE/JetResults/Jets_EMC_pp_823_824_825_826/LHC10d/merge/AnalysisResults.root");
  interface->AddInputFileName("/Volumes/DATA/ALICE/JetResults/Jets_EMC_pp_823_824_825_826/LHC10e/merge/AnalysisResults.root");
  interface->SetInputTreename("AliAnalysisTaskDmesonJets_AnyINT_D0");
  interface->SetInputDBranchname("DmesonJet");
  interface->SetInputJetBranchname("Jet_AKTChargedR040_pt_scheme");

  interface->SetMassEdgesAndBinWidthForMassPlot(1.5664,2.1664,0.006);  
  interface->SetDmesonPtBins(nDbins,ptDbins);
  interface->SetJetPtBins(nJetbins,ptJetbins);
  interface->SetDmesonEfficiency(DMesonEff);

  interface->SetSigmaForSignalRegion(2.); //only for SB method: sigma range of signal region (usually 3 sigma, also 2 is fine if low S/B)
  interface->SetSigmaToFix(sigmafixed);
  interface->SetChi2Cut(chi2cut);
  interface->SetMeanSigmaVariations(meansigmaVar);
  interface->SetBkgVariations(bkgVar);
  interface->SetRebinSteps(nRebinSteps,rebinStep);
  interface->SetMinMassSteps(nMinMassSteps,minMassStep);
  interface->SetMaxMassSteps(nMaxMassSteps,maxMassStep);
  interface->SetSigmaBinCounting(nStepsBC,nSigmasBC);
  interface->SetMaskOfVariations(nmask,mask);

  return;
}


//________________________________________
void SetInputParametersDstar(AliDJetRawYieldUncertainty *interface){

  //Dstar cfg
  Int_t nDbins = 8;
  Double_t ptDbins[9] = {3,4,5,6,7,8,10,12,24};
  Int_t nJetbins = 6;
  Double_t ptJetbins[7] = {4,6,8,10,12,16,24}; //used for eff.scale approach only (for sideband approach, jet bins are hardcoded in the THnSparses)
  Double_t DMesonEff[8] = {0.0274342, 0.06084, 0.104142, 0.164893, 0.209574, 0.288254, 0.316152, 0.372691};

  Double_t sigmafixed=0.0006;
  Double_t chi2cut=3;
  Bool_t meansigmaVar[6] = {kTRUE,kTRUE,kTRUE,kTRUE,kTRUE,kTRUE}; //set mean/sigma variations: fixedS, fixedS+15%, fixedS+15%, freeS&M, freeS/fixedM, fixedS&M
  Bool_t bkgVar[8] = {kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kTRUE,kTRUE}; //set bgk variations: exp, lin, pol2, pol3, pol4, pol5, PowLaw, PowLaw*Exp
  Int_t nRebinSteps=1;
  Int_t rebinStep[1]={1};
  Int_t nMinMassSteps=2;
  Double_t minMassStep[2]={0.140,0.142};
  Int_t nMaxMassSteps=2;
  Double_t maxMassStep[2]={0.158,0.160};
  Int_t nStepsBC=2;
  Double_t nSigmasBC[2]={3.5,4.0};
  //WARNING! set nmask value to active mean/sigma*active bkg variations! 
  //And adjust consequently the following matrix (put an entry for each variation, with value: 0=don't consider it, 1=consider it in the final syst eval) 
  Int_t nmask = 12;
  Bool_t mask[12] =    {1,1,   // fixed sigma (Expo, Lin, Pol2, Pol3, Pol4, Pol5, PowLaw, PowLaw*Exp)
			1,1,   // fixed sigma+15%
			1,1,   // fixed sigma-15%
			1,1,   // free sigma, free mean
			1,1,   // free sigma, fixed mean
			1,1};  // fixed mean, fixed sigma

  interface->SetInputFilename("./AnalysisResults_Djets_pPb.root");
  interface->SetInputDirname("DmesonsForJetCorrelations");
  interface->SetInputListname("histosDStarMBN");
  interface->SetInputObjectname("hsDphiz");

  interface->SetDmesonPtBins(nDbins,ptDbins);
  interface->SetJetPtBins(nJetbins,ptJetbins);
  interface->SetDmesonEfficiency(DMesonEff);

  interface->SetRebinSpectrumIfSBApproach(kTRUE); //kTRUE=rebin the jet spectrum with ptJetbins[] vals, otherwise use the binning from THnSparse projection

  interface->SetSigmaForSignalRegion(3.); //only for SB method: sigma range of signal region (usually 3 sigma, also 2 is fine if low S/B)
  interface->SetSigmaToFix(sigmafixed);
  interface->SetChi2Cut(chi2cut);
  interface->SetMeanSigmaVariations(meansigmaVar);
  interface->SetBkgVariations(bkgVar);
  interface->SetRebinSteps(nRebinSteps,rebinStep);
  interface->SetMinMassSteps(nMinMassSteps,minMassStep);
  interface->SetMaxMassSteps(nMaxMassSteps,maxMassStep);
  interface->SetSigmaBinCounting(nStepsBC,nSigmasBC);
  interface->SetMaskOfVariations(nmask,mask);

  return;
}

