#include <Riostream.h>
#include <TFile.h>
#include <AliRDHFCutsD0toKpi.h>
#include <TClonesArray.h>
#include <TParameter.h>

/*
 * Generates the RDHF cuts for D0 mesons (v2 analysis).
 * Obtained from Andrea Festanti on March 17th, 2017.
 * Period: LHC15o
 * Pass: pass1
 * Centrality: 30-50%
 * Used for the results shown at Quark Matter 2017
 */

AliRDHFCutsD0toKpi* MakeD0toKpiCuts_PbPb2015_Cent30_50_v2_QM17_AFestanti(Bool_t saveFile = kTRUE)
{
  Float_t minc = 30;
  Float_t maxc = 50;

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  RDHFD0toKpi->SetName("D0toKpiCuts");
  RDHFD0toKpi->SetTitle("Cuts for D0 analysis");

  AliESDtrackCuts* esdTrackCuts=new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  //esdTrackCuts->SetMinNClustersITS(4); // default is 5
  //esdTrackCuts->SetMinNClustersTPC(90);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
					 AliESDtrackCuts::kAny);
 // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.4,1.e10);
  esdTrackCuts->SetMaxDCAToVertexXY(1.);
  esdTrackCuts->SetMaxDCAToVertexZ(1.);
  esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0025*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");

  //esdTrackCuts->SetMinNCrossedRowsTPC(80);
  //esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.9);
  //esdTrackCuts->SetTPCNclsF(80);

  RDHFD0toKpi->AddTrackCuts(esdTrackCuts);
  //RDHFD0toKpi->SetSelectCandTrackSPDFirst(kTRUE, 3);

  RDHFD0toKpi->SetTriggerClass("");
  //RDHFD0toKpi->ResetMaskAndEnableSemiCentralTrigger();
  //RDHFD0toKpi->EnableCentralTrigger();
  //RDHFD0toKpi->EnableMBTrigger();
  //RDHFD0toKpi->ResetMaskAndEnableMBTrigger();

  //RDHFD0toKpi->SetTriggerClass("CINT7");
  RDHFD0toKpi->SetTriggerMask(AliVEvent::kINT7);

  const Int_t nvars=11;

  const Int_t nptbins=28;
  Float_t* ptbins;
  ptbins=new Float_t[nptbins+1];
  ptbins[0]=0.;
  ptbins[1]=0.5;
  ptbins[2]=1.;
  ptbins[3]=1.5;
  ptbins[4]=2.;
  ptbins[5]=2.5;
  ptbins[6]=3.;
  ptbins[7]=3.5;
  ptbins[8]=4.;
  ptbins[9]=4.5;
  ptbins[10]=5.;
  ptbins[11]=5.5;
  ptbins[12]=6.;
  ptbins[13]=6.5;
  ptbins[14]=7.;
  ptbins[15]=7.5;
  ptbins[16]=8.;
  ptbins[17]=8.5;
  ptbins[18]=9.;
  ptbins[19]=10.;
  ptbins[20]=12.;
  ptbins[21]=14.;
  ptbins[22]=16.;
  ptbins[23]=20.;
  ptbins[24]=24.;
  ptbins[25]=36.;
  ptbins[26]=50.;
  ptbins[27]=100.;
  ptbins[28]=9999.;

  RDHFD0toKpi->SetPtBins(nptbins+1,ptbins);

  //setting my cut values
    //cuts order
    //       printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
    //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
    //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
    //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
    //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
    //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
    //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
    //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
    //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);
    //     printf("    |cosThetaPointXY| < %f\n",fD0toKpiCuts[9]);
    //     printf("    NormDecayLenghtXY    > %f\n",fD0toKpiCuts[10]);

    const Int_t nvary=3;
    Float_t varyd0xd0[nptbins][nvary]={{-35000.*1E-8,-40000.*1E-8,-5000.*1E-8},/* pt<0.5*/
				       {-35000.*1E-8,-40000.*1E-8,-25000.*1E-8},/* 0.5<pt<1*/
				       {-30000.*1E-8,-36000.*1E-8,-43000.*1E-8},/* 1<pt<1.5 */
				       {-30000.*1E-8,-36000.*1E-8,-43000.*1E-8},/* 1.5<pt<2 */
				       {-28000.*1E-8,-31000.*1E-8,-40000.*1E-8},/* 2<pt<2.5 */
				       {-28000.*1E-8,-31000.*1E-8,-40000.*1E-8},/* 2.5<pt<3 */
				       {-25000.*1E-8,-28000.*1E-8,-36000.*1E-8},/* 3<pt<3.5 */
				       {-25000.*1E-8,-28000.*1E-8,-36000.*1E-8},/* 3.5<pt<4 */
				       {-12000.*1E-8,-18000.*1E-8,-27000.*1E-8},/* 4<pt<4.5 */
				       {-12000.*1E-8,-18000.*1E-8,-27000.*1E-8},/* 4.5<pt<5 */
				       {-5000.*1E-8,-15000.*1E-8,-21000.*1E-8},/* 5<pt<5.5 */
				       {-5000.*1E-8,-15000.*1E-8,-21000.*1E-8},/* 5.5<pt<6 */
				       {-5000.*1E-8,-10000.*1E-8,-14000.*1E-8},/* 6<pt<6.5 */
				       {-5000.*1E-8,-10000.*1E-8,-14000.*1E-8},/* 6.5<pt<7 */
				       {-5000.*1E-8,-10000.*1E-8,-14000.*1E-8},/* 7<pt<7.5 */
				       {-5000.*1E-8,-10000.*1E-8,-14000.*1E-8},/* 7.5<pt<8 */
				       {-0.*1E-8,-5000.*1E-8,-5000.*1E-8},/* 8<pt<8.5 */
				       {-0.*1E-8,-5000.*1E-8,-5000.*1E-8},/* 8.5<pt<9*/
				       {-0.*1E-8,-5000.*1E-8,-5000.*1E-8},/* 9<pt<10*/
				       {-0.*1E-8,-5000.*1E-8,-5000.*1E-8},/* 10<pt<12 */
				       {0.*1E-8,-5000.*1E-8,0.*1E-8},/* 12<pt<14 */
				       {0.*1E-8,-5000.*1E-8,0.*1E-8},/* 14<pt<16 */
				       {0.*1E-8,-5000.*1E-8,0.*1E-8},/* 16<pt<20 */
				       {5000.*1E-8,-5000.*1E-8,0.*1E-8},/* 20<pt<24 */
				       {5000.*1E-8,-5000.*1E-8,0.*1E-8},/* 24<pt<36 */
				       {5000.*1E-8,-5000.*1E-8,0.*1E-8},/* 36<pt<50 */
				       {5000.*1E-8,-5000.*1E-8,0.*1E-8},/* 50<pt<100 */
				       {5000.*1E-8,-5000.*1E-8,0.*1E-8}};/* pt>100 */


    Float_t varyCosPoint[nptbins][nvary]={{0.75,0.80,0.80},	/* 0<pt<0.5 */
					  {0.75,0.80,0.85},	/* 0.5<pt<1*/
					  {0.75,0.80,0.92},	/* 1<pt<1.5 */
					  {0.75,0.80,0.92},	/* 1.5<pt<2 */
					  {0.92,0.94,0.95},	/* 2<pt<2.5 */
					  {0.92,0.94,0.95},	/* 2.5<pt<3 */
					  {0.90,0.92,0.95},	/* 3<pt<3.5 */
					  {0.90,0.92,0.95},	/* 3.5<pt<4 */
					  {0.85,0.90,0.95},	/* 4<pt<4.5 */
					  {0.85,0.90,0.95},	/* 4.5<pt<5 */
					  {0.88,0.88,0.92},	/* 5<pt<5.5 */
					  {0.88,0.88,0.92},	/* 5.5<pt<6 */
					  {0.83,0.85,0.88},	/* 6<pt<6.5 */
					  {0.83,0.85,0.88},	/* 6.5<pt<7 */
					  {0.83,0.85,0.88},	/* 7<pt<7.5 */
					  {0.83,0.85,0.88},	/* 7.5<pt<8 */
					  {0.8,0.83,0.85},/* 8<pt<8.5 */
					  {0.8,0.83,0.85},/* 8.5<pt<9 */
					  {0.8,0.83,0.85},/* 9<pt<10 */
					  {0.8,0.83,0.85},/* 10<pt<12 */
					  {0.78,0.8,0.83},/* 12<pt<14 */
					  {0.78,0.8,0.83},/* 14<pt<16 */
					  {0.78,0.8,0.83},/* 16<pt<20 */
					  {0.76,0.78,0.81},/* 20<pt<24 */
					  {0.75,0.8,0.8},/*24<pt<36 */
					  {0.75,0.8,0.8},/* 36<pt<50 */
					  {0.75,0.8,0.8},/* 50<pt<100 */
					  {0.75,0.8,0.8}};/* pt>100 */


    Float_t d0MeasMinusExpCut[nptbins]={1.,/* 0<pt<0.5*/
					1.,/* 0.5<pt<1*/
					1.,/* 1<pt<1.5*/
					1.,/* 1.5<pt<2*/
					1.,/*2<pt<2.5*/
					1.,/*2.5<pt<3*/
					1.,/*3<pt<3.5*/
					1.,/*3.5<pt<4*/
					1.,/*4<pt<4.5*/
					1.,/*4.5<pt<5*/
					1.,/*5<pt<5.5*/
					1.,/*5.5<pt<6*/
					1.5,/*6<pt<6.5*/
					1.5,/*6.5<pt<7*/
					1.5,/*7<pt<7.5*/
					1.5,/*7.5<pt<8*/
					1.5,/*8<pt<8.5*/
					1.5,/*8.5<pt<9*/
					1.5,/*9<pt<10*/
					1.5,/*10<pt<12*/
					1.5,/*12<pt<14*/
					1.5,/*14<pt<16*/
					1.5,/*16<pt<20*/
					1.5,/*20<pt<24*/
					1.5,/*24<pt<36*/
					1.5,/*36<pt<50*/
					1.5,/*50<pt<100*/
					1.5,/*>100*/};

    Int_t varycuts = 22;
    Int_t vcd0xd0=varycuts/10;
    Int_t vccospoint=varycuts%10;
    Float_t cutsMatrixD0toKpiStand[nptbins][nvars]={{0.400,400.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,varyd0xd0[0][vcd0xd0],varyCosPoint[0][vccospoint],0.99,6.},/* 0<pt<0.5*/
						    {0.400,400.*1E-4,0.8,0.3,0.3,1000.*1E-4,1000.*1E-4,varyd0xd0[1][vcd0xd0],varyCosPoint[1][vccospoint],0.995,6.},/* 0.5<pt<1*/
						    {0.400,300.*1E-4,0.8,0.4,0.4,1000.*1E-4,1000.*1E-4,varyd0xd0[2][vcd0xd0],varyCosPoint[2][vccospoint],0.993,6.},/* 1<pt<1.5 */
						    {0.400,300.*1E-4,0.8,0.4,0.4,1000.*1E-4,1000.*1E-4,varyd0xd0[3][vcd0xd0],varyCosPoint[3][vccospoint],0.993,6.},/* 1.5<pt<2 */
						    {0.400,250.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,varyd0xd0[4][vcd0xd0],varyCosPoint[4][vccospoint],0.991,5.},/* 2<pt<2.5 */
						    {0.400,250.*1E-4,0.8,0.5,0.5,1000.*1E-4,1000.*1E-4,varyd0xd0[5][vcd0xd0],varyCosPoint[5][vccospoint],0.991,5.},/* 2.5<pt<3 */
						    {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[6][vcd0xd0],varyCosPoint[6][vccospoint],0.993,5.},/* 3<pt<3.5 */
						    {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[7][vcd0xd0],varyCosPoint[7][vccospoint],0.993,5.},/* 3.5<pt<4 */
						    {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[8][vcd0xd0],varyCosPoint[8][vccospoint],0.995,4.},/* 4<pt<4.5*/
						    {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[9][vcd0xd0],varyCosPoint[9][vccospoint],0.995,4.},/* 4.5<pt<5*/
						    {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[10][vcd0xd0],varyCosPoint[10][vccospoint],0.995,4.},/* 5.<pt<5.5 */
						    {0.400,250.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[11][vcd0xd0],varyCosPoint[11][vccospoint],0.995,4.},/* 5.5<pt<6 */
						    {0.400,270.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[12][vcd0xd0],varyCosPoint[12][vccospoint],0.994,4.},/* 6.<pt<6.5 */
						    {0.400,270.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[13][vcd0xd0],varyCosPoint[13][vccospoint],0.994,4.},/* 6.5<pt<7 */
						    {0.400,270.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[14][vcd0xd0],varyCosPoint[14][vccospoint],0.994,4.},/* 7.<pt<7.5 */
						    {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[15][vcd0xd0],varyCosPoint[15][vccospoint],0.994,4.},/* 7.5<pt<8 */
						    {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[16][vcd0xd0],varyCosPoint[16][vccospoint],0.994,4.},/* 8.<pt<8.5 */
						    {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[17][vcd0xd0],varyCosPoint[17][vccospoint],0.994,4.},/* 8.5.<pt<9 */
						    {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[18][vcd0xd0],varyCosPoint[18][vccospoint],0.994,4.},/* 9.<pt<10 */
						    {0.400,300.*1E-4,0.8,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[19][vcd0xd0],varyCosPoint[19][vccospoint],0.994,4.},/* 10.<pt<12 */
						    {0.400,350.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[20][vcd0xd0],varyCosPoint[20][vccospoint],0.994,4.},/* 12.<pt<14 */
						    {0.400,350.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[21][vcd0xd0],varyCosPoint[21][vccospoint],0.994,4.},/* 14.<pt<16 */
						    {0.400,350.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[22][vcd0xd0],varyCosPoint[22][vccospoint],0.994,4.},/* 16.<pt<20 */
						    {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[23][vcd0xd0],varyCosPoint[23][vccospoint],0.994,4.},/* 20.<pt<24 */
						    {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[24][vcd0xd0],varyCosPoint[24][vccospoint],0.994,3.},/* 24.<pt<36 */
						    {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[25][vcd0xd0],varyCosPoint[25][vccospoint],0.994,3.},/* 36.<pt<50 */
						    {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[26][vcd0xd0],varyCosPoint[26][vccospoint],0.994,3.},/* 50.<pt<100 */
						    {0.400,400.*1E-4,1.0,0.7,0.7,1000.*1E-4,1000.*1E-4,varyd0xd0[27][vcd0xd0],varyCosPoint[27][vccospoint],0.994,3.}};/* 100.<pt<9999 */


	Float_t **cutsMatrixTransposeStand=new Float_t*[nvars];
	for(Int_t iv=0;iv<nvars;iv++)cutsMatrixTransposeStand[iv]=new Float_t[nptbins];

	for (Int_t ibin=0;ibin<nptbins;ibin++){
		for (Int_t ivar = 0; ivar<nvars; ivar++){
			cutsMatrixTransposeStand[ivar][ibin]=cutsMatrixD0toKpiStand[ibin][ivar];
			//printf("cutsMatrixD0toKpi[%d][%d] = %f\n",ibin, ivar,cutsMatrixD0toKpiStand[ibin][ivar]);
		}
	}


  RDHFD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);
  RDHFD0toKpi->Setd0MeasMinusExpCut(nptbins,d0MeasMinusExpCut);

  Bool_t pidflag=kTRUE;
  RDHFD0toKpi->SetUsePID(pidflag);
  if(pidflag) cout<<"PID is used"<<endl;
  else cout<<"PID is not used"<<endl;

    //pid settings
  AliAODPidHF* pidObj=new AliAODPidHF();
  //pidObj->SetName("pid4D0");
  Int_t mode=1;
  const Int_t nlims=2;
  Double_t plims[nlims]={0.6,0.8}; //TPC limits in momentum [GeV/c]
  Bool_t compat=kTRUE; //effective only for this mode
  Bool_t asym=kTRUE;
  Double_t sigmas[5]={2.,1.,0.,3.,0.}; //to be checked and to be modified with new implementation of setters by Rossella
  pidObj->SetAsym(asym);// if you want to use the asymmetric bands in TPC
  pidObj->SetMatch(mode);
  pidObj->SetPLimit(plims,nlims);
  pidObj->SetSigma(sigmas);
  pidObj->SetCompat(compat);
  pidObj->SetTPC(kTRUE);
  pidObj->SetTOF(kTRUE);

  pidObj->SetPCompatTOF(2.);
  pidObj->SetSigmaForTPCCompat(3.);
  pidObj->SetSigmaForTOFCompat(3.);

  pidObj->SetOldPid(kFALSE);

  RDHFD0toKpi->SetLowPt(kFALSE);
  RDHFD0toKpi->SetUseDefaultPID(kFALSE); //to use the AliAODPidHF

  //RDHFD0toKpi->SetMaximumPforPID(4.);
  RDHFD0toKpi->SetPidHF(pidObj);

  //activate pileup rejection (for pp)
  RDHFD0toKpi->SetOptPileup(AliRDHFCuts::kNoPileupSelection);

  RDHFD0toKpi->SetUseSpecialCuts(kFALSE);
  //Do not recalculate the vertex
  RDHFD0toKpi->SetRemoveDaughtersFromPrim(kFALSE); //activate for pp

  //RDHFD0toKpi->SetUseAOD049(kTRUE);
  RDHFD0toKpi->SetMaxVtxZ(10.);

  //RDHFD0toKpi->SetRemoveTrackletOutliers(kTRUE);
  RDHFD0toKpi->SetCutOnzVertexSPD(2);

  TString cent="";
  //centrality selection (Pb-Pb)

  RDHFD0toKpi->SetMinCentrality(minc);
  RDHFD0toKpi->SetMaxCentrality(maxc);
  cent=Form("%.0f%.0f",minc,maxc);
  RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentV0M); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid

//  RDHFD0toKpi->SetMaxDifferenceTRKV0Centraltity(5.);

  //temporary
  //RDHFD0toKpi->SetFixRefs();

  //RDHFD0toKpi->SetSelectCandTrackSPDFirst(kTRUE, 3.);

  RDHFD0toKpi->PrintAll();

  if (saveFile) {
    TFile* fout = new TFile("RDHFCuts_v2_D0_LHC15o_pass1_Cent2.root", "recreate");
    fout->cd();
    RDHFD0toKpi->Write();
    fout->Close();
  }

  return RDHFD0toKpi;
}
