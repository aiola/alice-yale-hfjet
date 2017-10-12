// ROOT
#include "Riostream.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TParameter.h"

// ALiRoot
#include "AliESDtrackCuts.h"
#include "AliVEvent.h"

// AliPhysics
#include "AliRDHFCutsD0toKpi.h"

/*
 * Generates the RDHF cuts for D0 mesons (RAA analysis).
 * Based on a file provided by Cristina Terrevoli on March 14th, 2017.
 * Period: LHC10b,c,d,e
 * Pass: pass4
 */

enum pp2010pass4_CutVariation {
  kLoosePointingLoosed0d0 = 1,
  kLoosed0d0 = 2,
  kLoosePointing = 3,
  kTopoOnlyNSigma1 = 11,
  kTopoOnlyNSigma2 = 12,
  kTopoOnlyNSigma3 = 13,
  kTopoOnlyNSigma4 = 14
};

void SetCutsLoosePointingLoosed0d0(AliRDHFCutsD0toKpi* RDHFD0toKpi);
void SetCutsLoosed0d0(AliRDHFCutsD0toKpi* RDHFD0toKpi);
void SetCutsLoosePointing(AliRDHFCutsD0toKpi* RDHFD0toKpi);
void SetCutsTopoOnly(AliRDHFCutsD0toKpi* RDHFD0toKpi, Float_t nsigma);

AliRDHFCutsD0toKpi* MakeD0toKpiCuts_pp2010pass4_Variations(pp2010pass4_CutVariation cutIndex, Bool_t pidflag = kTRUE)
{
  AliRDHFCutsD0toKpi* RDHFD0toKpi = new AliRDHFCutsD0toKpi();
  RDHFD0toKpi->SetName("D0toKpiCuts");
  RDHFD0toKpi->SetTitle("Cuts for D0 analysis");

  RDHFD0toKpi->SetSelectCandTrackSPDFirst(kTRUE,3.);
  // PILE UP REJECTION
  RDHFD0toKpi->SetOptPileup(AliRDHFCuts::kRejectPileupEvent);

  // EVENT CUTS
  RDHFD0toKpi->SetMinVtxContr(1);
  // MAX Z-VERTEX CUT
  RDHFD0toKpi->SetMaxVtxZ(10.);

  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts();
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  //default
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetEtaRange(-0.8,0.8);
  esdTrackCuts->SetMinNClustersTPC(70);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  // default is kBoth, otherwise kAny
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetPtRange(0.3,1.e10);
  esdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  RDHFD0toKpi->AddTrackCuts(esdTrackCuts);

  switch(cutIndex) {
  case kLoosePointingLoosed0d0:
    SetCutsLoosePointingLoosed0d0(RDHFD0toKpi);
    break;
  case kLoosed0d0:
    SetCutsLoosed0d0(RDHFD0toKpi);
    break;
  case kLoosePointing:
    SetCutsLoosePointing(RDHFD0toKpi);
    break;
  case kTopoOnlyNSigma1:
    SetCutsTopoOnly(RDHFD0toKpi, 1);
    break;
  case kTopoOnlyNSigma2:
    SetCutsTopoOnly(RDHFD0toKpi, 2);
    break;
  case kTopoOnlyNSigma3:
    SetCutsTopoOnly(RDHFD0toKpi, 3);
    break;
  case kTopoOnlyNSigma4:
    SetCutsTopoOnly(RDHFD0toKpi, 4);
    break;
  default:
    Printf("Error: cut index %d not defined!", cutIndex);
    return 0;
  }

  //Event selection
  RDHFD0toKpi->SetUsePhysicsSelection(kTRUE);
  //Trigger mask
  RDHFD0toKpi->SetTriggerMask(0);
  RDHFD0toKpi->SetTriggerMask(AliVEvent::kAnyINT);

  RDHFD0toKpi->SetUseSpecialCuts(kTRUE);
  RDHFD0toKpi->SetRemoveDaughtersFromPrim(kTRUE); //Do not recalculate the vertex

  RDHFD0toKpi->SetUsePID(pidflag);
  if (pidflag) {
    Printf("PID is used");
    //pid settings
    AliAODPidHF* pidObj = new AliAODPidHF();
    Int_t mode = 1;
    const Int_t nlims=2;
    Double_t plims[nlims] = {0.6,0.8}; //TPC limits in momentum [GeV/c]
    Bool_t compat = kTRUE; //effective only for this mode
    Bool_t asym = kTRUE;
    Double_t sigmas[5] = {2.,1.,0.,3.,0.}; //to be checked and to be modified with new implementation of setters by Rossella
    pidObj->SetAsym(asym);// if you want to use the asymmetric bands in TPC
    pidObj->SetMatch(mode);
    pidObj->SetPLimit(plims,nlims);
    pidObj->SetSigma(sigmas);
    pidObj->SetCompat(compat);
    pidObj->SetTPC(kTRUE);
    pidObj->SetTOF(kTRUE);
    pidObj->SetPCompatTOF(1.5);
    pidObj->SetSigmaForTPCCompat(3.);
    pidObj->SetSigmaForTOFCompat(3.);
    pidObj->SetOldPid(kFALSE);
    RDHFD0toKpi->SetPidHF(pidObj);
    RDHFD0toKpi->SetUseDefaultPID(kFALSE);
  }
  else {
    Printf("PID is not used");
  }

  RDHFD0toKpi->SetLowPt(kFALSE);
  RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentOff); //kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid

  RDHFD0toKpi->PrintAll();

  return RDHFD0toKpi;
}

void SetCutsLoosePointingLoosed0d0(AliRDHFCutsD0toKpi* RDHFD0toKpi)
{
  const Int_t nvars = 11;
  const Int_t nptbins = 16;
  Float_t ptbins[nptbins+1] = {0., 0.5, 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 20., 24., 36., 9999.};

  RDHFD0toKpi->SetPtBins(nptbins,ptbins);
  RDHFD0toKpi->SetGlobalIndex(nvars,nptbins);

  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]=  {
      //m     dca        cost*  ptk  ptpi d0k         d0pi        d0d0          cosp   cosxy  normdxy
      {0.400, 350.*1E-4, 0.8,   0.5, 0.5, 1000.*1E-4, 1000.*1E-4,  15000. *1E-8, 0.70,  0.,    0.},/* pt<0.5*/
      {0.400, 350.*1E-4, 0.8,   0.5, 0.5, 1000.*1E-4, 1000.*1E-4,  15000. *1E-8, 0.70,  0.,    0.},/* 0.5<pt<1*/
      {0.400, 300.*1E-4, 0.8,   0.4, 0.4, 1000.*1E-4, 1000.*1E-4,  -5000. *1E-8, 0.70,  0.,    0.},/* 1<pt<2 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,      0. *1E-8, 0.80,  0.,    0.},/* 2<pt<3 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,   8000. *1E-8, 0.80,  0.,    0.},/* 3<pt<4 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  12000. *1E-8, 0.75,  0.,    0.},/* 4<pt<5 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  12000. *1E-8, 0.75,  0.,    0.},/* 5<pt<6 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  12000. *1E-8, 0.75,  0.,    0.},/* 6<pt<7 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  13000. *1E-8, 0.75,  0.,    0.},/* 7<pt<8 */
      {0.400, 300.*1E-4, 0.9,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  15000. *1E-8, 0.75,  0.,    0.},/* 8<pt<10 */
      {0.400, 300.*1E-4, 0.9,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  15000. *1E-8, 0.75,  0.,    0.},/* 10<pt<12 */
      {0.400, 300.*1E-4, 1.0,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  30000. *1E-8, 0.75,  0.,    0.},/* 12<pt<16 */
      {0.400, 300.*1E-4, 1.0,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.75,  0.,    0.},/* 16<pt<20 */
      {0.400, 300.*1E-4, 1.0,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.75,  0.,    0.},/* 20<pt<24 */
      {0.400, 300.*1E-4, 1.0,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.75,  0.,    0.},/* 24<pt<36 */
      {0.400, 300.*1E-4, 1.0,   0.6, 0.6, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.70,  0.,    0.}};/* pt>36 */

  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand = new Float_t*[nvars];
  for(Int_t iv = 0; iv < nvars; iv++) cutsMatrixTransposeStand[iv] = new Float_t[nptbins];

  for (Int_t ibin = 0; ibin < nptbins; ibin++) {
    for (Int_t ivar = 0; ivar < nvars; ivar++) {
      cutsMatrixTransposeStand[ivar][ibin] = cutsMatrixD0toKpiStand[ibin][ivar];
    }
  }
  //new cut
  Float_t cutvalTopo[nptbins] = {2., 2., 2., 2., 2., 2., 3., 3., 3., 3., 3., 4., 3., 2.5, 2.5, 2.};
  RDHFD0toKpi->Setd0MeasMinusExpCut(nptbins,cutvalTopo);
  RDHFD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);

  for (Int_t iv = 0; iv < nvars; iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand = nullptr;
}

void SetCutsLoosePointing(AliRDHFCutsD0toKpi* RDHFD0toKpi)
{
  const Int_t nvars = 11;
  const Int_t nptbins = 16;
  Float_t ptbins[nptbins+1] = {0., 0.5, 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 20., 24., 36., 9999.};

  RDHFD0toKpi->SetPtBins(nptbins,ptbins);
  RDHFD0toKpi->SetGlobalIndex(nvars,nptbins);

  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]=  {
      //m    dca      cost*  ptk ptpi  d0k          d0pi       d0d0          cosp  cosxy normdxy
      {0.400,350.*1E-4, 0.8, 0.5, 0.5, 1000.*1E-4, 1000.*1E-4,  -5000. *1E-8, 0.70,  0.,0.},/* pt<0.5*/
      {0.400,350.*1E-4, 0.8, 0.5, 0.5, 1000.*1E-4, 1000.*1E-4,  -5000. *1E-8, 0.70,  0.,0.},/* 0.5<pt<1*/
      {0.400,300.*1E-4, 0.8, 0.4, 0.4, 1000.*1E-4, 1000.*1E-4, -25000. *1E-8, 0.70,  0.,0.},/* 1<pt<2 */
      {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -20000. *1E-8, 0.80,  0.,0.},/* 2<pt<3 *///d0d0 e cosp
      {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, -12000. *1E-8, 0.80,  0.,0.},/* 3<pt<4 */
      {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  -8000. *1E-8, 0.75,  0.,0.},/* 4<pt<5 */
      {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  -8000. *1E-8, 0.75,  0.,0.},/* 5<pt<6 */
      {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  -8000. *1E-8, 0.75,  0.,0.},/* 6<pt<7 */
      {0.400,300.*1E-4, 0.8, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  -7000. *1E-8, 0.75,  0.,0.},/* 7<pt<8 */
      {0.400,300.*1E-4, 0.9, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  -5000. *1E-8, 0.75,  0.,0.},/* 8<pt<10 */
      {0.400,300.*1E-4, 0.9, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  -5000. *1E-8, 0.75,  0.,0.},/* 10<pt<12 */
      {0.400,300.*1E-4, 1.0, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  10000. *1E-8, 0.75,  0.,0.},/* 12<pt<16 */
      {0.400,300.*1E-4, 1.0, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.75,  0.,0.},/* 16<pt<20 */
      {0.400,300.*1E-4, 1.0, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.75,  0.,0.},/* 20<pt<24 */
      {0.400,300.*1E-4, 1.0, 0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.75,  0.,1.},/* 24<pt<36 */
      {0.400,300.*1E-4, 1.0, 0.6, 0.6, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.70,  0.,0.}};/* pt>24 */

  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand = new Float_t*[nvars];
  for(Int_t iv = 0; iv < nvars; iv++) cutsMatrixTransposeStand[iv] = new Float_t[nptbins];

  for (Int_t ibin = 0; ibin < nptbins; ibin++) {
    for (Int_t ivar = 0; ivar < nvars; ivar++) {
      cutsMatrixTransposeStand[ivar][ibin] = cutsMatrixD0toKpiStand[ibin][ivar];
    }
  }
  //new cut
  Float_t cutvalTopo[nptbins] = {2., 2., 2., 2., 2., 2., 3., 3., 3., 3., 3., 4., 3., 2.5, 2.5, 2.};
  RDHFD0toKpi->Setd0MeasMinusExpCut(nptbins,cutvalTopo);
  RDHFD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);

  for (Int_t iv = 0; iv < nvars; iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand = nullptr;
}

void SetCutsLoosed0d0(AliRDHFCutsD0toKpi* RDHFD0toKpi)
{
  const Int_t nvars = 11;
  const Int_t nptbins = 16;
  Float_t ptbins[nptbins+1] = {0., 0.5, 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 20., 24., 36., 9999.};

  RDHFD0toKpi->SetPtBins(nptbins,ptbins);
  RDHFD0toKpi->SetGlobalIndex(nvars,nptbins);

  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]=  {
      //m     dca        cost*  ptk  ptpi d0k         d0pi        d0d0          cosp   cosxy  normdxy
      {0.400, 350.*1E-4, 0.8,   0.5, 0.5, 1000.*1E-4, 1000.*1E-4,  15000. *1E-8, 0.80,  0.,    0.},/* pt<0.5*/
      {0.400, 350.*1E-4, 0.8,   0.5, 0.5, 1000.*1E-4, 1000.*1E-4,  15000. *1E-8, 0.80,  0.,    0.},/* 0.5<pt<1*/
      {0.400, 300.*1E-4, 0.8,   0.4, 0.4, 1000.*1E-4, 1000.*1E-4,  -5000. *1E-8, 0.80,  0.,    0.},/* 1<pt<2 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,      0. *1E-8, 0.90,  0.,    0.},/* 2<pt<3 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,   8000. *1E-8, 0.90,  0.,    0.},/* 3<pt<4 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  12000. *1E-8, 0.85,  0.,    0.},/* 4<pt<5 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  12000. *1E-8, 0.85,  0.,    0.},/* 5<pt<6 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  12000. *1E-8, 0.85,  0.,    0.},/* 6<pt<7 */
      {0.400, 300.*1E-4, 0.8,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  13000. *1E-8, 0.85,  0.,    0.},/* 7<pt<8 */
      {0.400, 300.*1E-4, 0.9,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  15000. *1E-8, 0.85,  0.,    0.},/* 8<pt<10 */
      {0.400, 300.*1E-4, 0.9,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  15000. *1E-8, 0.85,  0.,    0.},/* 10<pt<12 */
      {0.400, 300.*1E-4, 1.0,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4,  30000. *1E-8, 0.85,  0.,    0.},/* 12<pt<16 */
      {0.400, 300.*1E-4, 1.0,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.85,  0.,    0.},/* 16<pt<20 */
      {0.400, 300.*1E-4, 1.0,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.85,  0.,    0.},/* 20<pt<24 */
      {0.400, 300.*1E-4, 1.0,   0.7, 0.7, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.85,  0.,    0.},/* 24<pt<36 */
      {0.400, 300.*1E-4, 1.0,   0.6, 0.6, 1000.*1E-4, 1000.*1E-4, 999999. *1E-8, 0.80,  0.,    0.}};/* pt>36 */

  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand = new Float_t*[nvars];
  for(Int_t iv = 0; iv < nvars; iv++) cutsMatrixTransposeStand[iv] = new Float_t[nptbins];

  for (Int_t ibin = 0; ibin < nptbins; ibin++) {
    for (Int_t ivar = 0; ivar < nvars; ivar++) {
      cutsMatrixTransposeStand[ivar][ibin] = cutsMatrixD0toKpiStand[ibin][ivar];
    }
  }
  //new cut
  Float_t cutvalTopo[nptbins] = {2., 2., 2., 2., 2., 2., 3., 3., 3., 3., 3., 4., 3., 2.5, 2.5, 2.};
  RDHFD0toKpi->Setd0MeasMinusExpCut(nptbins,cutvalTopo);
  RDHFD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);

  for (Int_t iv = 0; iv < nvars; iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand = nullptr;
}

void SetCutsTopoOnly(AliRDHFCutsD0toKpi* RDHFD0toKpi, Float_t nsigma)
{
  const Int_t nvars = 11;
  const Int_t nptbins = 16;
  Float_t ptbins[nptbins+1] = {0., 0.5, 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 20., 24., 36., 9999.};

  RDHFD0toKpi->SetPtBins(nptbins,ptbins);
  RDHFD0toKpi->SetGlobalIndex(nvars,nptbins);

  Float_t cutsMatrixD0toKpiStand[nptbins][nvars]=  {
      //m     dca        cost*  ptk  ptpi d0k         d0pi        d0d0          cosp   cosxy  normdxy
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* pt<0.5*/
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 0.5<pt<1*/
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 1<pt<2 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 2<pt<3 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 3<pt<4 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 4<pt<5 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 5<pt<6 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 6<pt<7 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 7<pt<8 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 8<pt<10 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 10<pt<12 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 12<pt<16 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 16<pt<20 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 20<pt<24 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.},/* 24<pt<36 */
      {0.3,     999999., 1.1,    0., 0.,     999999.,    999999.,       999999.,   0.,  -1,    0.}};/* pt>36 */

  //CREATE TRANSPOSE MATRIX...REVERSE INDICES as required by AliRDHFCuts
  Float_t **cutsMatrixTransposeStand = new Float_t*[nvars];
  for(Int_t iv = 0; iv < nvars; iv++) cutsMatrixTransposeStand[iv] = new Float_t[nptbins];

  for (Int_t ibin = 0; ibin < nptbins; ibin++) {
    for (Int_t ivar = 0; ivar < nvars; ivar++) {
      cutsMatrixTransposeStand[ivar][ibin] = cutsMatrixD0toKpiStand[ibin][ivar];
    }
  }
  //new cut
  Float_t cutvalTopo[nptbins] = {nsigma, nsigma, nsigma, nsigma, nsigma, nsigma, nsigma, nsigma, nsigma, nsigma, nsigma, nsigma, nsigma, nsigma, nsigma, nsigma};
  RDHFD0toKpi->Setd0MeasMinusExpCut(nptbins,cutvalTopo);
  RDHFD0toKpi->SetCuts(nvars,nptbins,cutsMatrixTransposeStand);

  for (Int_t iv = 0; iv < nvars; iv++) delete [] cutsMatrixTransposeStand[iv];
  delete [] cutsMatrixTransposeStand;
  cutsMatrixTransposeStand = nullptr;
}
