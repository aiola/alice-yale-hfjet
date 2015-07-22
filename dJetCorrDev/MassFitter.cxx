// This class provides tools to fit invariant mass distributions
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include "MassFitter.h"
#include <TF1.h>
#include <TH1.h>
#include <TDatabasePDG.h>
#include <TFitter.h>

const Double_t MassFitter::fgkEpsilon = 1e-15;

//____________________________________________________________________________________
MassFitter::MassFitter() :
  TNamed("MassFitter", "MassFitter"),
  fMassFitTypeSig(kGaus),
  fMassFitTypeBkg(kExpo),
  fMean(0.),
  fMeanError(0.),
  fWidth(0.),
  fWidthError(0.),
  fSignal(0.),
  fSignalError(0.),
  fBackground(0.),
  fBackgroundError(0.),
  fDisableBkg(kFALSE),
  fDisableSig(kFALSE),
  fNParSig(0),
  fNParBkg(0),
  fMinMass(0.),
  fMaxMass(0.),
  fFunction(0),
  fFunctionBkg(0),
  fHistogram(0),
  fPionMass(TDatabasePDG::Instance()->GetParticle(211)->Mass())
{
  // Default constructor.

  Reset();
}

//____________________________________________________________________________________
MassFitter::MassFitter(const char* name, EMassFitTypeSig ts, EMassFitTypeBkg tb, Double_t minMass, Double_t maxMass) :
  TNamed(name, name),
  fMassFitTypeSig(ts),
  fMassFitTypeBkg(tb),
  fMean(0.),
  fMeanError(0.),
  fWidth(0.),
  fWidthError(0.),
  fSignal(0.),
  fSignalError(0.),
  fBackground(0.),
  fBackgroundError(0.),
  fDisableBkg(kFALSE),
  fDisableSig(kFALSE),
  fNParSig(0),
  fNParBkg(0),
  fMinMass(minMass),
  fMaxMass(maxMass),
  fFunction(0),
  fFunctionBkg(0),
  fHistogram(0),
  fPionMass(TDatabasePDG::Instance()->GetParticle(211)->Mass())
{
  // Standard constructor.

  Reset();
}

//____________________________________________________________________________________
void MassFitter::Reset(TH1* histo)
{
  // Reset fitter

  switch (fMassFitTypeSig) {
  case kGaus:
    {
      fNParSig = 3;
      break;
    }
  default:
    {
      Printf("Error: signal fit type %d not recognized! Using Gaussian fit.", fMassFitTypeSig);
      fMassFitTypeSig = kGaus;
      fNParSig = 3;
    }
  }

  switch (fMassFitTypeBkg) {
  case kExpo:
    {
      fNParBkg = 2;
      break;
    }
  case kExpoPower:
    {
      if (fMinMass < fPionMass) {
        Printf("MassFitter: setting min mass to %.3f", fPionMass);
        fMinMass = fPionMass;
      }
      fNParBkg = 2;
      break;
    }
  default:
    {
      Printf("Error: background fit type %d not recognized! Using exponential fit.", fMassFitTypeBkg);
      fMassFitTypeBkg = kExpo;
      fNParBkg = 2;
    }
  }
  
  fHistogram = histo;

  fSignal = 0;
  fSignalError = 0;
  fBackground = 0;
  fBackgroundError = 0;

  fMean = 0;
  fMeanError = 0;
  fWidth = 0;
  fWidthError = 0;
  
  delete fFunction;
  delete fFunctionBkg;

  TString fname(Form("%s_fit", GetName()));
  fFunction = new TF1(fname, this, &MassFitter::FunctionSigBkg, fMinMass, fMaxMass, fNParSig+fNParBkg);
  fFunction->SetLineColor(kBlue+1);
  fFunction->SetLineWidth(2);
    
  TString fnameBkg(Form("%sBkg_fit", GetName()));
  fFunctionBkg = new TF1(fnameBkg, this, &MassFitter::FunctionBkg, fMinMass, fMaxMass, fNParBkg);
  fFunctionBkg->SetLineColor(kRed+1);
  fFunctionBkg->SetLineWidth(2);
}

//____________________________________________________________________________________
TFitResultPtr MassFitter::Fit(TH1* histo, Option_t* opt)
{
  TFitter::SetPrecision(0.1);
  
  if (!fFunction) {
    Reset(histo);
  }
  else {
    fHistogram = histo;
  }

  return Fit(opt);
}

//____________________________________________________________________________________
TFitResultPtr MassFitter::Fit(Option_t* opt)
{
  if (!fHistogram) {
    Printf("Error: no histogram provided!");
    return 0;
  }

  if (!fFunction) {
    Reset(fHistogram);
  }
  
  TFitResultPtr r = fHistogram->Fit(fFunction, opt);

  for (Int_t i = 1; i < fNParBkg; i++) {
    fFunctionBkg->SetParameter(i, fFunction->GetParameter(i));
  }
  fFunctionBkg->SetParameter(0, fFunction->GetParameter(0) - fFunction->GetParameter(fNParBkg));

  fMean = fFunction->GetParameter(fNParBkg+1);
  fMeanError = fFunction->GetParError(fNParBkg+1);
  fWidth = fFunction->GetParameter(fNParBkg+2);
  fWidthError = fFunction->GetParError(fNParBkg+2);
  
  fSignal = fFunction->GetParameter(fNParBkg);
  fSignalError = fFunction->GetParError(fNParBkg);

  Double_t minMass = fMean - fWidth*3.0;
  if (minMass < fMinMass) minMass = fMinMass;
  Double_t maxMass = fMean + fWidth*3.0;
  if (maxMass > fMaxMass) maxMass = fMaxMass;
  
  Double_t bkgScalingFactor = 1.;
  switch (fMassFitTypeBkg) {
  case kExpo:
    {
      // Expo = a * exp(bx) -> integral = a/b*(TMath::Exp(b*x2) - TMath::Exp(b*x1))
      bkgScalingFactor = ((TMath::Exp(fFunction->GetParameter(1)*maxMass)  - TMath::Exp(fFunction->GetParameter(1)*minMass))/
                          (TMath::Exp(fFunction->GetParameter(1)*fMaxMass) - TMath::Exp(fFunction->GetParameter(1)*fMinMass)));
      break;
    }
  case kExpoPower:
    {
      // ExpoPower = a * sqrt(x - mpi) * exp(-b*(x-mpi)) -> integral = a / sqrt(b^3*) * TMath::Gamma(3/2) * (TMath::Gamma(3/2, b*(x2-mpi)) - TMath::Gamma(3/2, b*(x1-mpi)))
      bkgScalingFactor = ((TMath::Gamma(1.5) * (TMath::Gamma(1.5, fFunction->GetParameter(1)*(maxMass-fPionMass)) - TMath::Gamma(1.5, fFunction->GetParameter(1)*(minMass-fPionMass)))) /
                          (TMath::Gamma(1.5) * (TMath::Gamma(1.5, fFunction->GetParameter(1)*(fMaxMass-fPionMass)) - TMath::Gamma(1.5, fFunction->GetParameter(1)*(fMinMass-fPionMass)))));
      break;
    }
  }

  fBackground = (fFunction->GetParameter(0) - fSignal)*bkgScalingFactor;
  fBackgroundError = TMath::Sqrt(fFunction->GetParError(0)*fFunction->GetParError(0) + fSignalError*fSignalError)*bkgScalingFactor;
  
  return r;
}

//____________________________________________________________________________________
TString MassFitter::GetSignalString() const
{
  TString r;
  
  Double_t sigErrLog10 = TMath::Log10(fSignalError);
  if (sigErrLog10 < 0) {
    Int_t sigPrec = TMath::CeilNint(-sigErrLog10);
    r = Form("S=%%.%df #pm %%.%df", sigPrec, sigPrec);
    r = Form(r.Data(), fSignal, fSignalError);
  }
  else {
    r = Form("S=%.0f #pm %.0f", fSignal, fSignalError);
  }

  return r;
}

//____________________________________________________________________________________
TString MassFitter::GetBackgroundString() const
{
  TString r;
  
  Double_t bkgErrLog10 = TMath::Log10(fBackgroundError);
  if (bkgErrLog10 < 0) {
    Int_t bkgPrec = TMath::CeilNint(-bkgErrLog10);
    r = Form("B=%%.%df #pm %%.%df", bkgPrec, bkgPrec);
    r = Form(r.Data(), fBackground, fBackgroundError);
  }
  else {
    r = Form("B=%.0f #pm %.0f", fBackground, fBackgroundError);
  }

  return r;
}

//____________________________________________________________________________________
TString MassFitter::GetSignalOverBackgroundString() const
{
  TString r;
  
  Double_t v = GetSignalOverBackground();
  Double_t vlog10 = TMath::Log10(v);
  if (vlog10 > -1) vlog10 = -1;
  
  Int_t prec = TMath::CeilNint(-vlog10);
  r = Form("S/B=%%.%df", prec);
  r = Form(r.Data(), v);

  return r;
}

//____________________________________________________________________________________
TString MassFitter::GetSignalOverSqrtSignalBackgroundString() const
{
  TString r;
  
  Double_t v = GetSignalOverSqrtSignalBackgorund();
  Double_t vlog10 = TMath::Log10(v);
  if (vlog10 > -1) vlog10 = -1;

  Int_t prec = TMath::CeilNint(-vlog10);
  r = Form("S/#sqrt{S+B}=%%.%df", prec);
  r = Form(r.Data(), v);

  return r;
}

//____________________________________________________________________________________
TString MassFitter::GetSignalMeanString() const
{
  TString r;

  Double_t v = GetSignalMean();
  Double_t vErr = GetSignalMeanError();
  Double_t vLog10 = TMath::Log10(v);
  Double_t vErrLog10 = TMath::Log10(vErr);

  TString unit("GeV/#it{c}^{2}");
  
  if (vLog10 <-4) {
    v *= 1000000.;
    vErr *= 1000000.;
    vLog10 += 3;
    vErrLog10 += 3;
    unit = "keV/#it{c}^{2}";
  }
  else if (vLog10 <-2) {
    v *= 1000.;
    vErr *= 1000.;
    vLog10 += 3;
    vErrLog10 += 3;
    unit = "MeV/#it{c}^{2}";
  }
  
  if (vErrLog10 < 0) {
    Int_t vPrec = TMath::CeilNint(-vErrLog10);
    r = Form("#mu=%%.%df #pm %%.%df GeV/#it{c}^{2}", vPrec, vPrec);
    r = Form(r.Data(), v, vErr);
  }
  else {
    r = Form("#mu=%.0f #pm %.0f GeV/#it{c}^{2}", v, vErr);
  }

  return r;
}

//____________________________________________________________________________________
TString MassFitter::GetSignalWidthString() const
{
  TString r;

  Double_t v = GetSignalWidth();
  Double_t vErr = GetSignalWidthError();
  Double_t vLog10 = TMath::Log10(v);
  Double_t vErrLog10 = TMath::Log10(vErr);

  TString unit("GeV/#it{c}^{2}");
  
  if (vLog10 <-4) {
    v *= 1000000.;
    vErr *= 1000000.;
    vLog10 += 3;
    vErrLog10 += 3;
    unit = "keV/#it{c}^{2}";
  }
  else if (vLog10 <-2) {
    v *= 1000.;
    vErr *= 1000.;
    vLog10 += 3;
    vErrLog10 += 3;
    unit = "MeV/#it{c}^{2}";
  }
  
  if (vErrLog10 < 0) {
    Int_t vPrec = TMath::CeilNint(-vErrLog10);
    r = Form("#sigma=%%.%df #pm %%.%df %s", vPrec, vPrec, unit.Data());
    r = Form(r.Data(), v, vErr);
  }
  else {
    r = Form("#sigma=%.0f #pm %.0f %s", v, vErr, unit.Data());
  }

  return r;
}

//____________________________________________________________________________________
TString MassFitter::GetBkgPar1String() const
{
  TString r;

  if (fNParBkg > 1) {
  
    Double_t v = fFunction->GetParameter(1);
    Double_t vErr = fFunction->GetParError(1);
    Double_t vErrLog10 = TMath::Log10(vErr);
    if (vErrLog10 < 0) {
      Int_t vPrec = TMath::CeilNint(-vErrLog10);
      r = Form("b=%%.%df #pm %%.%df", vPrec, vPrec);
      r = Form(r.Data(), v, vErr);
    }
    else {
      r = Form("b=%.0f #pm %.0f", v, vErr);
    }
  }

  return r;
}

//____________________________________________________________________________________
TString MassFitter::GetBkgPar2String() const
{
  TString r;

  if (fNParBkg > 2) {
  
    Double_t v = fFunction->GetParameter(2);
    Double_t vErr = fFunction->GetParError(2);
    Double_t vErrLog10 = TMath::Log10(vErr);
    if (vErrLog10 < 0) {
      Int_t vPrec = TMath::CeilNint(-vErrLog10);
      r = Form("c=%%.%df #pm %%.%df", vPrec, vPrec);
      r = Form(r.Data(), v, vErr);
    }
    else {
      r = Form("c=%.0f #pm %.0f", v, vErr);
    }
  }

  return r;
}

//____________________________________________________________________________________
TString MassFitter::GetChisquareString() const
{
  TString r;
  
  Double_t v = fFunction->GetChisquare() / fFunction->GetNDF();
  r = Form("#chi^{2}/NdF=%.2f", v);

  return r;
}

//____________________________________________________________________________________
double MassFitter::FunctionSig(double *x, double *p)
{
  Double_t r = 0.;
  
  switch (fMassFitTypeSig) {
  case kGaus:
    {
      // Gaus = [0] * exp(-1/2*((x-[1])^2/[2]^2))  -> integral is sqrt(2pi)*[2]
      if (p[2] != 0.) {
        r = p[0] * TMath::Exp(-0.5 * (x[0] - p[1])*(x[0] - p[1]) / (p[2]*p[2]));
        r *= 1. / (TMath::Sqrt(TMath::TwoPi()) * p[2]); // p[0] is now the integral between -inf and +inf
      }
      break;
    }
  default:
    {
      Printf("Error: signal fit type %d not recognized! Using Gaussian fit.", fMassFitTypeSig);
      fMassFitTypeSig = kGaus;
      fNParSig = 3;
      r = FunctionSig(x, p);
    }
  }

  return r;
}

//____________________________________________________________________________________
double MassFitter::FunctionBkg(double *x, double *p)
{
  Double_t r = 0.;
  
  switch (fMassFitTypeBkg) {
  case kExpo:
    {
      // Expo = a * exp(bx) -> integral = a/b*(TMath::Exp(b*x2) - TMath::Exp(b*x1))
      r = p[0] * TMath::Exp(p[1]*x[0]);
      r *= p[1] / (TMath::Exp(p[1]*fMaxMass) - TMath::Exp(p[1]*fMinMass)); // p[0] is now the integral between fMinMass and fMaxMass
      break;
    }
  case kExpoPower:
    {
      // ExpoPower = a * sqrt(x - mpi) * exp(-b*(x-mpi)) -> integral = a / sqrt(b^3*) * TMath::Gamma(3/2) * (TMath::Gamma(3/2, b*(x2-mpi)) - TMath::Gamma(3/2, b*(x1-mpi)))
      if (x[0] < fPionMass) {
        r = 0;
        break;
      }
      r = p[0] * TMath::Sqrt(x[0] - fPionMass) * TMath::Exp(-p[1] * (x[0] - fPionMass));
      r *= TMath::Sqrt(p[1]*p[1]*p[1]) / (TMath::Gamma(1.5) * (TMath::Gamma(1.5, p[1]*(fMaxMass-fPionMass)) - TMath::Gamma(1.5, p[1]*(fMinMass-fPionMass)))); // p[0] is now the integral between fMinMass and fMaxMass
      break;
    }
  default:
    {
      Printf("Error: background fit type %d not recognized! Using exponential fit.", fMassFitTypeBkg);
      fMassFitTypeBkg = kExpo;
      fNParBkg = 2;
      r = FunctionBkg(x, p);
    }
  }

  return r;
}

//____________________________________________________________________________________
double MassFitter::FunctionSigBkg(double *x, double *p)
{
  Double_t r = 0;

  // p[0] = Total integral
  // p[1..fNParBkg-1] = Bkg pars
  // p[fNParBkg] = Sig integral
  // p[fNParBkg+1..fNParBkg+fNParSig] = Sig pars
  
  if (!fDisableBkg) {
    Double_t bkgPars[10] = {0.};
    for (Int_t i = 0; i < fNParBkg; i++) {
      bkgPars[i] = p[i];
    }
    bkgPars[0] -= p[fNParBkg];
    
    r += FunctionBkg(x, bkgPars);
  }
  
  if (!fDisableSig) {
    Double_t sigPars[10] = {0.};
    for (Int_t i = 0; i < fNParSig; i++) {
      sigPars[i] = p[i+fNParBkg];
    }
    
    r += FunctionSig(x, sigPars);
  }

  return r;
}

//____________________________________________________________________________________
void MassFitter::Draw(Option_t* opt)
{
  TString optAndSame(opt);
  if (!optAndSame.Contains("same")) optAndSame += " same";

  for (Int_t i = 1; i < fNParBkg; i++) {
    fFunctionBkg->SetParameter(i, fFunction->GetParameter(i));
  }
  fFunctionBkg->SetParameter(0, fFunction->GetParameter(0) - fFunction->GetParameter(fNParBkg));
  
  fFunctionBkg->Draw(opt);
  fFunction->Draw(optAndSame);
}

//____________________________________________________________________________________
void MassFitter::SetMassRange(Double_t min, Double_t max)
{
  if (fMinMass < fMaxMass) {
    fMinMass = min;
    fMaxMass = max;
  }
  else {
    Printf("Error: min mass %.3f must be smaller then mass max %.3f!", min, max);
  }
}
