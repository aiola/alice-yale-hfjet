// This class provides tools to fit invariant mass distributions
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#ifndef MASSFITTER
#define MASSFITTER

#include <TNamed.h>
#include <TMath.h>
#include <TFitResultPtr.h>
#include <TH1.h>
#include <TF1.h>

class MassFitter : public TNamed {
  
 public:
  enum EMassFitTypeSig { kGaus };
  // Gaus = a / (sqrt(2pi)*c) * exp(-1/2*((x-b)^2/c^2))  -> integral is a
  enum EMassFitTypeBkg { kExpo, kExpoPower };  //
  // Expo = ab * exp(bx) -> integral = a*(TMath::Exp(b*x2) - TMath::Exp(b*x1))
  // ExpoPower = a * sqrt(b^3*(x - mpi)) * exp(-b*(x-mpi)) -> integral = a * TMath::Gamma(3/2) * (TMath::Gamma(3/2, b*(x2-mpi)) - TMath::Gamma(3/2, b*(x1-mpi)))

  static const Double_t fgkEpsilon;
  
  MassFitter();
  MassFitter(const char* name, EMassFitTypeSig ts, EMassFitTypeBkg tb, Double_t minMass, Double_t maxMass);

  ~MassFitter();

  void SetHistogram(TH1* histo);
  void SetMassFitTypeSig(EMassFitTypeSig ts)                       { fMassFitTypeSig = ts                      ; Reset(); }
  void SetMassFitTypeBkg(EMassFitTypeBkg tb)                       { fMassFitTypeBkg = tb                      ; Reset(); }
  void SetMassFitTypes(EMassFitTypeSig ts, EMassFitTypeBkg tb)     { fMassFitTypeSig = ts; fMassFitTypeBkg = tb; Reset(); }
  
  void Reset(TH1* histo=0);
  TFitResultPtr Fit(TH1* histo, Option_t* opt="");
  TFitResultPtr Fit(Option_t* opt="");

  void DisableBkg(Bool_t d = kTRUE);
  void DisableSig(Bool_t d = kTRUE);
  void SetMassRange(Double_t min, Double_t max);
  void SetFitRange(Double_t min, Double_t max);

  Double_t GetSignalMean()                           const { return fDisableSig == kFALSE ? fMean : 0.            ; }
  Double_t GetSignalMeanError()                      const { return fDisableSig == kFALSE ? fMeanError : 0.       ; }
  Double_t GetSignalWidth()                          const { return fDisableSig == kFALSE ? fWidth : 0.           ; }
  Double_t GetSignalWidthError()                     const { return fDisableSig == kFALSE ? fWidthError : 0.      ; }
  Double_t GetSignal()                               const { return fDisableSig == kFALSE ? fSignal : 0.          ; }
  Double_t GetSignalError()                          const { return fDisableSig == kFALSE ? fSignalError : 0.     ; }
  Double_t GetBackground(Double_t sigmas=3.0)        const;
  Double_t GetBackgroundError(Double_t sigmas=3.0)   const;
  Double_t GetBackgroundAndError(Double_t& bkgErr,
      Double_t sigmas=3.0)                           const;
  Double_t GetBackgroundBinCountAndError(Double_t& error,
      Double_t minNSigmas, Double_t maxNSigmas)      const;
  Double_t GetBackgroundBinCount(Double_t minNSigmas,
      Double_t maxNSigmas)                           const { Double_t error = 0.; return GetBackgroundBinCountAndError(error, minNSigmas, maxNSigmas); }
  Double_t GetSignalOverBackground()                 const;
  Double_t GetSignalOverSqrtSignalBackgorund()       const;

  Double_t GetTotalEntries()                         const;
  Double_t GetTotalEntriesError()                    const;

  TString  GetSignalString()                         const;
  TString  GetBackgroundString()                     const;
  TString  GetSignalOverBackgroundString()           const;
  TString  GetSignalOverSqrtSignalBackgroundString() const;
  TString  GetChisquareString()                      const;
  TString  GetSignalMeanString()                     const;
  TString  GetSignalWidthString()                    const;
  TString  GetBkgPar1String()                        const;
  TString  GetBkgPar2String()                        const;
  TString  GetTotalEntriesString()                   const;
  
  TF1*     GetFitFunction()                          const { return fFunction     ; }
  TF1*     GetBkgFunction()                          const { return fFunctionBkg  ; }
  TH1*     GetFitHistogram()                         const { return fHistogram    ; }

  TFitResultPtr GetFitStatus()                       const { return fFitResult    ; }

  void     DivideByBinWidth();
  void     NormalizeBackground();

  void     Draw(Option_t* opt = "");
  
  double   FunctionSig(double *x, double *p);
  double   FunctionBkg(double *x, double *p);
  double   FunctionSigBkg(double *x, double *p);

  static TString GetValueString(Double_t value, Double_t err);

 protected:
  EMassFitTypeSig   fMassFitTypeSig    ;//  Mass fit type for the signal
  EMassFitTypeBkg   fMassFitTypeBkg    ;//  Mass fit type for the background
  Double_t          fMean              ;//  Signal mean extracted from the fit
  Double_t          fMeanError         ;//  Signal mean error extracted from the fit
  Double_t          fWidth             ;//  Signal width extracted from the fit
  Double_t          fWidthError        ;//  Signal width error extracted from the fit
  Double_t          fSignal            ;//  Signal extracted from the fit
  Double_t          fSignalError       ;//  Signal error extracted from the fit
  Double_t          fBackground        ;//  Background extracted from the fit
  Double_t          fBackgroundError   ;//  Background error extracted from the fit
  Bool_t            fDisableBkg        ;//  Disable background component of the fit function
  Bool_t            fDisableSig        ;//  Disable signal component of the fit function
  Int_t             fNParSig           ;//  Number of parameters of the signal function
  Int_t             fNParBkg           ;//  Number of parameters of the background function
  Double_t          fMinMass           ;//  Minimum mass value
  Double_t          fMaxMass           ;//  Maximum mass value
  Double_t          fMinFitRange       ;//  Minimum fit range
  Double_t          fMaxFitRange       ;//  Maximum fit range
  Double_t          fScaleFactor       ;//  Scale factor
  TFitResultPtr     fFitResult         ;//  Fit result
  
  TF1*              fFunction          ;//  Fit function
  TF1*              fFunctionBkg       ;//  Bkg function
  TH1*              fHistogram         ;//  Histogram to be fitted
  const Double_t    fPionMass          ;//! Pion mass 
  
 private:
   
  MassFitter(const MassFitter &source);
  MassFitter& operator=(const MassFitter& source);

  ClassDef(MassFitter, 1);
};

#endif
