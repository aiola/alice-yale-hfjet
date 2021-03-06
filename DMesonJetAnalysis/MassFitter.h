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
  enum EMassFitTypeBkg { kExpo, kExpoPower };
  // Expo = ab * exp(bx) -> integral = a*(TMath::Exp(b*x2) - TMath::Exp(b*x1))
  // ExpoPower = a * sqrt(b^3*(x - mpi)) * exp(-b*(x-mpi)) -> integral = a * TMath::Gamma(3/2) * (TMath::Gamma(3/2, b*(x2-mpi)) - TMath::Gamma(3/2, b*(x1-mpi)))

  enum EMeson { kDzeroKpi, kDstarKpipi };

  static const Double_t fgkEpsilon;
  
  MassFitter();
  MassFitter(const char* name, EMeson m, Double_t minMass, Double_t maxMass);

  ~MassFitter();

  void SetHistogram(TH1* histo);
  void SetMassFitTypeSig(EMassFitTypeSig ts)                       { fMassFitTypeSig = ts                        ; Reset(); }
  void SetMassFitTypeBkg(EMassFitTypeBkg tb)                       { fMassFitTypeBkg = tb                        ; Reset(); }
  void SetMassFitTypes(EMassFitTypeSig ts, EMassFitTypeBkg tb)     { fMassFitTypeSig = ts; fMassFitTypeBkg = tb  ; Reset(); }
  void SetReflectionTemplate(TH1* refl, Double_t ros);
  
  void Reset(TH1* histo=0);
  TFitResultPtr Fit(TH1* histo, Option_t* opt="");
  TFitResultPtr Fit(Option_t* opt="");

  void DisableBkg(Bool_t d = kTRUE);
  void DisableRefl(Bool_t d = kTRUE);
  void DisableSig(Bool_t d = kTRUE);
  void SetMassRange(Double_t min, Double_t max);
  void SetFitRange(Double_t min, Double_t max);

  Double_t GetSignalMean()                           const { return fDisableSig == kFALSE ? fMean : 0.            ; }
  Double_t GetSignalMeanError()                      const { return fDisableSig == kFALSE ? fMeanError : 0.       ; }
  Double_t GetSignalWidth()                          const { return fDisableSig == kFALSE ? fWidth : 0.           ; }
  Double_t GetSignalWidthError()                     const { return fDisableSig == kFALSE ? fWidthError : 0.      ; }
  Double_t GetSignal()                               const { return fDisableSig == kFALSE ? fSignal : 0.          ; }
  Double_t GetSignalError()                          const { return fDisableSig == kFALSE ? fSignalError : 0.     ; }
  Double_t GetSignal(Double_t sigmas)                const;
  Double_t GetSignalError(Double_t sigmas)           const;
  Double_t GetBackground(Double_t sigmas)            const;
  Double_t GetBackgroundError(Double_t sigmas)       const;
  Double_t GetBackgroundAndError(Double_t& bkgErr,
      Double_t sigmas)                               const;
  Double_t GetBackgroundAndError(Double_t& bkgErr,
      Double_t minMass, Double_t maxMass)            const;
  Double_t GetBackgroundAndErrorFullRange(Double_t& bkgErr) const;
  Double_t GetSignalOverBackground()                 const;
  Double_t GetSignificance()                         const;
  Double_t GetSignificanceW()                        const;
  Double_t GetChisquare()                            const;
  Double_t GetChisquareW()                           const;
  Double_t GetReflOverSign()                         const { return fReflOverSign; }

  Double_t GetTotalEntriesAndError(Double_t& err)    const;

  TString  GetSignalString()                         const;
  TString  GetBackgroundString()                     const;
  TString  GetSignalOverBackgroundString()           const;
  TString  GetSignificanceString()                   const;
  TString  GetSignificanceWString()                  const;
  TString  GetChisquareString()                      const;
  TString  GetChisquareWString()                     const;
  TString  GetSignalMeanString()                     const;
  TString  GetSignalWidthString()                    const;
  TString  GetBkgPar1String()                        const;
  TString  GetBkgPar2String()                        const;
  TString  GetTotalEntriesString()                   const;
  TString  GetReflOverSignString()                   const;
  
  TF1*     GetFitFunction()                          const { return fFunction           ; }
  TF1*     GetBkgFunction()                          const { return fFunctionBkg        ; }
  TF1*     GetBkgFunctionNoRefl()                    const { return fFunctionBkgNoRefl  ; }
  TH1*     GetFitHistogram()                         const { return fHistogram          ; }

  TFitResultPtr GetFitStatus()                       const { return fFitResult          ; }
  Bool_t        FitSuccessfull()                     const { return fFitSuccessfull     ; }

  void     Draw(Option_t* opt = "");
  
  double   FunctionSig(double *x, double *p);
  double   FunctionBkg(double *x, double *p);
  double   FunctionBkgNoRefl(double *x, double *p);
  double   FunctionSigBkg(double *x, double *p);
  double   FunctionRefl(double *x, double*);

  static TString GetValueString(Double_t value, Double_t err);

 protected:
  EMassFitTypeSig   fMassFitTypeSig    ;//  Mass fit type for the signal
  EMassFitTypeBkg   fMassFitTypeBkg    ;//  Mass fit type for the background
  TH1*              fReflectionTempl   ;//  MC reflection template for the D0
  Double_t          fReflOverSign      ;//  Reflection over signal ratio
  Double_t          fMean              ;//  Signal mean extracted from the fit
  Double_t          fMeanError         ;//  Signal mean error extracted from the fit
  Double_t          fWidth             ;//  Signal width extracted from the fit
  Double_t          fWidthError        ;//  Signal width error extracted from the fit
  Double_t          fSignal            ;//  Signal extracted from the fit
  Double_t          fSignalError       ;//  Signal error extracted from the fit
  Double_t          fBackground        ;//  Background extracted from the fit in the full mass range
  Double_t          fBackgroundError   ;//  Background error extracted from the fit in the full mass range
  Bool_t            fDisableBkg        ;//  Disable background component of the fit function
  Bool_t            fDisableRefl       ;//  Disable reflection component of the fit function
  Bool_t            fDisableSig        ;//  Disable signal component of the fit function
  Int_t             fNParSig           ;//  Number of parameters of the signal function
  Int_t             fNParBkg           ;//  Number of parameters of the background function
  Double_t          fMinMass           ;//  Minimum mass value
  Double_t          fMaxMass           ;//  Maximum mass value
  Double_t          fMinFitRange       ;//  Minimum fit range
  Double_t          fMaxFitRange       ;//  Maximum fit range
  TFitResultPtr     fFitResult         ;//  Fit result
  
  TF1*              fFunction          ;//  Fit function
  TF1*              fFunctionBkgNoRefl ;//  Bkg function (w/o reflections)
  TF1*              fFunctionBkg       ;//  Bkg function
  TH1*              fHistogram         ;//  Histogram to be fitted
  Bool_t            fFitSuccessfull    ;//  Whether the fit was successful
  Double_t          fPDGMass           ;//  PDG mass
  Double_t          fMaxAllowedWidth   ;//  Maximum allowed width of the signal peak
  Double_t          fMaxAllowedMeanShift;// Maximum allowed mean shift of the signal peak
  Double_t          fDefaultNSigmas    ;//  Default number of sigmas (for S/B, etc.)
  const Double_t    fPionMass          ;//! Pion mass 
  
 private:
   
  MassFitter(const MassFitter &source);
  MassFitter& operator=(const MassFitter& source);

  ClassDef(MassFitter, 1);
};

#endif
