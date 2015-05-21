// This class provides tools to fit invariant mass distributions
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TNamed.h>
#include <TMath.h>

class TF1;
class TH1;

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

  void SetHistogram(TH1* histo)                                    { fHistogram = histo                                 ; }
  void SetMassFitTypeSig(EMassFitTypeSig ts)                       { fMassFitTypeSig = ts                      ; Reset(); }
  void SetMassFitTypeBkg(EMassFitTypeBkg tb)                       { fMassFitTypeBkg = tb                      ; Reset(); }
  void SetMassFitTypes(EMassFitTypeSig ts, EMassFitTypeBkg tb)     { fMassFitTypeSig = ts; fMassFitTypeBkg = tb; Reset(); }
  
  void Reset(TH1* histo=0);
  TF1* Fit(TH1* histo);
  TF1* Fit();

  void DisableBkg(Bool_t d = kTRUE) { fDisableBkg = d; }
  void DisableSig(Bool_t d = kTRUE) { fDisableSig = d; }
  void SetMassRange(Double_t min, Double_t max) { fMinMass = min; fMaxMass = max; }

  Double_t GetMean()                             const { return fMean                 ; }
  Double_t GetSigma()                            const { return fSigma                ; }
  Double_t GetMeanError()                        const { return fMeanError            ; }
  Double_t GetSigmaError()                       const { return fSigmaError           ; }
  Double_t GetSignal()                           const { return fSignal               ; }
  Double_t GetBackground()                       const { return fBackground           ; }
  Double_t GetSignalError()                      const { return fSignalError          ; }
  Double_t GetBackgroundError()                  const { return fBackgroundError      ; }
  Double_t GetSignalOverBackground()             const { return fBackground > fgkEpsilon ? fSignal/fBackground : 0                              ; }
  Double_t GetSignalOverSqrtSignalBackgorund()   const { return fBackground+fSignal > fgkEpsilon ? fSignal/TMath::Sqrt(fSignal+fBackground) : 0 ; }

  TString  GetSignalString()                     const;
  TString  GetBackgorundString()                 const;
  TString  GetSignalOverBackgroundString()       const;
  TString  GetSignalOverSqrtSignalBackgorundString()   const;

  TF1*     GetFitFunction()                      const { return fFunction     ; }
  TF1*     GetBkgFunction()                      const { return fFunctionBkg  ; }
  TH1*     GetFitHistogram()                     const { return fHistogram    ; }

  void     Draw(Option_t* opt = "");
  
  double   FunctionSig(double *x, double *p);
  double   FunctionBkg(double *x, double *p);
  double   FunctionSigBkg(double *x, double *p);

 protected:
  EMassFitTypeSig   fMassFitTypeSig    ;//  Mass fit type for the signal
  EMassFitTypeBkg   fMassFitTypeBkg    ;//  Mass fit type for the background
  Double_t          fMean              ;//  Signal extracted from the fit
  Double_t          fMeanError         ;//  Signal error extracted from the fit
  Double_t          fSigma             ;//  Background extracted from the fit
  Double_t          fSigmaError        ;//  Background error extracted from the fit
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
  
  TF1*              fFunction          ;//! Fit function
  TF1*              fFunctionBkg       ;//! Bkg function
  TH1*              fHistogram         ;//! Histogram to be fitted
  const Double_t    fPionMass          ;//! Pion mass 
  
 private:
   
  MassFitter(const MassFitter &source);
  MassFitter& operator=(const MassFitter& source);

  ClassDef(MassFitter, 1);
};
