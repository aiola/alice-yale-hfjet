#ifndef ALIANALYSISTASKEMCALHFtests_H
#define ALIANALYSISTASKEMCALHFtests_H

class TH1;
class TH2;
class TH3;
class AliParticleContainer;
class AliClusterContainer;

#include <AliAnalysisTaskEmcal.h>

class AliAnalysisTaskEmcalHFtests : public AliAnalysisTaskEmcal {
 public:

  AliAnalysisTaskEmcalHFtests();
  AliAnalysisTaskEmcalHFtests(const char *name);

  void                        UserCreateOutputObjects();

  void                        SetGenNameMustContain(const char* s) { fGenNameMustContain = s; }

 protected:
  Bool_t                      IsEventSelected()          ;
  Bool_t                      FillHistograms()           ;
  Bool_t                      Run()                      ;
  Bool_t                      DoesEventContain(Int_t pdg);
  
  TString                     fGenNameMustContain ;
  TH1                        *fHistPtHardAll      ; //!
  TH1                        *fHistPtHardCharm    ; //!
  TH1                        *fHistPtHardBottom   ; //!

 private:
  AliAnalysisTaskEmcalHFtests(const AliAnalysisTaskEmcalHFtests&);            // not implemented
  AliAnalysisTaskEmcalHFtests &operator=(const AliAnalysisTaskEmcalHFtests&); // not implemented

  ClassDef(AliAnalysisTaskEmcalHFtests, 1) // emcal HF tests analysis task
};
#endif
