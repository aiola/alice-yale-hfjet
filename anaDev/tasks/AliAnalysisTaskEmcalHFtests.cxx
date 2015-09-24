//
// Emcal HF tests analysis task.
//
// Author: S.Aiola

#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TH1F.h>
#include <TList.h>
#include <TMath.h>

#include "AliLog.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliParticleContainer.h"

#include "AliAnalysisTaskEmcalHFtests.h"

ClassImp(AliAnalysisTaskEmcalHFtests)

//________________________________________________________________________
AliAnalysisTaskEmcalHFtests::AliAnalysisTaskEmcalHFtests() : 
  AliAnalysisTaskEmcal("AliAnalysisTaskEmcalSample", kTRUE),
  fGenNameMustContain(""),
  fHistPtHardAll(0),
  fHistPtHardCharm(0),
  fHistPtHardBottom(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalHFtests::AliAnalysisTaskEmcalHFtests(const char *name) : 
  AliAnalysisTaskEmcal(name, kTRUE),
  fGenNameMustContain(""),
  fHistPtHardAll(0),
  fHistPtHardCharm(0),
  fHistPtHardBottom(0)
{
  // Standard constructor.
  
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalHFtests::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  fHistPtHardAll = new TH1F("fHistPtHardAll", "fHistPtHardAll", 500, 0, 500);
  fOutput->Add(fHistPtHardAll);

  fHistPtHardCharm = new TH1F("fHistPtHardCharm", "fHistPtHardCharm", 500, 0, 500);
  fOutput->Add(fHistPtHardCharm);

  fHistPtHardBottom = new TH1F("fHistPtHardBottom", "fHistPtHardBottom", 500, 0, 500);
  fOutput->Add(fHistPtHardBottom);
  
  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHFtests::IsEventSelected()
{
  if (!fPythiaHeader) {
    Printf("Could not find PYTHIA generator event header!");
    return kTRUE;
  }

  TString pname(fPythiaHeader->GetName());

  if (!fGenNameMustContain.IsNull() && !pname.Contains(fGenNameMustContain)) return kFALSE;

  return AliAnalysisTaskEmcal::IsEventSelected();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHFtests::FillHistograms()
{
  // Fill histograms.

  if (!fPythiaHeader) {
    Printf("Could not find PYTHIA generator event header!");
    return kTRUE;
  }

  TString pname(fPythiaHeader->GetName());

  TObjArray* elements = pname.Tokenize('_');
  TObjString* pthardObjStr = static_cast<TObjString*>(elements->At(elements->GetEntriesFast()-2));
  TString pthardStr = pthardObjStr->GetString();
  Int_t pthard = pthardStr.Atoi();

  fHistPtHardAll->Fill(pthard);

  if (DoesEventContain(4)) fHistPtHardCharm->Fill(fPtHard);
  if (DoesEventContain(5)) fHistPtHardBottom->Fill(fPtHard);

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHFtests::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;  // If return kFALSE FillHistogram() will NOT be executed.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalHFtests::DoesEventContain(Int_t pdg)
{
  AliParticleContainer* partCont = GetParticleContainer(0);
  if (!partCont) return kFALSE;
  
  AliVParticle* part = 0;

  partCont->ResetCurrentID();
  while ((part = partCont->GetNextParticle())) {
    if (TMath::Abs(part->PdgCode()) == pdg) return kTRUE;
  }

  return kFALSE;
}
