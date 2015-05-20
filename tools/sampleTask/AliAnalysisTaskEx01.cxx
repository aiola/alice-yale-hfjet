/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnalysisTaskEx01.cxx
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Instructions for adding histograms can be found below, starting with NEW HISTO
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 */
#include "AliAnalysisTaskEx01.h"

#include <Riostream.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TClonesArray.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliVParticle.h"

ClassImp(AliAnalysisTaskEx01)

//________________________________________________________________________
AliAnalysisTaskEx01::AliAnalysisTaskEx01() :
AliAnalysisTaskSE(),
  fTracksName(),
  fOutput(0),
  fHistPt(0), 
  fHistEta(0)
{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskEx01::AliAnalysisTaskEx01(const char *name) :
  AliAnalysisTaskSE(name),
  fTracksName("tracks"),
  fOutput(0),
  fHistPt(0), 
  fHistEta(0)
{
  // Constructor
  // Define input and output slots here (never in the dummy constructor)

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskEx01::~AliAnalysisTaskEx01()
{
  
}

//________________________________________________________________________
void AliAnalysisTaskEx01::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)
        
  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
    
  // Create histograms
  Int_t ptbins = 15;
  Float_t ptlow = 0.1, ptup = 3.1;
  fHistPt = new TH1F("fHistPt", "P_{T} distribution for reconstructed", ptbins, ptlow, ptup);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
        
  Int_t etabins = 40;
  Float_t etalow = -2.0, etaup = 2.0;
  fHistEta = new TH1F("fHistEta","#eta distribution for reconstructed",etabins, etalow, etaup);
  fHistEta->GetXaxis()->SetTitle("#eta");
  fHistEta->GetYaxis()->SetTitle("counts");
        
  // NEW HISTO should be defined here, with a sensible name,
        
  fOutput->Add(fHistPt);
  fOutput->Add(fHistEta);
  // NEW HISTO added to fOutput here
  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskEx01::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
        
  // Create pointer to reconstructed event
  AliVEvent *event = InputEvent();
  if (!event) { Printf("ERROR: Could not retrieve event"); return; }

  TClonesArray* tracks = static_cast<TClonesArray*>(event->FindListObject(fTracksName));
  if (!tracks) { Printf("ERROR: Could not retrieve tracks"); return; }
  
  // Track loop for reconstructed event
  Int_t ntracks = tracks->GetEntriesFast();
  for(Int_t i = 0; i < ntracks; i++) {
    AliVParticle* esdtrack = static_cast<AliVParticle*>(tracks->At(i)); // pointer to reconstructed to track
                
    fHistPt->Fill(esdtrack->Pt());
    fHistEta->Fill(esdtrack->Eta());
  }
  // NEW HISTO should be filled before this point, as PostData puts the
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
}
