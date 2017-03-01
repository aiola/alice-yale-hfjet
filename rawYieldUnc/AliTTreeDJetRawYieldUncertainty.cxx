/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TChain.h>
#include <Riostream.h>

#include <AliAnalysisTaskDmesonJets.h>

#include "AliTTreeDJetRawYieldUncertainty.h"

/// \cond CLASSIMP
ClassImp(AliTTreeDJetRawYieldUncertainty);
/// \endcond

/**
 * Default constructor.
 */
AliTTreeDJetRawYieldUncertainty::AliTTreeDJetRawYieldUncertainty():
  AliDJetRawYieldUncertainty(),
  fInputFileNames(),
  fTreeName(),
  fDBranchName(),
  fJetBranchName(),
  fmassmin(0),
  fmassmax(0),
  fmasswidth(0)
{
}

/**
 * Copy constructor.
 * @param[in] source Const reference to an object to copy from
 */
AliTTreeDJetRawYieldUncertainty::AliTTreeDJetRawYieldUncertainty(const AliTTreeDJetRawYieldUncertainty &source):
  AliDJetRawYieldUncertainty(source),
  fInputFileNames(source.fInputFileNames),
  fTreeName(source.fTreeName),
  fDBranchName(source.fDBranchName),
  fJetBranchName(source.fJetBranchName),
  fmassmin(source.fmassmin),
  fmassmax(source.fmassmax),
  fmasswidth(source.fmasswidth)
{
}

/**
 * Destructor
 */
AliTTreeDJetRawYieldUncertainty::~AliTTreeDJetRawYieldUncertainty()
{
}

/**
 * Generates the chain containing the D-tagged jet candidates.
 * @return Pointer to a new chain
 */
TChain* AliTTreeDJetRawYieldUncertainty::GenerateChain()
{
  Printf("Generating the chain...");
  TChain* chain = new TChain(fTreeName);
  for (std::string fname : fInputFileNames) {
    chain->Add(fname.c_str());
  }
  return chain;
}

/**
 * Extract the input mass plots for the efficiency scaled method.
 */
Bool_t AliTTreeDJetRawYieldUncertainty::ExtractInputMassPlotEffScale()
{
  std::cout << "Extracting input mass plot: " << fpTmin << " to " << fpTmax << std::endl;

  TTree *tree = GenerateChain();
  if (!tree) {
    std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
    return kFALSE;
  }

  AliAnalysisTaskDmesonJets::AliD0InfoSummary *brD = 0;
  AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
  tree->SetBranchAddress(fDBranchName,&brD);
  tree->SetBranchAddress(fJetBranchName,&brJet);

  TString hname = TString::Format("invMass_JetPt_%.0f_%.0f", fpTmin, fpTmax);
  fMassPlot = new TH1D(hname,hname, (fmassmax-fmassmin) / fmasswidth, fmassmin, fmassmax);
  fMassPlot->Sumw2();

  for (int k = 0; k < tree->GetEntries(); k++) {
    tree->GetEntry(k);
    if (brJet->fEta < -0.5 || brJet->fEta >= 0.5) continue;
    if (brJet->fPt < fpTmin || brJet->fPt >= fpTmax) continue;
    for (int j = 0; j < fnDbins; j++) {
      if (brD->fPt < fDbinpTedges[j] || brD->fPt >= fDbinpTedges[j + 1]) continue;
      fMassPlot->Fill(brD->fInvMass, 1. / fDEffValues[j]);
    }//end of D-meson pT bin loop
  }

  return kTRUE;
}

/**
 * Extract the input mass plots for the side band method.
 */
Bool_t AliTTreeDJetRawYieldUncertainty::ExtractInputMassPlotSideband()
{
  double jetmin = 5, jetmax = 30;

  std::cout << "Extracting input mass plot: " << fpTmin << " to " << fpTmax << std::endl;

  TTree *tree = GenerateChain();
  if (!tree) {
    std::cout << "Error in setting the tree/branch names! Exiting..." << std::endl;
    return kFALSE;
  }

  AliAnalysisTaskDmesonJets::AliD0InfoSummary *brD = 0;
  AliAnalysisTaskDmesonJets::AliJetInfoSummary *brJet = 0;
  tree->SetBranchAddress(fDBranchName,&brD);
  tree->SetBranchAddress(fJetBranchName,&brJet);

  TString hname = TString::Format("invMass_DPt_%.0f_%.0f", fpTmin, fpTmax);
  fMassPlot = new TH1D(hname,hname, (fmassmax-fmassmin) / fmasswidth, fmassmin, fmassmax);
  fMassPlot->Sumw2();

  fMassVsJetPtPlot = new TH2D("hInvMassJetPt", "hInvMassJetPt", (fmassmax-fmassmin) / fmasswidth, fmassmin, fmassmax, fnJetbins, fJetbinpTedges);
  fMassVsJetPtPlot->Sumw2();

  for (int k = 0; k < tree->GetEntries(); k++) {
    tree->GetEntry(k);
    if (brJet->fEta < -0.5 || brJet->fEta >= 0.5) continue;
    if (brJet->fPt < jetmin || brJet->fPt >= jetmax) continue;
    if (brD->fPt < fpTmin || brD->fPt >= fpTmax) continue;
    fMassPlot->Fill(brD->fInvMass);
    fMassVsJetPtPlot->Fill(brD->fInvMass, brJet->fPt);
  }

  return kTRUE;
}