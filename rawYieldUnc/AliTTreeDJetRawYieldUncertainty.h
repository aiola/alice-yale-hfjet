/**
 * \file AliTTreeDJetRawYieldUncertainty.h
 * \brief Declaration of class AliTTreeDJetRawYieldUncertainty
 *
 * In this header file the class AliTTreeDJetRawYieldUncertainty is declared.
 * Class to extract jet Pt spectrum yield uncertainty via multi-trial approach,
 * specialized for the case of a TTree input.
 *
 * \author Fabio Colamaria <fabio.colamaria@cern.ch>, INFN Bari
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Feb 28, 2016
 */

#ifndef ALITTREEDJETRAWYIELDUNCERTAINTY_H
#define ALITTREEDJETRAWYIELDUNCERTAINTY_H

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <vector>
#include <string>
#include <TString.h>

#include "AliDJetRawYieldUncertainty.h"

/**
 * \class AliTTreeDJetRawYieldUncertainty
 * \brief Implementation of a class to extract jet Pt spectrum yield uncertainty via multi-trial approach.
 *
 * Implementation of a class to extract jet Pt spectrum yield uncertainty via multi-trial approach,
 * specialized for the case of a TTree input.
 */
class AliTTreeDJetRawYieldUncertainty : public AliDJetRawYieldUncertainty {

public:
  AliTTreeDJetRawYieldUncertainty();
  AliTTreeDJetRawYieldUncertainty(const AliTTreeDJetRawYieldUncertainty &source);
  virtual ~AliTTreeDJetRawYieldUncertainty();

  void SetInputTreename(TString treename)           { fTreeName       = treename ; }
  void SetInputDBranchname(TString dname)           { fDBranchName    = dname    ; }
  void SetInputJetBranchname(TString jetname)       { fJetBranchName  = jetname  ; }
  void AddInputFileName(std::string filename)       { fInputFileNames.push_back(filename);}
  void SetMassEdgesAndBinWidthForMassPlot(Double_t mmin, Double_t mmax, Double_t mbinwidth) { fmassmin = mmin; fmassmax = mmax; fmasswidth = mbinwidth; }

  TChain* GenerateChain();

  Bool_t ExtractInputMassPlotEffScale();
  Bool_t ExtractInputMassPlotSideband();

protected:

  std::vector<std::string> fInputFileNames ; ///< Name of input file
  TString                  fTreeName       ; ///< Name of input TTree
  TString                  fDBranchName    ; ///< Name of input branch for D meson
  TString                  fJetBranchName  ; ///< Name of input branch for jet
  Double_t                 fmassmin        ; ///< Mass lower edge of inv.mass plots
  Double_t                 fmassmax        ; ///< Mass upper edge of inv.mass plots
  Double_t                 fmasswidth      ; ///< Mass plots bin width

private:
  ClassDef(AliTTreeDJetRawYieldUncertainty,1);
};

#endif
