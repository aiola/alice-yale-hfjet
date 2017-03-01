/**
 * \file AliTHnDJetRawYieldUncertainty.h
 * \brief Declaration of class AliTHnDJetRawYieldUncertainty
 *
 * In this header file the class AliTHnDJetRawYieldUncertainty is declared.
 * Class to extract jet Pt spectrum yield uncertainty via multi-trial approach,
 * specialized for the case of a THnBase input.
 *
 * \author Fabio Colamaria <fabio.colamaria@cern.ch>, INFN Bari
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Feb 28, 2016
 */

#ifndef ALITHNDJETRAWYIELDUNCERTAINTY_H
#define ALITHNDJETRAWYIELDUNCERTAINTY_H

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliDJetRawYieldUncertainty.h"

#include <TString.h>
class TFile;

/**
 * \class AliTHnDJetRawYieldUncertainty
 * \brief Implementation of a class to extract jet Pt spectrum yield uncertainty via multi-trial approach.
 *
 * Implementation of a class to extract jet Pt spectrum yield uncertainty via multi-trial approach,
 * specialized for the case of a THnBase input.
 */
class AliTHnDJetRawYieldUncertainty : public AliDJetRawYieldUncertainty {
public:

  AliTHnDJetRawYieldUncertainty();
  AliTHnDJetRawYieldUncertainty(const AliTHnDJetRawYieldUncertainty &source);
  virtual ~AliTHnDJetRawYieldUncertainty();

  void SetInputFilename(TString filename)   {fFileNameInput  = filename ; }
  void SetInputDirname(TString dirname)     {fDirName        = dirname  ; }
  void SetInputListname(TString listname)   {fListName       = listname ; }
  void SetInputObjectname(TString objname)  {fObjectName     = objname  ; }

  Bool_t ExtractInputMassPlotEffScale();
  Bool_t ExtractInputMassPlotSideband();

private:

  TString     fFileNameInput ; ///< Name of input file
  TString     fDirName       ; ///< Name of input directory in the root file
  TString     fListName      ; ///< Name of input list
  TString     fObjectName    ; ///< Name of input container to extract the mass plot
  TFile      *fFileInput     ; //!<!File containing the task output

  ClassDef(AliTHnDJetRawYieldUncertainty,1);
};

#endif
