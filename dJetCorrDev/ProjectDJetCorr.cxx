// This class projects the THnSparse results of the AliAnalysisTaskDmesonJetCorrelations into TH1/TH2/TH3
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMap.h>
#include <Riostream.h>
#include <TParameter.h>
#include <TObjString.h>
#include <TString.h>
#include <TDirectoryFile.h>
#include <TSystem.h>

#include "ProjectDJetCorr.h"

const Double_t ProjectDJetCorr::fgkEpsilon = 1E-6;

//____________________________________________________________________________________
ProjectDJetCorr::ProjectDJetCorr() :
  fUseTestResults(kFALSE),
  fTrainName(),
  fInputPath(),
  fInputFileName(),
  fInputDirFileName(),
  fInputListName(),
  fOutputPath(),
  fOutputFileName(),
  fOverwrite(kFALSE),
  fTHnAxisMap(),
  fDPtAxisTitle(),
  fDEtaAxisTitle(),
  fDPhiAxisTitle(),
  fDInvMassAxisTitle(),
  fD2prongInvMassAxisTitle(),
  fDDeltaInvMassAxisTitle(),
  fDSoftPionPtAxisTitle(),
  fDzAxisTitle(),
  fDeltaRAxisTitle(),
  fDeltaEtaAxisTitle(),
  fDeltaPhiAxisTitle(),
  fJetPtAxisTitle(),
  fJetEtaAxisTitle(),
  fJetPhiAxisTitle(),
  fJetLeadPtAxisTitle(),
  fJetAreaAxisTitle(),
  fJetConstAxisTitle(),
  fTHnSparseMapGenerated(kFALSE),
  fInputFile(0),
  fInputDirectoryFile(0),
  fInputList(0),
  fDmesons(0),
  fOutputList(new TList())
{
  // Default ctr.

  fTHnAxisMap.SetOwnerKeyValue();
}

//____________________________________________________________________________________
ProjectDJetCorr::ProjectDJetCorr(const char* train, const char* path) :
  fUseTestResults(kFALSE),
  fTrainName(train),
  fInputPath(path),
  fInputFileName("AnalysisResults.root"),
  fInputDirFileName("SA_DmesonJetCorr"),
  fInputListName(""),
  fOutputPath("../data/"),
  fOutputFileName("<train>/DJetCorr.root"),
  fOverwrite(kFALSE),
  fTHnAxisMap(),
  fDPtAxisTitle("#it{p}_{T,D} (GeV/#it{c})"),
  fDEtaAxisTitle("#eta_{D}"),
  fDPhiAxisTitle("#phi_{D} (rad)"),
  fDInvMassAxisTitle("#it{M}_{D} (GeV/#it{c}^{2})"),
  fD2prongInvMassAxisTitle("#it{M}_{2-prong} (GeV/#it{c}^{2})"),
  fDDeltaInvMassAxisTitle("#it{M}_{D*} - #it{M}_{D_{0}} (GeV/#it{c}^{2})"),
  fDSoftPionPtAxisTitle("#it{p}_{T,#pi} (GeV/#it{c})"),
  fDzAxisTitle("#it{z}_{D}"),
  fDeltaRAxisTitle("#Delta R_{D-jet}"),
  fDeltaEtaAxisTitle("#eta_{D} - #eta_{jet}"),
  fDeltaPhiAxisTitle("#phi_{D} - #phi_{jet} (rad)"),
  fJetPtAxisTitle("#it{p}_{T,jet} (GeV/#it{c})"),
  fJetEtaAxisTitle("#eta_{jet}"),
  fJetPhiAxisTitle("#phi_{jet} (rad)"),
  fJetLeadPtAxisTitle("#it{p}_{T,particle}^{leading} (GeV/#it{c})"),
  fJetAreaAxisTitle("#it{A}_{jet}"),
  fJetConstAxisTitle("No. of constituents"),
  fTHnSparseMapGenerated(kFALSE),
  fInputFile(0),
    fInputDirectoryFile(0),
  fInputList(0),
  fDmesons(0),
  fOutputList(new TList())
{
  // Standard ctr.

  fTHnAxisMap.SetOwnerKeyValue();
}

//____________________________________________________________________________________
void ProjectDJetCorr::SetInputListName(const char* dmeson, const char* jetType, const char* jetRadius)
{
  // Set input list name using the D meson type, jet type and jet radius provided.
  
  TString lname(Form("AliAnalysisTaskDmesonJetCorrelations_%s_Jet_AKT%s%s_AODFilterTracks_pT0150_pt_scheme_TPC_histos", dmeson, jetType, jetRadius));
  SetInputListName(lname);
}

//____________________________________________________________________________________
void ProjectDJetCorr:: GenerateAxisMap(THnSparse* hn)
{
  fTHnAxisMap.Delete();

  if (!LoadTHnSparse()) {
    return;
  }
  
  for (Int_t i = 0; i < hn->GetNdimensions(); i++) {
    TObjString* key = new TObjString(hn->GetAxis(i)->GetTitle());
    TParameter<Int_t>* value = new TParameter<Int_t>(hn->GetAxis(i)->GetTitle(), i);
    fTHnAxisMap.Add(key, value);
  }

  fTHnSparseMapGenerated = kTRUE;
}

//____________________________________________________________________________________
Bool_t ProjectDJetCorr::Run()
{
  // Execute the task.
  
  Bool_t result = kFALSE;
  
  result = OpenInputFile();
  if (!result) return kFALSE;
  
  result = ProjectNonCorrD();
  if (!result) return kFALSE;

  result = ProjectCorrD();
  if (!result) return kFALSE;
  
  result = SaveOutputFile();
  if (!result) return kFALSE;

  return kTRUE;
}

//____________________________________________________________________________________
Int_t ProjectDJetCorr::GetAxisIndex(TString title)
{
  if (!fTHnSparseMapGenerated) {
    GenerateAxisMap(fDmesons);
    if (!fTHnSparseMapGenerated) return -1;
  }

  TParameter<Int_t>* par = static_cast<TParameter<Int_t>*>(fTHnAxisMap.GetValue(title));
  if (!par) {
    Printf("Error-ProjectDJetCorr::GetAxisIndex : could not find axis with title '%s'", title.Data());
    return -1;
  }
  
  return par->GetVal();
}

//____________________________________________________________________________________
Bool_t ProjectDJetCorr::ProjectNonCorrD()
{
  // Project histograms related to the D meson not correlated to a jet.
  
  return ProjectDmesonJetCorr("NonCorr", 0, 1, 0, 100, kFALSE);
}

//____________________________________________________________________________________
Bool_t ProjectDJetCorr::ProjectCorrD()
{
  // Project histograms related to the D meson correlated to a jet.
  
  return ProjectDmesonJetCorr("Corr", 1, 250, 0, 100, kTRUE);
}
  
//____________________________________________________________________________________
Bool_t ProjectDJetCorr::ProjectDmesonJetCorr(const char* name, Double_t minJetPt, Double_t maxJetPt, Double_t minDPt, Double_t maxDPt, Bool_t /*doCorrPlots*/)
{
  // Project histograms related to the D meson with specified cuts.

  if (!LoadTHnSparse()) {
    return kFALSE;
  }

  Int_t dPtAxis = GetAxisIndex(fDPtAxisTitle);
  if (dPtAxis < 0) return kFALSE;
  
  Int_t dEtaAxis = GetAxisIndex(fDEtaAxisTitle);
  if (dEtaAxis < 0) return kFALSE;
  
  Int_t dPhiAxis = GetAxisIndex(fDPhiAxisTitle);
  if (dPhiAxis < 0) return kFALSE;
  
  Int_t jetPtAxis = GetAxisIndex(fJetPtAxisTitle);
  if (jetPtAxis < 0) return kFALSE;
    
  Int_t jetEtaAxis = GetAxisIndex(fJetEtaAxisTitle);
  if (jetEtaAxis < 0) return kFALSE;
  
  Int_t jetPhiAxis = GetAxisIndex(fJetPhiAxisTitle);
  if (jetPhiAxis < 0) return kFALSE;

  fDmesons->GetAxis(jetPtAxis)->SetRangeUser(minJetPt + fgkEpsilon, maxJetPt - fgkEpsilon);
  fDmesons->GetAxis(dPtAxis)->SetRangeUser(minDPt + fgkEpsilon, maxDPt - fgkEpsilon);

  TH1* hDpt = fDmesons->Projection(dPtAxis);
  hDpt->SetName(Form("hDmesonPt_%s", name));
  fOutputList->Add(hDpt);

  TH2* hDpos = fDmesons->Projection(dPhiAxis, dEtaAxis);
  hDpos->SetName(Form("hDmesonPhiVsEta_%s", name));
  fOutputList->Add(hDpos);
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t ProjectDJetCorr::OpenInputFile()
{
  // Open the input file and set the input list.
  
  TString fname;
  
  fname = fInputPath;
  fname += fInputFileName;

  if (fInputFile && fInputFileName != fInputFile->GetName()) {
    CloseInputFile();
  }

  if (!fInputFile) {
    fInputFile = TFile::Open(fname);
  }

  if (!fInputFile || fInputFile->IsZombie()) {
    Printf("Error-ProjectDJetCorr::OpenInputFile : could not open file '%s'", fname.Data()); 
    fInputFile = 0;
    return kFALSE;
  }

  fInputDirectoryFile = dynamic_cast<TDirectoryFile*>(fInputFile->Get(fInputDirFileName));

  if (!fInputDirectoryFile) {
    Printf("Error-ProjectDJetCorr::OpenInputFile : could not get directory '%s' from file '%s'", fInputDirFileName.Data(), fname.Data()); 
    return kFALSE;
  }
  
  fInputList = dynamic_cast<TList*>(fInputDirectoryFile->Get(fInputListName));
  
  if (!fInputList) {
    Printf("Error-ProjectDJetCorr::OpenInputFile : could not get list '%s' from directory '%s' of file '%s'", fInputListName.Data(), fInputDirFileName.Data(), fname.Data()); 
    return kFALSE;
  }
    
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t ProjectDJetCorr::LoadTHnSparse()
{
  if (!fDmesons) {
    fDmesons = static_cast<THnSparse*>(fInputList->FindObject("fDmesons"));
  }
  if (!fDmesons) {
    Printf("Error-ProjectDJetCorr::LoadTHnSparse : could not open find THnSparse 'fDmesons'"); 
    return kFALSE;
  }

  fTHnSparseMapGenerated = kFALSE;

  return kTRUE;
}

//____________________________________________________________________________________
void ProjectDJetCorr::CloseInputFile()
{
  // Close the input file.
  
  if (fInputFile) {
    fInputFile->Close();
    delete fInputFile;
    fInputFile = 0;
  }

  if (fInputList) {
    delete fInputList;
    fInputList = 0;
  }
}

//____________________________________________________________________________________
Bool_t ProjectDJetCorr::SaveOutputFile()
{
  // Save the content of fOutputList in a file.
  
  TString opt("create");
  if (fOverwrite) opt = "recreate";
  
  TString fname(fOutputPath);
  fname += fOutputFileName;
  
  if (fUseTestResults) {
    fname.ReplaceAll("<train>", "test");
  }
  else {
    fname.ReplaceAll("<train>", fTrainName);
  }

  TString path(fname);
  path.Remove(path.Last('/'));
  if (gSystem->AccessPathName(path)) {
    Printf("Info-ProjectDJetCorr::SaveOutputFile : creating directory '%s'.", path.Data()); 
    gSystem->mkdir(path.Data(), kTRUE);
  }
  
  TFile* outputFile = TFile::Open(fname, opt);

  if (!outputFile || outputFile->IsZombie()) {
    Printf("Error-ProjectDJetCorr::SaveOutputFile : could not open file '%s' to write (if file exists, you may need to set the overwrite option).", fname.Data()); 
    outputFile = 0;
    return kFALSE;
  }

  outputFile->cd();

  fOutputList->Write();

  outputFile->Close();
  delete outputFile;
  outputFile = 0;
  
  return kTRUE;
}
