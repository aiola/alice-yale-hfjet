// This class projects the THnSparse results of the AliAnalysisTaskDmesonJetCorrelations into TH1/TH2/TH3
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMap.h>
#include <Riostream.h>
#include <TParameter.h>
#include <TObjString.h>
#include <TString.h>
#include <TDirectoryFile.h>
#include <TSystem.h>

#include "DJetCorrAnalysis.h"

const Double_t DJetCorrAnalysis::fgkEpsilon = 1E-6;

//____________________________________________________________________________________
DJetCorrAnalysis::DJetCorrAnalysis() :
  fUseTestResults(kFALSE),
  fTrainName(),
  fInputPath(),
  fInputFileName(),
  fInputDirFileName(),
  fInputListName(),
  fOutputPath(),
  fOutputFileName(),
  fOverwrite(kFALSE),
  fJetType(),
  fJetRadius(),
  fDmesonName(),
  fQAListName(),
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
  fNDPtBins(0),
  fDPtBins(0),
  fNJetPtBins(0),
  fJetPtBins(0),
  fTHnAxisMap(),
  fTHnSparseMapGenerated(kFALSE),
  fInputFile(0),
  fInputDirectoryFile(0),
  fInputList(0),
  fInputQAList(0),
  fDmesons(0),
  fOutputList(0)
{
  // Default ctr.

  fTHnAxisMap.SetOwnerKeyValue();
}

//____________________________________________________________________________________
DJetCorrAnalysis::DJetCorrAnalysis(const char* train, const char* path) :
  fUseTestResults(kFALSE),
  fTrainName(train),
  fInputPath(path),
  fInputFileName("AnalysisResults.root"),
  fInputDirFileName("SA_DmesonJetCorr"),
  fInputListName("AliAnalysisTaskDmesonJetCorrelations_<dmeson>_Jet_AKT<jettype><jetradius>_AODFilterTracks_pT0150_pt_scheme_TPC_histos"),
  fOutputPath("../data/"),
  fOutputFileName("<train>/DJetCorr.root"),
  fOverwrite(kFALSE),
  fJetType("Charged"),
  fJetRadius("R040"),
  fDmesonName("D0"),
  fQAListName("AliAnalysisTaskSAQA_AODFilterTracks_TPC_histos"),
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
  fNDPtBins(8),
  fDPtBins(0),
  fNJetPtBins(10),
  fJetPtBins(0),
  fTHnAxisMap(),
  fTHnSparseMapGenerated(kFALSE),
  fInputFile(0),
  fInputDirectoryFile(0),
  fInputList(0),
  fInputQAList(0),
  fDmesons(0),
  fOutputList(0)
{
  // Standard ctr.

  fTHnAxisMap.SetOwnerKeyValue();

  fDPtBins = new Double_t[fNDPtBins+1];
  fDPtBins[ 0] =  4.0;
  fDPtBins[ 1] =  5.0;
  fDPtBins[ 2] =  6.0;
  fDPtBins[ 3] =  7.0;
  fDPtBins[ 4] =  8.0;
  fDPtBins[ 5] = 10.0;
  fDPtBins[ 6] = 15.0;
  fDPtBins[ 7] = 20.0;
  fDPtBins[ 8] = 30.0;
  fDPtBins[ 9] = 50.0;

  fJetPtBins = new Double_t[fNJetPtBins+1];
  fJetPtBins[ 0] =   1.0;
  fJetPtBins[ 1] =   3.0;
  fJetPtBins[ 2] =   5.0;
  fJetPtBins[ 3] =  10.0;
  fJetPtBins[ 4] =  15.0;
  fJetPtBins[ 5] =  20.0;
  fJetPtBins[ 6] =  25.0;
  fJetPtBins[ 7] =  30.0;
  fJetPtBins[ 8] =  40.0;
  fJetPtBins[ 9] =  60.0;
  fJetPtBins[10] =  80.0;
  fJetPtBins[11] = 100.0;
}

//____________________________________________________________________________________
void DJetCorrAnalysis::SetAnalysisParams(const char* dmeson, const char* jetType, const char* jetRadius)
{
  // Set input list name using the D meson type, jet type and jet radius provided.

  fJetType = jetType;
  fJetRadius = jetRadius;
  fDmesonName = dmeson;
  
  TString lname(Form("AliAnalysisTaskDmesonJetCorrelations_%s_Jet_AKT%s%s_AODFilterTracks_pT0150_pt_scheme_TPC_histos", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data()));
  SetInputListName(lname);
}

//____________________________________________________________________________________
void DJetCorrAnalysis:: GenerateAxisMap(THnSparse* hn)
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
Bool_t DJetCorrAnalysis::Init()
{
  // Init class.

  Printf("Info-DJetCorrAnalysis::Init : Now initializing.");
  
  if (fOutputList) {
    delete fOutputList;
    fOutputList = 0;
  }
  fOutputList = new TList();

  Printf("Info-DJetCorrAnalysis::Init : Initialization done.");

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrAnalysis::ClearInputData()
{
  // Clear the input data.
  
  fDmesons = 0;
  fTHnAxisMap.Delete();
  fTHnSparseMapGenerated = kFALSE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateQAHistograms()
{
  // Generate QA histograms.

  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;
  }
  
  result = LoadQAList();
  if (!result) return kFALSE;
  
  result = ProjectQA();
  if (!result) return kFALSE;

  TH1::AddDirectory(addDirStatus);

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::GenerateDJetCorrHistograms(const char* dmeson, const char* jetType, const char* jetRadius)
{
  // Generate D-jet correlation histograms.

  SetAnalysisParams(dmeson, jetType, jetRadius);
  
  Bool_t addDirStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  Bool_t result = kFALSE;

  if (!fOutputList) {
    result = Init();
    if (!result) return kFALSE;
  }
  
  result = LoadInputList();
  if (!result) return kFALSE;

  result = ProjectCorrD();
  if (!result) return kFALSE;

  TH1::AddDirectory(addDirStatus);

  return kTRUE;
}

//____________________________________________________________________________________
Int_t DJetCorrAnalysis::GetAxisIndex(TString title)
{
  if (!fTHnSparseMapGenerated) {
    GenerateAxisMap(fDmesons);
    if (!fTHnSparseMapGenerated) return -1;
  }

  TParameter<Int_t>* par = static_cast<TParameter<Int_t>*>(fTHnAxisMap.GetValue(title));
  if (!par) {
    Printf("Warning-DJetCorrAnalysis::GetAxisIndex : could not find axis with title '%s'", title.Data());
    return -1;
  }
  
  return par->GetVal();
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectQA()
{
  // Project QA histograms.

  TH1* evCount = dynamic_cast<TH1*>(fInputQAList->FindObject("fHistEventCount"));
  evCount->SetName("hEventCount");
  fOutputList->Add(evCount);

  TString labels[4] = {"Global", "Constr", "ConstrNoITS", "Hybrid"};
  
  for (Int_t i = 0; i <= 3; i++) {
    TH3* hTrPhiEtaPt = dynamic_cast<TH3*>(fInputQAList->FindObject(Form("fHistTrPhiEtaPt_0_%d",i)));
    if (!hTrPhiEtaPt) {
      Printf("Error-DJetCorrAnalysis::ProjectQA : Could not find track histogram for track type %d!", i);
      continue;
    }
    TH1* hTrEta = hTrPhiEtaPt->ProjectionX(Form("hTracks_%s_Eta", labels[i].Data()));
    TH1* hTrPhi = hTrPhiEtaPt->ProjectionY(Form("hTracks_%s_Phi", labels[i].Data()));
    TH1* hTrPt = hTrPhiEtaPt->ProjectionZ(Form("hTracks_%s_Pt", labels[i].Data()));

    fOutputList->Add(hTrEta);
    fOutputList->Add(hTrPhi);
    fOutputList->Add(hTrPt);
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectCorrD()
{
  // Project histograms related to the D meson correlated to a jet.

  TString prefix(Form("%s_%s_%s", fDmesonName.Data(), fJetType.Data(), fJetRadius.Data()));
  TString suffix("");

  ProjectDJetCorr(prefix, suffix, kTRUE, 
                  0, 1,
                  fDPtBins[0], fDPtBins[fNDPtBins+1], -0.5, 0.5);

  for (Int_t i = 0; i < fNDPtBins; i++) {
    ProjectDJetCorr(prefix, suffix, kTRUE, 
                    0, 1,
                    fDPtBins[i], fDPtBins[i+1], -0.5, 0.5);
  }
  
  ProjectDJetCorr(prefix, suffix, kTRUE, 
                  fJetPtBins[0], fJetPtBins[fNJetPtBins+1],
                  fDPtBins[0], fDPtBins[fNDPtBins+1], -0.5, 0.5);

  for (Int_t i = 0; i < fNDPtBins; i++) {
    ProjectDJetCorr(prefix, suffix, kTRUE, 
                    fJetPtBins[0], fJetPtBins[fNJetPtBins+1],
                    fDPtBins[i], fDPtBins[i+1], -0.5, 0.5);
  }

  for (Int_t i = 0; i < fNJetPtBins; i++) {
    ProjectDJetCorr(prefix, suffix, kTRUE, 
                    fJetPtBins[i], fJetPtBins[i+1],
                    fDPtBins[0], fDPtBins[fNDPtBins+1], -0.5, 0.5);
  }

  return kTRUE;
}
  
//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::ProjectDJetCorr(TString prefix, TString suffix, Bool_t doCorrPlots, 
                                        Double_t minJetPt, Double_t maxJetPt,
                                        Double_t minDPt, Double_t maxDPt, Double_t minDEta, Double_t maxDEta)
{
  // Project histograms related to the D meson with specified cuts.

  if (!LoadTHnSparse()) {
    return kFALSE;
  }

  TString cuts(Form("JetPt_%03.0f_%03.0f_DPt_%02.0f_%02.0f", minJetPt, maxJetPt, minDPt, maxDPt));
  cuts.ReplaceAll(".", "");
  if (!suffix.IsNull()) suffix.Prepend("_");
  suffix.Prepend(cuts);

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

  Int_t dInvMassAxis = GetAxisIndex(fDInvMassAxisTitle);
  Int_t d2ProngInvMassAxis = GetAxisIndex(fD2prongInvMassAxisTitle);
  Int_t dDeltaInvMassAxis = GetAxisIndex(fDDeltaInvMassAxisTitle);
  Int_t dSoftPionPtAxis = GetAxisIndex(fDSoftPionPtAxisTitle);
  
  Int_t deltaRAxis = GetAxisIndex(fDeltaRAxisTitle);
  Int_t deltaEtaAxis = GetAxisIndex(fDeltaEtaAxisTitle);
  Int_t deltaPhiAxis = GetAxisIndex(fDeltaPhiAxisTitle);
  Int_t dzAxis = GetAxisIndex(fDzAxisTitle);

  fDmesons->GetAxis(jetPtAxis)->SetRangeUser(minJetPt * (1+fgkEpsilon), maxJetPt * (1-fgkEpsilon));
  fDmesons->GetAxis(dPtAxis)->SetRangeUser(minDPt * (1+fgkEpsilon), maxDPt * (1-fgkEpsilon));
  fDmesons->GetAxis(dEtaAxis)->SetRangeUser(minDEta * (1+fgkEpsilon), maxDEta * (1-fgkEpsilon));

  if (!cuts.Contains("DPt")) {
    TH1* hDpt = fDmesons->Projection(dPtAxis);
    hDpt->SetName(Form("h%s_MesonPt_%s", prefix.Data(), suffix.Data()));
    fOutputList->Add(hDpt);
  }
  
  TH2* hDpos = fDmesons->Projection(dPhiAxis, dEtaAxis);
  hDpos->SetName(Form("h%s_MesonPhiVsEta_%s", prefix.Data(), suffix.Data()));
  fOutputList->Add(hDpos);

  if (dInvMassAxis >= 0) {
    TH1* hdinvmass = fDmesons->Projection(dInvMassAxis);
    hdinvmass->SetName(Form("h%s_InvMass_%s", prefix.Data(), suffix.Data()));
    fOutputList->Add(hdinvmass);
  }

  if (d2ProngInvMassAxis >= 0) {
    TH1* hd2pronginvmass = fDmesons->Projection(d2ProngInvMassAxis);
    hd2pronginvmass->SetName(Form("h%s_D0InvMass_%s", prefix.Data(), suffix.Data()));
    fOutputList->Add(hd2pronginvmass);
  }

  if (dDeltaInvMassAxis >= 0) {
    TH1* hddeltainvmass = fDmesons->Projection(dDeltaInvMassAxis);
    hddeltainvmass->SetName(Form("h%s_DeltaInvMass_%s", prefix.Data(), suffix.Data()));
    fOutputList->Add(hddeltainvmass);
  }

  if (dSoftPionPtAxis >= 0) {
    TH1* hdsoftpionpt = fDmesons->Projection(dSoftPionPtAxis);
    hdsoftpionpt->SetName(Form("h%s_SoftPion_%s", prefix.Data(), suffix.Data()));
    fOutputList->Add(hdsoftpionpt);
  }

  if (doCorrPlots) {
    if (deltaRAxis >= 0) {
      TH1* hDeltaR = fDmesons->Projection(deltaRAxis);
      hDeltaR->SetName(Form("h%s_DeltaR_%s", prefix.Data(), suffix.Data()));
      fOutputList->Add(hDeltaR);
    }

    if (deltaEtaAxis >= 0) {
      TH1* hDeltaEta = fDmesons->Projection(deltaEtaAxis);
      hDeltaEta->SetName(Form("h%s_DeltaEta_%s", prefix.Data(), suffix.Data()));
      fOutputList->Add(hDeltaEta);
    }

    if (deltaPhiAxis >= 0) {
      TH1* hDeltaPhi = fDmesons->Projection(deltaPhiAxis);
      hDeltaPhi->SetName(Form("h%s_DeltaPhi_%s", prefix.Data(), suffix.Data()));
      fOutputList->Add(hDeltaPhi);
    }

    if (dzAxis >= 0) {
      TH1* hdz = fDmesons->Projection(dzAxis);
      hdz->SetName(Form("h%s_MesonZ_%s", prefix.Data(), suffix.Data()));
      fOutputList->Add(hdz);
    }
  }
  
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::OpenInputFile()
{
  // Open the input file and set the input list.
  
  TString fname;
  
  fname = fInputPath;

  if (!fUseTestResults) {
    fname += "/";
    fname += fTrainName;
  }

  fname += "/";
  fname += fInputFileName;

  if (!fInputFile) {
    Printf("Info-DJetCorrAnalysis::OpenInputFile : Opening file '%s'", fname.Data()); 
    fInputFile = TFile::Open(fname);
  }
  else {
    Printf("Info-DJetCorrAnalysis::OpenInputFile : File '%s' already open", fInputFile->GetName()); 
  }

  if (!fInputFile || fInputFile->IsZombie()) {
    Printf("Error-DJetCorrAnalysis::OpenInputFile : Could not open file '%s'", fname.Data()); 
    fInputFile = 0;
    return kFALSE;
  }

  if (!fInputDirectoryFile || fInputDirFileName != fInputDirectoryFile->GetName()) {
    ClearInputData();
    
    delete fInputDirectoryFile;
    Printf("Info-DJetCorrAnalysis::OpenInputFile : Getting directory '%s' from file '%s'", fInputDirFileName.Data(), fInputFile->GetName()); 
    fInputDirectoryFile = dynamic_cast<TDirectoryFile*>(fInputFile->Get(fInputDirFileName));

    delete fInputList;
    fInputList = 0;
  }

  if (!fInputDirectoryFile) {
    Printf("Error-DJetCorrAnalysis::OpenInputFile : Could not get directory '%s' from file '%s'", fInputDirFileName.Data(), fInputFile->GetName()); 
    return kFALSE;
  }

  Printf("Info-DJetCorrAnalysis::OpenInputFile : Success.");

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::LoadInputList()
{
  // Load the input list.

  if (!OpenInputFile()) return kFALSE;
  
  if (!fInputList || fInputListName != fInputList->GetName()) {
    ClearInputData();
    
    delete fInputList;
    Printf("Info-DJetCorrAnalysis::OpenInputFile : Getting list '%s' from directory '%s' of file '%s'", fInputListName.Data(), fInputDirectoryFile->GetName(), fInputFile->GetName()); 
    fInputList = dynamic_cast<TList*>(fInputDirectoryFile->Get(fInputListName));
  }
  
  if (!fInputList) {
    Printf("Error-DJetCorrAnalysis::OpenInputFile : Could not get list '%s' from directory '%s' of file '%s'", fInputListName.Data(), fInputDirectoryFile->GetName(), fInputFile->GetName());
    return kFALSE;
  }

  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::LoadQAList()
{
  // Load the input QA list.

  if (!OpenInputFile()) return kFALSE;
  
  if (!fInputQAList || fQAListName != fInputQAList->GetName()) {
    delete fInputQAList;
    Printf("Info-DJetCorrAnalysis::OpenInputFile : Getting list '%s' from file '%s'", fQAListName.Data(), fInputFile->GetName()); 
    fInputQAList = dynamic_cast<TList*>(fInputFile->Get(fQAListName));
  }
  
  if (!fInputQAList) {
    Printf("Error-DJetCorrAnalysis::OpenInputFile : Could not get list '%s' from file '%s'", fQAListName.Data(), fInputFile->GetName());
    return kFALSE;
  }
    
  return kTRUE;
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::LoadTHnSparse()
{
  if (!fDmesons) {
    fDmesons = static_cast<THnSparse*>(fInputList->FindObject("fDmesons"));
  }
  if (!fDmesons) {
    Printf("Error-DJetCorrAnalysis::LoadTHnSparse : could not open find THnSparse 'fDmesons'"); 
    return kFALSE;
  }

  fTHnSparseMapGenerated = kFALSE;

  return kTRUE;
}

//____________________________________________________________________________________
void DJetCorrAnalysis::CloseInputFile()
{
  // Close the input file.

  ClearInputData();
  
  if (fInputFile) {
    fInputFile->Close();
    delete fInputFile;
    fInputFile = 0;
  }

  if (fInputDirectoryFile) {
    delete fInputDirectoryFile;
    fInputDirectoryFile = 0;
  }

  if (fInputList) {
    delete fInputList;
    fInputList = 0;
  }
}

//____________________________________________________________________________________
Bool_t DJetCorrAnalysis::SaveOutputFile()
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
    Printf("Info-DJetCorrAnalysis::SaveOutputFile : creating directory '%s'.", path.Data()); 
    gSystem->mkdir(path.Data(), kTRUE);
  }
  
  TFile* outputFile = TFile::Open(fname, opt);

  if (!outputFile || outputFile->IsZombie()) {
    Printf("Error-DJetCorrAnalysis::SaveOutputFile : could not open file '%s' to write (if file exists, you may need to set the overwrite option).", fname.Data()); 
    outputFile = 0;
    return kFALSE;
  }

  outputFile->cd();

  Printf("Info-DJetCorrAnalysis::SaveOutputFile : Now streming results."); 
  fOutputList->Write();

  Printf("Info-DJetCorrAnalysis::SaveOutputFile : Closing the output file."); 
  outputFile->Close();
  delete outputFile;
  outputFile = 0;
  
  return kTRUE;
}
