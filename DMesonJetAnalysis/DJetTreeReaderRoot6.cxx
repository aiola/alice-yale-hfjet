#include <TTree.h>

#include "DJetTreeReaderRoot6.h"

template<class D>
DJetTreeReaderRoot6<D>::DJetTreeReaderRoot6(TTree* tree) : fDMeson(), fJets(), fTree(tree), fReader(tree), fDMesonReader(fReader, "DmesonJet"), fJetReaders() {
}

template<class D>
void DJetTreeReaderRoot6<D>::AddJetDef(const char* jetName) {
  auto it = fJets.insert(std::pair<std::string, DJetInfo>(jetName, DJetInfo()));
  fJetReaders.insert(std::pair<std::string, std::pair<TTreeReaderValue<J>*, DJetInfo* > >(jetName, std::pair<TTreeReaderValue<J>*, DJetInfo* >(new TTreeReaderValue<J>(fReader, jetName), &(it.first->second))));
}

template<class D>
bool DJetTreeReaderRoot6<D>::Next() {
  bool success = fReader.Next();
  if (!success) return false;
  FillDMesonInfo<D>(fDMesonReader.Get(), &fDMeson);
  for (auto jetReaderIt : fJetReaders) {
    FillDJetInfo(jetReaderIt.second.first->Get(), jetReaderIt.second.second);
  }
  return true;
}

template<class D>
void DJetTreeReaderRoot6<D>::Restart() {
  fReader.Restart();
}

template<class D>
void FillDMesonInfo(const D* source, DMesonInfo* dest) {
  dest->fPt = source->fPt;
  dest->fEta = source->fEta;
  dest->fPhi = source->fPhi;
}

template<>
void FillDMesonInfo(const AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary* source, DMesonInfo* dest) {
  dest->fPt = source->fPt;
  dest->fEta = source->fEta;
  dest->fPhi = source->fPhi;
  dest->fFirstPartonType = source->fFirstPartonType;
  dest->fLastPartonType = source->fLastPartonType;
  dest->fFirstPartonPt = source->fFirstPartonPt;
  dest->fLastPartonPt = source->fLastPartonPt;
}

template<>
void FillDMesonInfo(const AliAnalysisTaskDmesonJets::AliD0InfoSummary* source, DMesonInfo* dest) {
  dest->fPt = source->fPt;
  dest->fEta = source->fEta;
  dest->fPhi = source->fPhi;
  dest->fInvMass = source->fInvMass;
  dest->fSelectionType = source->fSelectionType;
}

template<>
void FillDMesonInfo(const AliAnalysisTaskDmesonJets::AliDStarInfoSummary* source, DMesonInfo* dest) {
  dest->fPt = source->fPt;
  dest->fEta = source->fEta;
  dest->fPhi = source->fPhi;
  dest->fDeltaInvMass = source->fDeltaInvMass;
  dest->f2ProngInvMass = source->f2ProngInvMass;
}

void FillDJetInfo(const J* source, DJetInfo* dest) {
  dest->fPt = source->fPt;
  dest->fEta = source->fEta;
  dest->fPhi = source->fPhi;
  dest->fR = source->fR;
  dest->fZ = source->fZ;
  dest->fN = source->fN;
}
