#include <TTree.h>

#include "DJetTreeReaderRoot6.h"

template<class S, class D>
DJetTreeReaderRoot6<S, D>::DJetTreeReaderRoot6(TTree* tree) : fDMeson(), fJets(), fTree(tree), fReader(tree), fDMesonReader(fReader, "DmesonJet"), fJetReaders() {
}

template<class S, class D>
void DJetTreeReaderRoot6<S, D>::AddJetDef(const char* jetName) {
  auto it = fJets.insert(std::pair<std::string, DJetObjectInfo>(jetName, DJetObjectInfo()));
  fJetReaders.insert(std::pair<std::string, std::pair<TTreeReaderValue<J>*, DJetObjectInfo* > >(jetName,
      std::pair<TTreeReaderValue<J>*, DJetObjectInfo* >(new TTreeReaderValue<J>(fReader, jetName),
          &(it.first->second))));
}

template<class S, class D>
bool DJetTreeReaderRoot6<S, D>::Next() {
  bool success = fReader.Next();
  if (!success) return false;
  FillDJetObjectInfo(fDMesonReader.Get(), &fDMeson);
  for (auto jetReaderIt : fJetReaders) {
    FillDJetObjectInfo(jetReaderIt.second.first->Get(), jetReaderIt.second.second);
  }
  return true;
}

template<class S, class D>
void DJetTreeReaderRoot6<S, D>::Restart() {
  fReader.Restart();
}

template<class S, class D>
void FillDJetObjectInfo(const S* source, D* dest) {
  dest->fPt = source->fPt;
  dest->fEta = source->fEta;
  dest->fPhi = source->fPhi;
}

template<>
void FillDJetObjectInfo(const AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary* source, DJetObjectInfo* dest) {
  dest->fPt = source->fPt;
  dest->fEta = source->fEta;
  dest->fPhi = source->fPhi;
  dest->fFirstPartonType = source->fFirstPartonType;
  dest->fLastPartonType = source->fLastPartonType;
  dest->fFirstPartonPt = source->fFirstPartonPt;
  dest->fLastPartonPt = source->fLastPartonPt;
}

template<>
void FillDJetObjectInfo(const AliAnalysisTaskDmesonJets::AliD0InfoSummary* source, DJetObjectInfo* dest) {
  dest->fPt = source->fPt;
  dest->fEta = source->fEta;
  dest->fPhi = source->fPhi;
  dest->fInvMass = source->fInvMass;
  dest->fSelectionType = source->fSelectionType;
}

template<>
void FillDJetObjectInfo(const AliAnalysisTaskDmesonJets::AliDStarInfoSummary* source, DJetObjectInfo* dest) {
  dest->fPt = source->fPt;
  dest->fEta = source->fEta;
  dest->fPhi = source->fPhi;
  dest->fDeltaInvMass = source->fDeltaInvMass;
  dest->f2ProngInvMass = source->f2ProngInvMass;
}

void FillDJetObjectInfo(const J* source, DJetObjectInfo* dest) {
  dest->fPt = source->fPt;
  dest->fEta = source->fEta;
  dest->fPhi = source->fPhi;
  dest->fR = source->fR;
  dest->fZ = source->fZ;
  dest->fN = source->fN;
}

template<>
void FillDJetObjectInfo(const AliAnalysisTaskDmesonJetsDetectorResponse::AliD0MatchInfoSummary* source,
    DJetObjectInfoResponse* dest) {
  FillDJetObjectInfo(&(source->fReconstructed), &(dest->fReconstructed));
  FillDJetObjectInfo(&(source->fGenerated), &(dest->fGenerated));
}

template<>
void FillDJetObjectInfo(const AliAnalysisTaskDmesonJetsDetectorResponse::AliDStarMatchInfoSummary* source,
    DJetObjectInfoResponse* dest) {
  FillDJetObjectInfo(&(source->fReconstructed), &(dest->fReconstructed));
  FillDJetObjectInfo(&(source->fGenerated), &(dest->fGenerated));
}
