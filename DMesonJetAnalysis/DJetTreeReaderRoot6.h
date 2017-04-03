#include <map>
#include <string>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <AliAnalysisTaskDmesonJets.h>
#include <AliAnalysisTaskDmesonJetsDetectorResponse.h>

class TTree;

typedef AliAnalysisTaskDmesonJets::AliJetInfoSummary J;

struct DJetObjectInfo {
  double  fPt             ;
  double  fEta            ;
  double  fPhi            ;
  int     fFirstPartonType;
  int     fLastPartonType ;
  double  fFirstPartonPt  ;
  double  fLastPartonPt   ;
  double  fInvMass        ;
  int     fSelectionType  ;
  double  fDeltaInvMass   ;
  double  f2ProngInvMass  ;
  double  fR  ;
  double  fZ  ;
  int     fN  ;
};

struct DJetObjectInfoResponse {
  DJetObjectInfo fReconstructed;
  DJetObjectInfo fGenerated;
};

template<class S, class D=DJetObjectInfo>
class DJetTreeReaderRoot6 {
public:
  DJetTreeReaderRoot6(TTree* tree);

  void AddJetDef(const char* jetName);

  bool Next();

  void Restart();

  Long64_t GetEntries() const  { return  fReader.GetEntries(kTRUE); }

  D                                           fDMeson        ;
  std::map<std::string, DJetObjectInfo>       fJets          ;

protected:
  TTree                                      *fTree          ;
  TTreeReader                                 fReader        ;
  TTreeReaderValue<S>                         fDMesonReader  ;
  std::map<std::string,
  std::pair<TTreeReaderValue<J>*,
  DJetObjectInfo* > >                         fJetReaders    ;
};


template<class S, class D>
void FillDJetObjectInfo(const S* source, D* dest);

template<>
void FillDJetObjectInfo(const AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary* source, DJetObjectInfo* dest);

template<>
void FillDJetObjectInfo(const AliAnalysisTaskDmesonJets::AliD0InfoSummary* source, DJetObjectInfo* dest);

template<>
void FillDJetObjectInfo(const AliAnalysisTaskDmesonJets::AliDStarInfoSummary* source, DJetObjectInfo* dest);

template<>
void FillDJetObjectInfo(const J* source, DJetObjectInfo* dest);

template<>
void FillDJetObjectInfo(const AliAnalysisTaskDmesonJetsDetectorResponse::AliD0MatchInfoSummary* source,
    DJetObjectInfoResponse* dest);

template<>
void FillDJetObjectInfo(const AliAnalysisTaskDmesonJetsDetectorResponse::AliDStarMatchInfoSummary* source,
    DJetObjectInfoResponse* dest);
