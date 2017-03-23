#include <map>
#include <string>

#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <AliAnalysisTaskDmesonJets.h>

class TTree;

typedef AliAnalysisTaskDmesonJets::AliJetInfoSummary J;

struct DJetInfo {
  double  fPt ;
  double  fEta;
  double  fPhi;
  double  fR  ;
  double  fZ  ;
  int     fN  ;
};

struct DMesonInfo {
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
};

template<class D>
class DJetTreeReaderRoot6 {
public:
  DJetTreeReaderRoot6(TTree* tree);

  void AddJetDef(const char* jetName);

  bool Next();

  void Restart();

  DMesonInfo                                         fDMeson        ;
  std::map<std::string, DJetInfo>                    fJets          ;

protected:
  TTree                                             *fTree          ;
  TTreeReader                                        fReader        ;
  TTreeReaderValue<D>                                fDMesonReader  ;
  std::map<std::string, std::pair<TTreeReaderValue<J>*,  DJetInfo* > > fJetReaders;
};


template<class D>
void FillDMesonInfo(const D* source, DMesonInfo* dest);

template<>
void FillDMesonInfo(const AliAnalysisTaskDmesonJets::AliDmesonMCInfoSummary* source, DMesonInfo* dest);

template<>
void FillDMesonInfo(const AliAnalysisTaskDmesonJets::AliD0InfoSummary* source, DMesonInfo* dest);

template<>
void FillDMesonInfo(const AliAnalysisTaskDmesonJets::AliDStarInfoSummary* source, DMesonInfo* dest);

void FillDJetInfo(const J* source, DJetInfo* dest);
