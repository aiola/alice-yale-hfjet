// Helper class to compare results from DJetCorr*
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

#include "DJetCorrBase.h"

class TObjArray;
class TArrayI;

class DJetCorrAnalysisComparer : public DJetCorrBase {
  
 public:
  enum ECompareTask { kNone=0, kCompareTruth=1<<1, kCompareMeasured=1<<2 };
  enum ENormalizationType { kNone=0, kIntegral=1, kEvents=2 };

  DJetCorrAnalysisComparer();
  DJetCorrAnalysisComparer(UInt_t task);
  DJetCorrAnalysisComparer(UInt_t task, DJetCorrBase* ana1, DJetCorrBase* ana2, Int_t ipar1=0, Int_t ipar2=0);

  void SetCompareTask(UInt_t task)        { fCompareTask       = task; }
  void SetMakeRatios(Bool_t r)            { fMakeRatios        = r   ; }
  void SetNormalizationType(UInt_t task)  { fNormalizationType = task; }

  void Start();
  Bool_t Prepare();
  Bool_t ExecuteTasks();

  void CompareTruth();
  void CompareMeasured();
  void CompareZvsJetPt(const char* name, TObjArray& array);
  void CompareDPt(const char* name, TObjArray& array);

  DJetCorrBase* AddAnalysis(DJetCorrBase* ana, Int_t ipar=0);
  
 protected:
  TArrayI              fParamIndexes     ; //
  Bool_t               fForceRegeneration; //
  Bool_t               fMakeRatios       ; //
  UInt_t               fCompareTask      ; //
  ENormalizationType   fNormalizationType; //

  TObjArray*        fAnalysisArray    ; //! analysis array
  
 private: 
  DJetCorrAnalysisComparer(const DJetCorrAnalysisComparer &source);
  DJetCorrAnalysisComparer& operator=(const DJetCorrAnalysisComparer& source); 

  ClassDef(DJetCorrAnalysisComparer, 1);
};
