// Helper class to compare results from DJetCorrAnalysis
// Author: Salvatore Aiola, Yale University (salvatore.aiola@cern.ch)
// Copyright (c) 2015 Salvatore Aiola

class DJetCorrAnalysis;

class DJetCorrAnalysisComparer : public TNamed {
  
 public:
  DJetCorrAnalysisComparer();
  DJetCorrAnalysisComparer(DJetCorrAnalysis* ana1, DJetCorrAnalysis* ana2);

  void Start();
  
 protected:
  DJetCorrAnalysis*        fAnalysis1; //! analysis 1
  DJetCorrAnalysis*        fAnalysis2; //! analysis 2
  
 private: 
  DJetCorrAnalysisComparer(const DJetCorrAnalysisComparer &source);
  DJetCorrAnalysisComparer& operator=(const DJetCorrAnalysisComparer& source); 

  ClassDef(DJetCorrAnalysisComparer, 1);
};
