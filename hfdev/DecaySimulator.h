class TLorentzVector;
class TTree;

class DecaySimulator : public TNamed {
  
 public:
  DecaySimulator();
  DecaySimulator(Double_t mass, Double_t massDaughter1, Double_t massDaughter2);

  void SetMass(Double_t mass)                             { fMass          = mass  ; }
  void SetDaughterMasses(Double_t mass1, Double_t mass2)  { fDaughterMass1 = mass1 ; fDaughterMass2 = mass2 ; }
  
  void SetPtBins(Int_t nbins, Double_t *bins);
  void Init();
  void StartSimulation(Int_t n = 1000);
  void SaveTree(const char* fileName);
  void CalculateDaughters(Double_t theta, Double_t phi);
  
protected:
  Double_t         fMass            ; //            Mass of the mother particle
  Double_t         fDaughterMass1   ; //            Mass of the daughter 1 particle
  Double_t         fDaughterMass2   ; //            Mass of the daughter 2 particle
  Int_t            fNPtBins         ; //            Number of pt bins of the simulation
  Double_t        *fPtBins          ; //[fNPtBins]  Pt bins of the simulation
  Double_t         fMinEta          ; //            Minimum simulated eta
  Double_t         fMaxEta          ; //            Maximum simulated eta

  Double_t         fCMMomentum      ; //            Momentum of daughters in the CM
  Double_t         fCMEnergy1       ; //            Energy of daughter 1 in the CM
  Double_t         fCMEnergy2       ; //            Energy of daughter 2 in the CM
  TTree           *fDecaySimulation ; //!           Tree of the simulated decays
  TLorentzVector  *fMother          ; //!           Mother particle
  TLorentzVector  *fDaughter1       ; //!           Daughter 1 particle
  TLorentzVector  *fDaughter2       ; //!           Daughter 1 particle
  
 private:
  DecaySimulator(const DecaySimulator &source);
  DecaySimulator& operator=(const DecaySimulator& source); 

  ClassDef(DecaySimulator, 1);
};
