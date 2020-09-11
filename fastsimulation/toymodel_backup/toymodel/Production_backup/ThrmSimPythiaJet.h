#ifndef __ThrmSimPythiaJet__hh
#define __ThrmSimPythiaJet__hh

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TPythia6.h"
#include "TClonesArray.h"

class TF1;
class TFile;
class TTree;
class TClonesArray;

class ThrmSimPythiaJet
{
 public:
  ThrmSimPythiaJet(TString wrkdir, bool jetonly=false, bool charge=false, bool efficorr=false, bool pTsmear=false);
  ~ThrmSimPythiaJet();
  
  void Run();

  // SETTERS
  void SetNevents(Int_t Nevents) {fNevents = Nevents;}
  void SetNbinaryColl(Int_t nbinary = 0) {fNbinaryColl = nbinary;}
  void SetMultiplicity(Int_t multiplicity = 0) {fMultiplicity = multiplicity;}

  void SetSigmaNbinaryColl(Int_t sigmaNbinary = 0) {fSigmaNbinaryColl = sigmaNbinary;}
  void SetSigmaMultiplicity(Int_t sigmaMultiplicity = 0) {fSigmaMultiplicity = sigmaMultiplicity;}
  void SetJetsOnly() {kJetOnly = kTRUE;}
  void SetThermalOnly() {kThrmOnly = kTRUE;}
  void SetKinematics(TString collider) {fkinematics = collider;}
  
 private:
  void ConfigEvent(Int_t &multiplicity);
  void ProduceParticles(Int_t multiplicity);
  
 protected:
  TTree *ftree;
  TFile *foutput;
  TFile *foutput_histo;
  TFile* efffile;

  TPythia6 *fpythia;

  TClonesArray fpartarr;

  TF1 *fjet;
  TF1 *fbkgd;
  TF1 *fspectrum;
  TF1* effL;
  TF1* effH;

  TH1D *fhalldNdpT;
  TH1D *fhjetpartdNdpT;
  TH1D *fhjetdNdpT;
  TH1D *fhjetreqpT;
  TH1D *fhboltzdNdpT;
  TH2D *fhpTpTleading;
  TH2D *fhpTpTleadingGen;

  Bool_t  kCentral;
  Bool_t kFixedSeed;
  Int_t kSeed;
  Int_t fNevents;
  Int_t fNbinaryColl;
  Int_t fMultiplicity;
  Int_t fSigmaNbinaryColl;
  Int_t fSigmaMultiplicity;
  Int_t hjettype;
  Float_t increment;
  Float_t fR;
  Float_t RAA;
  Float_t pTMAX;
  Float_t pTMIN;
  Float_t pTMINhard;
  Float_t meanpT;

  Bool_t kJetOnly;
  Bool_t kThrmOnly;
  Bool_t kCharged;
  Bool_t kEfficorr;
  Bool_t kpTsmear;
  Bool_t RAApTdep;

  TString fkinematics;
  TString frag;
  TString eff_type;
  TString eff_path;
  TString foutput_histo_path;
  TString foutput_path;

};

#endif
