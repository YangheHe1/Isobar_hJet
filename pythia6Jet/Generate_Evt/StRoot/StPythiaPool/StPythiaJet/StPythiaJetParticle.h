#ifndef STPYTHIAJETPARTICLE_H
#define STPYTHIAJETPARTICLE_H
#include "TObject.h"
#include "TLorentzVector.h"
//#include "TParticle.h"
//#include "TParticlePDG.h"
//#include "StJetVector.h"

class StPythiaJetParticle : public TObject{  
 public:
  
  // StPythiaJetParticle() : _momentum(0,0,0,0), _position(0, 0, 0, 0), _pdg(0), _status(0), _mother(0 ), _daughter1(0), _daughter2(0), _charge(0){} 
  StPythiaJetParticle(){}
/*
  StPythiaJetParticle(const StJetVector *particle){
    _momentum = particle->momentum();
    _position = particle->position();
    _pdg = particle->pdg();
    _status = particle->status();
    _mother = particle->mother();
    _daughter1 = particle->daughterfirst();
    _daughter2 = particle->daughtersecond();
    _charge = particle->charge();
  }
*/
~StPythiaJetParticle(){};
  double px() const { return _momentum.Px();}
  double py() const { return _momentum.Py();}
  double pz() const { return _momentum.Pz();}

  double pt() const { return _momentum.Pt(); }
  double phi() const { return _momentum.Phi();}
  double eta() const { 
    if(_momentum.Pt() > 0.) return _momentum.Eta();
    else if(_momentum.Pz() > 0.) return 999.;
    return -999.;
  }

  double eT() const { return _momentum.Et();}
  double e() const { return _momentum.E();}
  double mass() const { return _momentum.M();}

  double charge() const { return _charge;}
  TLorentzVector momentum() const { return _momentum;}
  TLorentzVector position() const { return _position;}

  int pdg() const { return _pdg;}
  int status() const { return _status;}
  int mother() const { return _mother;}
  int daughterfirst() const {return _daughter1;}
  int daughtersecond() const {return _daughter2;}

  void setMom(TLorentzVector t_mom) { _momentum = t_mom; }
  void setPos(TLorentzVector t_pos) { _position = t_pos; }
  void setPdg(int t_pdg) { _pdg = t_pdg; }
  void setStatus(int t_status) {  _status = t_status; }
  void setMother(int t_mother) { _mother = t_mother; }
  void setDaughterfirst (int t_df) { _daughter1 = t_df; }
  void setDaughtersecond(int t_ds) { _daughter2 = t_ds; }
  void setCharge(int t_charge) { _charge = t_charge; }
  
 private:
  TLorentzVector _momentum;
  TLorentzVector _position;
  int _pdg;
  int _status;
  int _mother;
  int _daughter1;
  int _daughter2;
  double _charge;

  ClassDef(StPythiaJetParticle, 1);
  
};
//ClassImp(StPythiaJetParticle);
#endif  
  
