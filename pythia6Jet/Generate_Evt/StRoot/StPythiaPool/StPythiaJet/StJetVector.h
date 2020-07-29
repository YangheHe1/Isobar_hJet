#ifndef STJETVECTOR_H
#define STJETVECTOR_H
#include "TObject.h"
#include "TLorentzVector.h"
//#include "TParticle.h"
#include "TParticlePDG.h"
#include "StJetFinder/AbstractFourVec.h"

class StJetVector : public AbstractFourVec{  
 public:
  
  // StJetVector() : _momentum(0,0,0,0), _position(0, 0, 0, 0), _pdg(0), _status(0), _mother(0), _daughter1(0), _daughter2(0), _charge(0){}
  StJetVector(){}
  StJetVector(int pp, int ss, int m, int d1, int d2, double px, double py, double pz, double e, double x, double y, double z, double t){
    _momentum.SetPxPyPzE(px, py, pz, e);
    _position.SetPxPyPzE(x, y, z, t);
    _pdg = pp;
    _status = ss;
    _mother = m;
    _daughter1 = d1;
    _daughter2 = d2;
    TParticlePDG tmp(_pdg);
    _charge = (tmp.Charge())/3; // in units of |e|/3
  }
~StJetVector(){};
  double px() const { return _momentum.Px();}
  double py() const { return _momentum.Py();}
  double pz() const { return _momentum.Pz();}

  double phi() const { return _momentum.Phi();}
  double eta() const { 
    //    return _momentum.Eta();
    if(_momentum.Pt() > 0.){ return _momentum.Eta();}
    else if(_momentum.Pz() > 0.){ return 999.;}
    return -999.0;
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

 private:
  TLorentzVector _momentum;
  TLorentzVector _position;
  int _pdg;
  int _status;
  int _mother;
  int _daughter1;
  int _daughter2;
  double _charge;

  //  ClassDef(StJetVector, 1);
  
};
//ClassImp(StJetVector);
#endif  
  
