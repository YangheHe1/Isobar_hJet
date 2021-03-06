#ifndef STPYTHIAJETCANDIDATE
#define STPYTHIAJETCANDIDATE
#include "TObject.h"
//#include "StJetFinder/StProtoJet.h"
#include "StPythiaJetParticle.h"
#include "TClonesArray.h"

class StPythiaJetCandidate: public TObject{
 public:
  StPythiaJetCandidate();
//  StPythiaJetCandidate(const StProtoJet *protojet);
  ~StPythiaJetCandidate(){}
  double pt() const{ return _pt;};
  double px() const{ return _px;};
  double py() const{ return _py;};
  double pz() const{ return _pz;};
  double phi() const{ return _phi;}
  double eta() const{ return _eta;}

  double area() const {return _area;}
  double area_error() const { return _areaError;}

  double e() const { return _e;}
  double eT() const { return _eT;}
  double mass() const { return _mass;}
  double charge() const { return _charge;}

  void setPt(double t_pt) { _pt = t_pt;}
  void setPx(double t_px) { _px = t_px;}
  void setPy(double t_py) { _py = t_py;}
  void setPz(double t_pz) { _pz = t_pz;}
  void setPhi(double t_phi) { _phi = t_phi;}
  void setEta(double t_eta) { _eta = t_eta;}

  void setArea(double t_area) { _area = t_area;}
  void setAreaError(double t_areaError) { _areaError = t_areaError;}

  void setE(double t_e) { _e = t_e;}
  void setEt(double t_eT) { _eT = t_eT;}
  void setMass(double t_mass) { _mass = t_mass;}
  void setCharge(double t_charge) { _charge = t_charge;}

  void clear();
  TClonesArray *getParticles() const{ return particles;}
  StPythiaJetParticle *constructParticle() { 
    int ii = particles->GetEntriesFast();
    return (StPythiaJetParticle*)particles->ConstructedAt(ii);
  }
  //  StPythia *GetParticles() const{ particles;}
  //  void fill(StProtoJet protojet);
 private:
  double _pt;
  double _px;
  double _py;
  double _pz;
  double _eta;
  double _phi;

  double _area;
  double _areaError;

  double _e;
  double _eT;
  double _mass;
  double _charge;

  TClonesArray *particles;
  ClassDef(StPythiaJetCandidate, 2);
};
//ClassImp(StPythiaJetCandidate);
#endif
