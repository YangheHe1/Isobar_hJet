#ifndef STPYTHIAJETSTATUSCUT
#define STPYTHIAJETSTATUSCUT

#include "TObject.h"
#include "TLorentzVector.h"
#include "TPythia6.h"
#include <iostream>
using std::cout;
using std::endl;
class StPythiaJetStatusCut : public TObject{
 public:
  virtual bool operator()(TPythia6 *pythia, int ip) const = 0;
 private:
  ClassDef(StPythiaJetStatusCut, 1);
};

class StPythiaJetParticleCut : public StPythiaJetStatusCut{
 public:
  StPythiaJetParticleCut(){}
  bool operator()(TPythia6 *pythia, int ip) const
  {
    //    TParticle *par = pythia->GetParticle(ip);
    int status = pythia->GetK(ip+1, 1);
    if(status != 1) return true;
    else return false;
  }

 private:
  ClassDef(StPythiaJetParticleCut, 1);
};

class StPythiaJetPartonCut : public StPythiaJetStatusCut{
 public:
  StPythiaJetPartonCut(){}
  bool operator()(TPythia6 *pythia, int ip) const
  {
    int mstu72 = pythia->GetMSTU(72);
    int mstu73 = pythia->GetMSTU(73);

    //    TParticle *par = pythia->GetParticle(ip);
    int fm = pythia->GetK(ip+1, 3);
    int status = pythia->GetK(ip+1, 1);
    TLorentzVector p(pythia->GetP(ip+1, 1), 
		   pythia->GetP(ip+1, 2), 
		   pythia->GetP(ip+1, 3), 
		   pythia->GetP(ip+1, 4));
    double pt = p.Pt();
    double eta = -999.;
    if(pt > 1e-9) eta = p.Eta();
//    if(p.Px()*p.Px()+p.Py()*p.Py() > 0.) eta = p.Eta();
    else if(p.Pz() > 0.) eta = 999;
    //if(eta > 10 || eta < -10) cout<<pt<<" "<<eta<<endl;

    bool flag = 0;
    
    flag =
      ip < mstu72
      || ip >= mstu73
      || fm > mstu72
      || fm == 0
      || fm == 1
      || fm == 2
      || status == 51
      || pt < 0.0001
      || (eta > 5. || eta < -5.);
      //if(flag) cout<<"ip:"<<ip<<" mstu72:"<<mstu72<<" mstu73:"<<mstu73<<" fm:"<<fm<<" pdg:"<<pythia->GetK(ip+1, 2)<<endl;
    return flag;
    
  }
 private:
  ClassDef(StPythiaJetPartonCut, 1);
};
#endif
