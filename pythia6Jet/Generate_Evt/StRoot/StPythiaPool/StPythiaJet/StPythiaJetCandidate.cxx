#include "StPythiaJetCandidate.h"
//#include "StJetVector.h"
//#include "StPythiaJetParticle.h"
//#include <iostream>
//using namespace std;
StPythiaJetCandidate::StPythiaJetCandidate()
{
  _pt = 0;
  _px = 0;
  _py = 0;
  _pz = 0;

  _eta = 0;
  _phi = 0;

  _area = 0.;
  _areaError = 0.;

  _e = 0;
  _eT = 0;
  _mass = 0;
  _charge = 0;

  particles = new TClonesArray("StPythiaJetParticle", 50);
}
/*
StPythiaJetCandidate::StPythiaJetCandidate(const StProtoJet *protojet)
{
  _pt = protojet->pt();
  _px = protojet->px();
  _py = protojet->py();
  _pz = protojet->pz();

  _eta = -999.;
  if(protojet->pt() > 0.) _eta = protojet->eta();
  else if(protojet->pz() > 0.) _eta = 999.;
   
  _phi = protojet->phi();

  _e = protojet->e();
  _eT = protojet->eT();
  _mass = protojet->mass();
  _charge = protojet->charge();

  particles = new TClonesArray("StPythiaJetParticle", 50);
  TClonesArray &array = *particles;
//  for (StProtoJet::FourVecList::const_iterator ipar = protojet->list().begin(); ipar != protojet->list().end(); ipar++){
//    array->Add(new StPythiaJetParticle((StJetVector *)*ipar));
//  }
  for(size_t ii = 0; ii < protojet->list().size(); ii++){
    StJetVector* vec = (StJetVector*)(protojet->list()[ii]);
    //    cout<<"px="<<vec->px()<<endl;
    //    StPythiaJetParticle *p = new StPythiaJetParticle(vec);
    //    array->At(ii) = p;
    new(array[ii]) StPythiaJetParticle(vec);
    //    cout<<"ientry:"<<array.GetEntriesFast()<<endl;  
  }
  //  cout<<"entry:"<<particles->GetEntries()<<endl;  
  //  cout<<"entry:"<<array.GetEntries()<<endl;  
}
*/
void StPythiaJetCandidate::clear()
{
  particles->Clear();
}
