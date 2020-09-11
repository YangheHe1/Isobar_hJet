#ifndef __ThrmFourVector__hh
#define __ThrmFourVector__hh

#include "Rtypes.h"
#include "TLorentzVector.h"

class ThrmFourVector : public TObject
{
 public:
  ThrmFourVector();
  ThrmFourVector(Double_t pT, Double_t eta, Double_t phi, Double_t M);
  ~ThrmFourVector();

  void SetPtEtaPhiM(Double_t pT, Double_t eta, Double_t phi, Double_t M);
  
  TLorentzVector GetTLorentzVector();
  
 protected:
  Double32_t fM;   //[0, 0, 16]
  Double32_t fpT;  //[0, 0, 16]
  Double32_t feta; //[0, 0,  8]
  Double32_t fphi; //[0, 0,  8]

  ClassDef(ThrmFourVector, 1);
};

#endif
