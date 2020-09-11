#include "ThrmFourVector.h"

//____________________________________________________________________________
ThrmFourVector::ThrmFourVector()
{
  fM = 0.;
  fpT = 0.;
  feta = 0.;
  fphi = 0.;
}

//____________________________________________________________________________
ThrmFourVector::ThrmFourVector(Double_t pT, Double_t eta, Double_t phi, Double_t M)
{
  fM = M;
  fpT = pT;
  feta = eta;
  fphi = phi;
}

//____________________________________________________________________________
ThrmFourVector::~ThrmFourVector()
{
}

//____________________________________________________________________________
void ThrmFourVector::SetPtEtaPhiM(Double_t pT, Double_t eta, Double_t phi, Double_t M)
{
  fM = M;
  fpT = pT;
  feta = eta;
  fphi = phi;
}

//____________________________________________________________________________
TLorentzVector ThrmFourVector::GetTLorentzVector()
{
  TLorentzVector lv;
  lv.SetPtEtaPhiM(fpT, feta, fphi, fM);
  
  return lv;
}
