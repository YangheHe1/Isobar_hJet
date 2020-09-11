#include "TLorentzVector.h"

#include "ThrmJet.h"

ClassImp(ThrmJet)

//________________________________________________________
ThrmJet::ThrmJet():
  jet_fv()
{
  // default constructor
  embeddedJet_idx = -1;
  Nconst = 0.;
  area = 0.0;
  pTleading = 0.0;
}

//________________________________________________________
ThrmJet::~ThrmJet()
{
  // default destructor
}

