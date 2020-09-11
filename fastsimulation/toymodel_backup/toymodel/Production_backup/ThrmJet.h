#ifndef __ThrmJet__hh
#define __ThrmJet__hh

#include "TClonesArray.h"
#include "ThrmFourVector.h"

#include "Rtypes.h"

class FourVector;

class ThrmJet : public TObject
{
 public:
  ThrmJet();
  ~ThrmJet();

  Int_t embeddedJet_idx;
  
  Int_t Nconst;
  Double32_t area;      //[0, 0, 16]
  Double32_t pTleading; //[0, 0, 16]

  ThrmFourVector jet_fv;

  ClassDef(ThrmJet, 1);
};

#endif
