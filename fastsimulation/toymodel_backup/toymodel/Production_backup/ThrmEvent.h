#ifndef __ThrmEvent__hh
#define __ThrmEvent__hh

#include "Rtypes.h"
#include "TClonesArray.h"

class TClonesArray;

class ThrmEvent
{
 public:
  ThrmEvent();
  ~ThrmEvent();
  void FreeMemory();
  void AddParticleAt(Double_t pT, Double_t eta, Double_t phi, Double_t M, Int_t i);

 protected:
  TClonesArray particles;

  ClassDef(ThrmEvent, 1);
};

#endif
