#include "ThrmEvent.h"

#include "ThrmFourVector.h"

ClassImp(ThrmEvent)

//________________________________________________________
ThrmEvent::ThrmEvent():
  particles("ThrmFourVector", 1000)
{
  // default constructor
}

//________________________________________________________
ThrmEvent::~ThrmEvent()
{
  // default destructor
  particles.Delete();
}

//________________________________________________________
void ThrmEvent::FreeMemory()
{
  particles.Delete();
}

//________________________________________________________
void ThrmEvent::AddParticleAt(Double_t pT, Double_t eta, Double_t phi, Double_t M, Int_t i)
{
  new (particles[i]) ThrmFourVector(pT, eta, phi, M);
}
