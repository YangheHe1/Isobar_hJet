#include "ThrmJetReco.h"
#include "JetReconstruction.h"

//==========================================================================
ThrmJetReco::ThrmJetReco(TString wrkdir, Double_t rparam, Double_t pTcut)
{
  fwrkdir = wrkdir;
  frparam = rparam;
  fpTcut = pTcut;
}

//==========================================================================
ThrmJetReco::~ThrmJetReco()
{
  
}

//==========================================================================
void ThrmJetReco::RunJetReco()
{
  JetReconstruction *reco = new JetReconstruction(fwrkdir, frparam, fpTcut);
  reco->Run();
  delete reco;
}

