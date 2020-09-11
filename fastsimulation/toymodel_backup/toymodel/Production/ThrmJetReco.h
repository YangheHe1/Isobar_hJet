#ifndef __ThrmJetReco__hh
#define __ThrmJetReco__hh

#include "TString.h"

class ThrmJetReco
{
 public:
  ThrmJetReco(TString wrkdir, Double_t rparam, Double_t pTcut);
  ~ThrmJetReco();

  void SetRparam(Double_t rparam) {frparam = rparam;}
  void SetPtCut(Double_t pTcut) {fpTcut = pTcut;}

  void RunJetReco();

 protected:
  Double_t frparam;
  Double_t fpTcut;

  TString fwrkdir;
};

#endif
