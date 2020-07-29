#ifndef HIST_PYTHIA
#define HIST_PYTHIA

#include "TObject.h"

#include "defPythia.h"

class TH1F;
class TH2F;
class TH3F;
class TProfile;

class histPythia : public TObject
{
 public:
  histPythia();
  void Add(const histPythia* hist, double w = 1.0);
  /******/    
  TH1F *h_partonic;
  TH1F *h[Npar];
  TH2F *h_2d[Npar];
  TH1F *h_matched;
  TProfile *h_matched_prof;
  TH2F *h_matched_2d;
  TH1F *hproc[Npar][Nproc];
  TH1F *hue_m[Ncone];
  TH1F *hue_unm[Ncone];
  TProfile *hue_m_prof[Ncone];
  TProfile *hue_unm_prof[Ncone];
  TH2F *hue_m_2d;
  TH2F *hue_unm_2d;
  TH1F *h_rmin;
  TH2F *h_rmin_pt;
  //jet
  TH1F *hjetpt;
  TH1F *hjeteta;
  TH1F *hjetphi;
  TH2F *hjetmult;
  TH2F *hjetht;
  TH2F *hjetparjt;
  TH2F *hjetparz;
  TH2F *hjetparmult;
  //pt shift
  TH1F *hjetmatched;
  TProfile *hjetmatchedratio;
  TProfile *hjetshift;
  TProfile *hjetmatchedratiocrr;
  TProfile *hjetshiftcrr;
  //sub process
  TH1F *hjetproc[Nproc];
  //ue
  TProfile *hjetuedpt;
  TH2F *hjetuemult;
  TH2F *hjetueht;
  TH2F *hjetueparjt;
  TH2F *hjetueparz;
  TH2F *hjetueparmult;
  TH3F *hjetuedpt_3d;
  /******/

 private:
  ClassDef(histPythia, 2);

};
#endif
