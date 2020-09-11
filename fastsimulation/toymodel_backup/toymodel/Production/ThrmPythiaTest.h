#ifndef __ThrmPythiaTest__hh
#define __ThrmPythiaTest__hh

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TPythia6.h"
#include "TClonesArray.h"

class TF1;
class TFile;
class TTree;
class TClonesArray;

class ThrmPythiaTest
{
 public:
  ThrmPythiaTest();
  ~ThrmPythiaTest();

	void Run(Int_t nevents);
	double getDistance(Double_t eta1, Double_t eta2,Double_t phi1, Double_t phi2);

 protected:
ThrmMemStat memstat;

Float_t fR; 
Float_t fmax_rap; 
Float_t fpTcut;
Float_t fAreaCut;
Bool_t kCharged; 
Int_t kEfficorr;
Int_t kpTsmear;
Float_t etaMinCut;
Float_t etaMaxCut;
TString frag;
float pTmin;
float pTmax;
int nptbins;
int fNjobs;
int fweight;

TFile *foutput;
TTree *ftreejets;

TClonesArray fembeddingarr;

TPythia6 *fpythia;

TF1 *effL;
TF1 *effH;
TF1 *fjet;

static const Int_t npTlead=11;
TH1D* hpTptlJet[npTlead];
TH1D* hpTptlJet0[npTlead];
TH1D* hpTptl4v[npTlead];
TH1D* hpTjet;
TH1D* hpt4v1;
TH1D* hpt4v2;
TH2D* hdpt_nparts;
TH2D* hdpt_pt;
TH2D* hdpt_z;
TH2D* hz_pT;
TH1D* hnJets;
TH2D* hpT_pTdete_pTl5;

};
#endif
