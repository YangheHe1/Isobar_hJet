#ifndef __ThrmPythiaEmb__hh
#define __ThrmPythiaEmb__hh

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TString.h"
#include "TPythia6.h"
#include "TClonesArray.h"

class TF1;
class TFile;
class TTree;
class TClonesArray;

class ThrmPythiaEmb
{
 public:
  ThrmPythiaEmb();
  ~ThrmPythiaEmb();

	void Run(Int_t nevents);
  //void read3dEff();
	double getDistance(Double_t eta1, Double_t eta2,Double_t phi1, Double_t phi2);

 protected:
ThrmMemStat memstat;

Float_t fR; 
Float_t fmax_rap; 
Float_t fpTcut;
Float_t fAreaCut;
Bool_t kCharged; 
Bool_t kFixedSeed;
Int_t kSeed;
Int_t kEfficorr;
Int_t kpTsmear;
Float_t etaMinCut;
Float_t etaMaxCut;
Bool_t kCen;
TString frag;
float pTmin;
float pTmax;
int nptbins;
int fNjobs;

TFile *foutput;
TTree *ftreejets;

TClonesArray fembeddingarr;

TPythia6 *fpythia;

TF1 *effL;
TF1 *effH;
TF1 *fjet;

TH3D* h3dEff;

static const Int_t npTlead=7;
//TH1D* heffi_par[npTlead];
//TH1D* heffi_dete[npTlead];
TH1D* hparticle;
TH1D* hTrackEff;
TH1D* hTrackEpsilon;
TH1D* hJetEff1D[npTlead];
TH1D* hRM_projectionY[npTlead];
TH2D* hJetEff2D[npTlead];
TH2D* hresponse[npTlead];
TH1D* hntrue_rm[npTlead];
TH1D* hndete_rm[npTlead];
TH1D* hnpart_rm[npTlead];
TH1D* hntrue_all[npTlead];
TH2D* hJER[npTlead];
TH2D* hntruerm_vs_eta[npTlead];
TH2D* hndeterm_vs_eta[npTlead];
TH2D* hntruerm_vs_phi[npTlead];
TH2D* hndeterm_vs_phi[npTlead];
TH1D* hntrue_all_eta[npTlead];
TH1D* hntrue_rm_eta[npTlead];
TH1D* hntrue_all_phi[npTlead];
TH1D* hntrue_rm_phi[npTlead];
//TH2D* hresponse_0;

};
#endif
