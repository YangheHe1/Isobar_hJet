#ifndef __ThrmRecPythia__hh
#define __ThrmRecPythia__hh

#include "TH1D.h"
#include "TPythia6.h"
#include "TClonesArray.h"

/*
Double_t Eff_track_rec_function(Double_t* x,Double_t* par);
Double_t RatioTsalis(Double_t* x_val, Double_t* par);
Double_t FitFunc(Double_t* x_val, Double_t* par);
Double_t efficiencyAlex(Double_t pt,Double_t increment=0,bool kcentral=true,TString ratio="pp");
Double_t efficiency11(Double_t pt, TF1* effLow, TF1* effHigh, Double_t increment=0, bool kcentral=true);
double efficiencyTOFvBEMC(double pt,double increment=0);
*/
//Double_t efficiency11H(Double_t pT, bool kcentral=true);

Double_t MakeCone(TPythia6 *pythia, TClonesArray *detearr, TClonesArray *pararr, Double_t pT, Double_t jetEta, Double_t jetPhi, Bool_t charged="kFALSE", Bool_t effi="kFALSE", Bool_t pTsmear="kFALSE", TString frag="u",bool fillEffHisto=false, TH1D* hTrackEff=NULL,TH1D* hparticle=NULL);

Double_t MakeCone2(TPythia6 *pythia, TClonesArray *detearr, Double_t pT, Double_t jetEta, Double_t jetPhi, Bool_t charged="kFALSE", Bool_t effi="kFALSE", Bool_t pTsmear="kFALSE", TString frag="u",bool fillEffHisto=false, TH1D* hTrackEff=NULL);
#endif
