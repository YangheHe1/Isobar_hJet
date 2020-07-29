#include<iostream>
using namespace std;

#include "histPythia.h"

#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

histPythia:: histPythia()
{
  TH1::StatOverflows(kTRUE);
  TH1::SetDefaultSumw2(kTRUE);

  const char *par[Npar] = {"Parton", "Particle"};
  const char *cone[Ncone] = {"minus", "plus", "sum"};
  const char *proc[Nproc] = {"gg", "qg", "qq"};

  for(int ipar = 0; ipar < Npar; ipar++){
    h[ipar] = new TH1F(Form("%s_pt", par[ipar]), Form("%s_pt", par[ipar]),  100, -0.5, 99.5);
    h_2d[ipar] = new TH2F(Form("%s_area_pt", par[ipar]), Form("%s_area_pt", par[ipar]), 100, -0.5, 99.5, 180, 0, 1.8);
  }
  h_matched = new TH1F(Form("%sMatched%s_pt", par[1], par[0]), Form("%sMatched%s_pt", par[1], par[0]),  100, -0.5, 99.5);
  h_matched_prof = new TProfile(Form("%sMatched%s_pt_prof", par[1], par[0]), Form("%sMatched%s_pt; jet p_{T} [GeV]; ratio", par[1], par[0]),  100, -0.5, 99.5);
  h_matched_2d = new TH2F(Form("%sMatched%s_area_pt", par[1], par[0]), Form("%sMatched%s_area_pt", par[1], par[0]), 100, -0.5, 99.5, 180, 0, 1.8);
  
  for(int icone = 0; icone < Ncone; icone++){
    hue_m[icone] = new TH1F(Form("%sMatched%s_ue_pt_%s", par[1], par[0], cone[icone]), Form("%sMatched%s_ue_pt_%s", par[1], par[0], cone[icone]), 160,-1, 15);
    hue_unm[icone] = new TH1F(Form("%sUnmatched%s_ue_pt_%s", par[1], par[0], cone[icone]), Form("%sUnmatched%s_ue_pt_%s", par[1], par[0], cone[icone]), 160,-1, 15);
    hue_m_prof[icone] = new TProfile(Form("%sMatched%s_ue_pt_prof_%s", par[1], par[0], cone[icone]), ";jet p_{T} [GeV];UE cone sum p_{T} [GeV]", nbins, ptbins);
    hue_unm_prof[icone] = new TProfile(Form("%sUnmatched%s_ue_pt_prof_%s", par[1], par[0], cone[icone]), ";jet p_{T} [GeV]; UE cone sum p_{T} [GeV]", nbins, ptbins);
  }
  hue_m_2d = new TH2F(Form("%sMatched%s_ue_pt_2d", par[1], par[0]), Form("%sMatched%s_ue_pt_2d;%s;%s", par[1], par[0], cone[0], cone[1]), 160, -1, 15, 160, -1, 15);
  hue_unm_2d = new TH2F(Form("%sUnmatched%s_ue_pt_2d", par[1], par[0]), Form("%sUnmatched%s_ue_pt_2d;%s;%s", par[1], par[0], cone[0], cone[1]), 160, -1, 15, 160, -1, 15);
  

  //pythia record partonic pT
  h_partonic = new TH1F("partonic_pt", "partonic_pt", 5000, 0, 100);

  h_rmin = new TH1F("particle_rmin", "particle_rmin", 1000, 0, 10);
  h_rmin_pt = new TH2F("particle_rmin_pt", "particle_rmin_pt", 30, 0, 30, 1000, 0, 10);
  //subprocess
  for(int ipar = 0; ipar < Npar; ipar++){
      for(int iproc = 0; iproc < Nproc; iproc++){
	hproc[ipar][iproc] = new TH1F(Form("hpt_%s_%s", par[ipar], proc[iproc]), Form("hpt_%s_%s;p_{T} [GeV]", par[ipar], proc[iproc]), 100, -0.5, 99.5);
      }
    }
    //jet dpt profile
    const char *name = "Particle";
    hjetpt = new TH1F(Form("%sJetPt", name), ";p_{T} [GeV]",nbins,ptbins);
    hjeteta = new TH1F(Form("%sJetEta", name), ";#eta", 40, -1., 1.);
    hjetphi = new TH1F(Form("%sJetPhi", name), ";#phi", 126, -3.15, 3.15);
    hjetmult = new TH2F(Form("%sJetMult", name), ";mult", nbins, ptbins, 100, 0, 100);
    
    hjetht = new TH2F(Form("%sJetHt", name), ";H_{T}", nbins, ptbins, 100, 0., 1.);
    hjetparjt = new TH2F(Form("%sJetParticleJt", name), ";j_{T}", nbins, ptbins, 100, 0., 10.);
    hjetparz = new TH2F(Form("%sJetParticleZ", name), ";z", nbins, ptbins, 100, 0., 1.);

    hjetparmult = new TH2F(Form("%sJetHadronMult", name), ";track mult", nbins, ptbins, 100, 0, 100);
    //pt shift
    hjetmatched = new TH1F(Form("%sJetMatchedPt", name), "; matched p_{T} [GeV]", nbins, ptbins);
    hjetmatchedratio = new TProfile(Form("%sJetMatchedRatio", name), "; p_{T} [GeV]; matched ratio", nbins, ptbins);
    hjetmatchedratiocrr = new TProfile(Form("%sJetMatchedRatioCrr", name), "; p_{T} [GeV]; matched ratio crr.", nbins, ptbins);
    hjetshift = new TProfile(Form("%sJetPtShift", name), "; p_{T} [GeV]; p_{T} shift", nbins, ptbins);
    hjetshiftcrr = new TProfile(Form("%sJetPtShiftCrr", name), "; p_{T} [GeV]; p_{T} shift crr.", nbins, ptbins);
    //sub process
    for(int iproc = 0; iproc < Nproc; iproc++){
      hjetproc[iproc] = new TH1F(Form("%sJet%sPt", name, proc[iproc]), ";p_{T} [GeV]", nbins, ptbins);
    }
    //ue
    hjetuedpt = new TProfile(Form("%sJetUeDpt", name), ";p_{T} [GeV]; ue dp_{T} [GeV]", nbins, ptbins);
    hjetuemult = new TH2F(Form("%sJetUeMult", name), ";mult", nbins, ptbins, 100, 0, 100);
    
    hjetueht = new TH2F(Form("%sJetUeHt", name), ";H_{T}", nbins, ptbins, 100, 0., 1.);
    hjetueparjt = new TH2F(Form("%sJetUeParticleJt", name), ";j_{T}", nbins, ptbins, 100, 0., 10.);
    hjetueparz = new TH2F(Form("%sJetUeParticleZ", name), ";z", nbins, ptbins, 100, 0., 1.);
    hjetueparmult = new TH2F(Form("%sJetUeHadronMult", name), ";track mult", nbins, ptbins, 100, 0, 100);
    double binT[40+1];
    for(int ib = 0; ib < 40+1; ib++){
      binT[ib] = -1. + ib*1.;
    }
    hjetuedpt_3d = new TH3F(Form("%sJetUeDpt3D", name), ";jet p_{T} [GeV]; dp_{T,1} [GeV]; dp_{T,2} [GeV]", nbins, ptbins, 40, binT, 40, binT);
}
void histPythia::Add(const histPythia *hist, double w)
{
  h_partonic->Add(hist->h_partonic, w);
  for(int i = 0; i < Npar; i++){
    h[i]->Add(hist->h[i], w);
    h_2d[i]->Add(hist->h_2d[i], w);
  }
  h_matched->Add(hist->h_matched, w);
  h_matched_prof->Add(hist->h_matched_prof, w);
  h_matched_2d->Add(hist->h_matched_2d);
  for(int i = 0; i < Npar; i++){
    for(int j = 0; j < Nproc; j++){
      hproc[i][j]->Add(hist->hproc[i][j], w);
    }
  }
  for(int i = 0; i < Ncone; i++){
    hue_m[i]->Add(hist->hue_m[i], w);
    hue_unm[i]->Add(hist->hue_unm[i], w);
    hue_m_prof[i]->Add(hist->hue_m_prof[i], w);
    hue_unm_prof[i]->Add(hist->hue_unm_prof[i], w);
  }
  hue_m_2d->Add(hist->hue_m_2d, w);
  hue_unm_2d->Add(hist->hue_unm_2d, w);
  h_rmin->Add(hist->h_rmin, w);
  h_rmin_pt->Add(hist->h_rmin_pt, w);

  //jet
   hjetpt->Add(hist->hjetpt, w);
   hjeteta->Add(hist->hjeteta, w);
   hjetphi->Add(hist->hjetphi, w);
   hjetmult->Add(hist->hjetmult, w);
   hjetht->Add(hist->hjetht, w);
   hjetparjt->Add(hist->hjetparjt, w);
   hjetparz->Add(hist->hjetparz, w);
   hjetparmult->Add(hist->hjetparmult, w);
  //pt shift
   hjetmatched->Add(hist->hjetmatched, w);
   hjetmatchedratio->Add(hist->hjetmatchedratio, w);
   hjetmatchedratiocrr->Add(hist->hjetmatchedratiocrr, w);
   hjetshift->Add(hist->hjetshift, w);
   hjetshiftcrr->Add(hist->hjetshiftcrr, w);
  //sub process
   for(int i = 0; i < Nproc; i++)
     hjetproc[i]->Add(hist->hjetproc[i], w);
  //ue
   hjetuedpt->Add(hist->hjetuedpt, w);
   hjetuemult->Add(hist->hjetuemult, w);
   hjetueht->Add(hist->hjetueht, w);
   hjetueparjt->Add(hist->hjetueparjt, w);
   hjetueparz->Add(hist->hjetueparz, w);
   hjetueparmult->Add(hist->hjetueparmult, w);
   hjetuedpt_3d->Add(hist->hjetuedpt_3d, w);
}
