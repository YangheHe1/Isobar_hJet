#include<iostream>
using namespace std;

#include "readPythia.h"
#include "utilPythia.h"

#include "TFile.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TH3F.h"
#include "TRefArray.h"
#include "TParticle.h"
#include "TMath.h"

#include "StPythiaPool/StPythiaJetHist/histPythia.h"
#include "StPythiaPool/StPythiaJetHist/histResponse.h"
#include "StPythiaPool/StPythiaJet/StPythiaJetEvent.h"
#include "StPythiaPool/StPythiaJet/StPythiaJetCandidate.h"
#include "StPythiaPool/StPythiaJet/StPythiaJetParticle.h"
#include "StSpinPool/StJetSkimEvent/StPythiaEvent.h"

readPythia::readPythia(StPythiaEvent *pythia, StPythiaJetEvent *parton, StPythiaJetEvent *particle)
{
  mPythia = pythia;
  mPythiaJet[0] = parton;
  mPythiaJet[1] = particle;
}
int readPythia::Init()
{
  mFile = new TFile(mFilename, "recreate");
  mHist = new histPythia;
  mHistResponse = new histResponse("Pythia");
  return 1;
}
int readPythia::ProcessIndex()
{
  int proc_id = mPythia->processId();
  int proc_index = -1;
  if(proc_id == 11 || proc_id == 12 || proc_id == 13){
      proc_index = 2;
  }else if(proc_id == 28){
      proc_index = 1;
  }else if(proc_id == 53 || proc_id == 68){
      proc_index = 0;
  }else{
      cout<<"Error: wrong sub-process id = "<<proc_id<<endl;
  }
  return proc_index;
}
int readPythia::Make(int iEvent)
{
  cout<<"Event: "<<iEvent<<endl;
  double pt_par = mPythia->pt();
  double ww = partonicWeight(pt_par);
  mHist->h_partonic->Fill(pt_par, ww);
  //Printf("partonic weight: %lf", ww);
  int proc_index = ProcessIndex();
 
  for(int ipar = 0; ipar < Npar; ipar++){
    TClonesArray *array = mPythiaJet[ipar]->getJets();
    //      cout<<par[ipar]<<endl;
    int Njet = array->GetEntriesFast();
    //     cout<<"Njet:"<<Njet<<endl;
    for(int ijet = 0; ijet < Njet; ijet++){
      StPythiaJetCandidate *jetcnd = (StPythiaJetCandidate *)array->At(ijet);
      //        cout<<"jet pt:"<<jetcnd->pt()<<endl;
      double pt = jetcnd->pt();
      double eta = jetcnd->eta();
      double area = jetcnd->area();
      //          cout<<"jet pt = "<<pt<<" area = "<<area<<endl;
      if(eta > -0.9 && eta < 0.9){
	mHist->h[ipar]->Fill(pt, ww);
	mHist->h_2d[ipar]->Fill(pt, area, ww);
	if(proc_index >= 0 && proc_index <= 2){
            mHist->hproc[ipar][proc_index]->Fill(pt, ww);
        }
      }
    }
  }
  //match particle jets to parton jets

  TClonesArray *arr_particle = mPythiaJet[1]->getJets();
  int Nparticle = arr_particle->GetEntriesFast();
  //      cout<<"particles: "<<Nparticle<<endl;
  //float ptmax = 0, dptmax = 0.;
  //bool flagmax = false;
  for(int ijet = 0; ijet < Nparticle; ijet++){
    StPythiaJetCandidate *jetcnd = dynamic_cast<StPythiaJetCandidate *>(arr_particle->At(ijet));
    double eta = jetcnd->eta();
    double pt = jetcnd->pt();
    double phi = jetcnd->phi();
    double area = jetcnd->area();
    if(pt < 6.) continue;
    if(eta < -0.9 || eta > 0.9) continue;
    //jet pt based on partonic pt
    if(pt_par < 3 && pt > 20.9) continue;
    if(pt_par < 4 && pt > 24.5) continue;
    if(pt_par < 7 && pt > 39.3) continue;
    TVector3 mom;
    mom.SetPtEtaPhi(pt, eta, phi);
    //double p = mom.Mag();
    mHist->hjetpt->Fill(pt, ww);
    if(proc_index >= 0 && proc_index <= 2){
      mHist->hjetproc[proc_index]->Fill(pt, ww);
    }
    mHist->hjeteta->Fill(eta, ww);
    mHist->hjetphi->Fill(phi, ww);
    double mult = jetcnd->getParticles()->GetEntriesFast();
    mHist->hjetmult->Fill(pt, mult, ww);
    double sumcharged = 0.;
    int multcharged = 0;
    for(int ip = 0; ip < mult; ip++){
      StPythiaJetParticle *ppar = (StPythiaJetParticle *)jetcnd->getParticles()->At(ip);
      if(!chargedHadron(ppar)) continue;
      //cout<<ip<<" "<<ppar->pdg()<<endl;
      //TVector3 pparmom(ppar->px(), ppar->py(), ppar->pz());
      double pparz = parFrag(ppar, mom);
      mHist->hjetparz->Fill(pt, pparz, ww);
      TVector3 pparlocal = parLocalMomentum(ppar, mom);
      double pparjt = pparlocal.Perp();
      mHist->hjetparjt->Fill(pt, pparjt, ww);
      sumcharged += ppar->pt();
      multcharged++;
    }
    mHist->hjetht->Fill(pt, sumcharged/pt, ww);
    mHist->hjetparmult->Fill(pt, multcharged, ww);
    bool flag_match = false;
    StPythiaJetCandidate *parton = FindParJet(jetcnd, mPythiaJet[0]);
    if(parton){
      mHist->h_matched->Fill(pt, ww);
      mHist->h_matched_2d->Fill(pt, area, ww);
      mHist->hjetmatched->Fill(pt, ww);
      //double parton_pt = parton->pt();
      //double ptshift = parton_pt - pt;
      //mHist->hjetshift->Fill(pt, ptshift);
      flag_match = true;
    }/*else{
       cout<<"---eventId:"<<id<<"---"<<endl;
           PrintR(jetcnd, jetevnt[0]);
           cout<<"------"<<endl;
           }*/

    mHist->h_matched_prof->Fill(pt, flag_match, ww);
    double ue_pt = UnderlyingEventCone(jetcnd, flag_match, 0.5);//, ue_pt[0], ue_pt[1]);
    //ue_pt[2] = UnderlyingEvent_Cone(jetcnd, pythiaevnt, ue_pt[0], ue_pt[1]);
    double density = (ue_pt)/(2*3.1416*0.5*0.5);
    double ue_dpt = density * area;
    cout<<"jet: "<<pt<<" ue dpt: "<<ue_dpt<<" area: "<<area<<" ue sum dpt: "<<ue_pt<<endl;
    mHist->hjetmatchedratio->Fill(pt, flag_match, ww);
    mHist->hjetmatchedratiocrr->Fill(pt-ue_dpt, flag_match, ww);
    mHistResponse->hparticle->Fill(pt-ue_dpt, TMath::Abs(eta), ww);
    if(parton){
      double parton_pt = parton->pt();
      double parton_eta = parton->eta();
      double ptshift = parton_pt - pt;
      mHist->hjetshift->Fill(pt, ptshift, ww);
      cout<<"jet pt: "<<pt<<" matched found with pt: "<<parton_pt<<endl;
      cout<<"jet pt: "<<pt<<" pt shift: "<<ptshift<<endl;
      mHist->hjetshiftcrr->Fill(pt-ue_dpt, ptshift+ue_dpt, ww);
      //response
      mHistResponse->FillResponse(pt-ue_dpt, TMath::Abs(eta), parton_pt, TMath::Abs(parton_eta),ww);
    }
    //mHist->hjetuedpt->Fill(pt, ue_dpt);
    double rmin = FindRmin(jetcnd, mPythiaJet[0]);
    //         cout<<"rmin="<<rmin<<endl;
    mHist->h_rmin->Fill(rmin, ww);
    double pt_copy = pt;
    if(pt > 30.) pt_copy = 30.;
    mHist->h_rmin_pt->Fill(pt_copy, rmin, ww);
  }
  //cout<<"max: "<<ptmax<<" "<<dptmax<<" "<<flagmax<<endl;
  
  TClonesArray *array = mPythiaJet[0]->getJets();
  //      cout<<par[ipar]<<endl;
  int Njet = array->GetEntriesFast();
  //     cout<<"Njet:"<<Njet<<endl;
  for(int ijet = 0; ijet < Njet; ijet++){
    StPythiaJetCandidate *jetcnd = (StPythiaJetCandidate *)array->At(ijet);
    //        cout<<"jet pt:"<<jetcnd->pt()<<endl;
    double pt = jetcnd->pt();
    double eta = jetcnd->eta();
    //double area = jetcnd->area();
    //          cout<<"jet pt = "<<pt<<" area = "<<area<<endl;
    if(eta > -0.9 && eta < 0.9){
      mHistResponse->hparton->Fill(pt, TMath::Abs(eta), ww);
    }
  }
  
  return 1;
}
int readPythia::Finish()
{
  mHist->Write("Hist");
  mHistResponse->Write("HistResponse");
  mFile->Write();
  mFile->Close();
  return 1;
}
double readPythia::UnderlyingEventCone(const StPythiaJetCandidate *jet, bool flag_match, double Rcone)
{
  double pt_par = mPythia->pt();
  double ww = partonicWeight(pt_par);
  //  double DeltaR(double etaA, float phiA, float etaB, float phiB);
  double pt = jet->pt();
  double eta = jet->eta();
  double phi = jet->phi();
  double area = jet->area();
  TVector3 mom;
  mom.SetPtEtaPhi(pt, eta, phi);

  const double PI = 3.1416;
  //const int Ncone = 2;
  double cone_eta[Ncone];
  double cone_phi[Ncone];
  for(int icone = 0; icone < Ncone; icone++){
    cone_eta[icone] = eta;
    if(icone <= 1) cone_phi[icone] = phi + (PI/2.)*(2*icone-1);
    else cone_phi[icone] = phi;
  }

  double factor = area/(2*PI*TMath::Power(Rcone, 2));

  int Npar = mPythia->numberOfParticles();
  //int mstu72 = mPythia->mstu72();
  int mstu73 = mPythia->mstu73();
  //cout<<"mstu72: "<<mstu72<<" mstu73: "<<mstu73<<" Npar:"<<Npar<<endl;

  double ue_pt[Ncone];
  for(int icone = 0; icone < Ncone; icone++)
    ue_pt[icone] = 0.;
  TRefArray conepar[2];
  for(int ipar = mstu73; ipar < Npar; ipar++){
    TParticle *particle = const_cast<TParticle *>(mPythia->particle(ipar));
    //particles with status 1
    if(particle->GetStatusCode() != 1) continue;
    double p_eta = particle->Eta();
    double p_phi = particle->Phi();
    double p_pt = particle->Pt();

    double dR = DeltaR(p_eta, p_phi, cone_eta[0], cone_phi[0]);
    if(dR < Rcone){
      ue_pt[0] += p_pt;
      conepar[0].Add(particle);
    }else{
      dR = DeltaR(p_eta, p_phi, cone_eta[1], cone_phi[1]);
      if(dR < Rcone){
	ue_pt[1] += p_pt;
	conepar[1].Add(particle);
      }
    }
    //    for(int icone = 0; icone < Ncone; icone++){
    //      double dR = DeltaR(p_eta, p_phi, cone_eta[0], cone_phi[0]);
    //      if(dR < Rcone)
    //	      ue_pt[icone] += p_pt;
    //    }
  }

  double pt_sum = ue_pt[0] + ue_pt[1];
  ue_pt[2] = ue_pt[0] + ue_pt[1];
  cout<<"pt1 = "<<ue_pt[0]<<" pt2 = "<<ue_pt[1]<<" pt_sum = "<<pt_sum<<" mult1 = "<<conepar[0].GetEntriesFast()<<" mult2 = "<<conepar[1].GetEntries()<<endl;
  int mult[Ncone];
  mult[0] = conepar[0].GetEntriesFast();
  mult[1] = conepar[1].GetEntriesFast();
  mult[2] = mult[0] + mult[1];
  mHist->hjetuedpt->Fill(pt, pt_sum*factor, ww);
  mHist->hjetuemult->Fill(pt, mult[2], ww);
  double sumcharged = 0.;
  int multcharged = 0;
  for(int ic = 0; ic < 2; ic++){
     for(int ip = 0; ip < mult[ic]; ip++){
       TParticle *ppar = dynamic_cast<TParticle *>(conepar[ic].At(ip));
       if(!chargedHadron(ppar)) continue;
       double pparz = parFrag(ppar, mom);
       mHist->hjetueparz->Fill(pt, pparz, factor*ww);
       TVector3 mom_cone;
       mom_cone.SetPtEtaPhi(pt, cone_eta[ic], cone_phi[ic]);
       TVector3 pparlocal = parLocalMomentum(ppar, mom_cone);
       double pparjt = pparlocal.Perp();
       mHist->hjetueparjt->Fill(pt, pparjt, factor*ww);
       sumcharged += ppar->Pt()*factor;
       multcharged++;
     }
  }
  mHist->hjetueht->Fill(pt, sumcharged/pt,ww);
  mHist->hjetueparmult->Fill(pt, multcharged, ww);
  mHist->hjetuedpt_3d->Fill(pt, ue_pt[0], ue_pt[1], ww); 
  //matched with parton jet
  
  if(flag_match){
    for(int icone = 0; icone < Ncone; icone++){
      mHist->hue_m[icone]->Fill(ue_pt[icone], ww);
      mHist->hue_m_prof[icone]->Fill(pt, ue_pt[icone], ww);
    }
    mHist->hue_m_2d->Fill(ue_pt[0], ue_pt[1], ww);
    //mHist->hjetuedpt_m->Fill(pt, ue_dpt);
  }else{
    for(int icone = 0; icone < Ncone; icone++){
      mHist->hue_unm[icone]->Fill(ue_pt[icone], ww);
      mHist->hue_unm_prof[icone]->Fill(pt, ue_pt[icone], ww);
    }
    mHist->hue_unm_2d->Fill(ue_pt[0], ue_pt[1], ww);
  }

  return pt_sum;
}
