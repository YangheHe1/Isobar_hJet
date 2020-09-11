#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TRandom.h"
#include "TPythia6.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TSystem.h"

#include "ThrmRecPythia.h"
#include "ThrmFourVector.h"

#include "Riostream.h"

#include "Efficiency.h"

//---------------------------------------------
Double_t MakeCone(TPythia6 *pythia, TClonesArray *detearr, TClonesArray *pararr,Double_t pT, Double_t jetEta, Double_t jetPhi, Bool_t charged, Bool_t effi, Bool_t pTsmear, TString frag, bool fillEffHisto, TH1D* hTrackEff,TH1D* hparticle)
{
  //Double_t mom_smear=atof(gSystem->Getenv("MOMRES")); //TPC momentum resolution
  short mom_smear=atoi(gSystem->Getenv("MOMRES")); //TPC momentum resolution; 0: sigma=0.01*pT^2 (global tracks) 1: sigma=0.005*pT^2 (primary tracks, simple) 2: sigma=a+b*pT+c*pT^2 (primary tracks, more accurate)
  //Double_t smear_err=0; //uncertainty error on momentum resolution
  Double_t eff_increment=atof(gSystem->Getenv("EFF_INCREMENT")); //increase/decrese efficiency for systematic studies
   eff_increment=eff_increment/100.0; //[-100%,100%]=>[-1,1]
  Bool_t kCentral=atoi(gSystem->Getenv("CENTRAL"));//Central / peripheral collisions
  Int_t doScale=atoi(gSystem->Getenv("DOSCALE")); //rescale charged jet pT to the value of full jet
  TString ratio=gSystem->Getenv("EFFTYPE"); //hardon ratio (pi/K/p): pp-like | AuAu-like
//  short cutset=atoi(gSystem->Getenv("CUTSET"));//set of track cuts which we are using
	bool TOF_effi=atoi(gSystem->Getenv("TOFEFFI")); //correct for TOF+BEMC matching efficiency
  	Double_t tof_increment=0;
	if(TOF_effi)tof_increment=atof(gSystem->Getenv("TOFEFF_INCREMENT")); //increase/decrese TOF efficiency for systematic studies
	tof_increment=tof_increment/100.0; //[-100%,100%]=>[-1,1]

  Bool_t found = kFALSE;
  Double_t jetpTpart = 0;
  Double_t jetpTdete = 0;
  Double_t pTlead_dete = 0;
  Double_t scaler=1;		
  if(charged && doScale==2) pT= 3.0*pT/2.0; //since we take only charge part of the jet, jet pT would be ~2/3 of desired value

  TClonesArray *simarr = new TClonesArray("TParticle", 1000);

  //Double_t radius = 0.4;
  //TF1 *fJetProf = new TF1("fJetProf", "1 - 1/(1 + TMath::Exp(-(x-[0])/[1]))", 0, 1.0);
  //fJetProf->SetParameters(0.5, 0.05);
  
  while(!found)
    {
      simarr->Delete();
      detearr->Delete();

      //Double_t phi_parton = gRandom->Uniform()*2.0*TMath::Pi();
      Double_t phi_parton = jetPhi;
      //Double_t eta_parton = gRandom->Uniform(-1., +1.);
      Double_t eta_parton = jetEta;
      Double_t theta = 2.0*TMath::ATan(TMath::Exp(-1*eta_parton));
      //Double_t theta = 3.1415/2.0;
      Double_t E = pT/TMath::Sin(theta);
      
		//cout<<"pT:"<<pT<<"E:"<<E<<endl;

		//cout<<"E: "<<E<<" jet eta: "<<eta_parton<<" jet phi: "<<phi_parton<<endl;
		if(frag=="u")
   	   pythia->Py1ent(0, 2, E, theta, phi_parton); //u->jet
		else if (frag=="g")
	      pythia->Py1ent(0, 21, E, theta, phi_parton); //g->jet

      pythia->Pyexec();

      TLorentzVector partonlv(0, 0, 0, 0);
      partonlv.SetPtEtaPhiE(pT, eta_parton, phi_parton, E);

      Int_t final = pythia->ImportParticles(simarr, "Final");
      Int_t nparticles = simarr->GetEntries();
     //cout<<"npart: "<<nparticles<<endl; 
      TLorentzVector conelvPart; //particle level jet vector
      TLorentzVector conelvDete; //detector level jet vector
  		TLorentzVector partlvPart; //particle level jet constituent vector
		TLorentzVector partlvDete; //detector level jet constituent vector

      Int_t goodparpart = 0;
      Int_t goodparticle = 0;
	
      for(Int_t ipart = 0; ipart < nparticles; ipart++)
	{
	  TParticle *particle = (TParticle*)simarr->At(ipart);

	  
	  particle->Momentum(partlvPart);
	 //cout<<ipart<<"| eta:"<<partlvPart.Eta()<<" phi: "<<partlvPart.Phi()<<endl; 
	  partlvPart.SetPtEtaPhiM(partlvPart.Pt(), partlvPart.Eta(), partlvPart.Phi(), 0); //set M=0, probably not necessary

	  Double_t pTpart = partlvPart.Pt();
    
 	  if(pTpart < 0.2) continue;
	 
	 if(charged)
         {
          Double_t charge = particle->GetPDG()->Charge();
	  	  Int_t pid = particle->GetPDG()->PdgCode();
		  
          if(!charge) continue;
		  if(pid==2212||pid==-2212)hparticle->Fill(0);
          if(pid==211||pid==-211)hparticle->Fill(1);
          if(pid==321||pid==-321)hparticle->Fill(2);
	      if(abs(pid)!=2212 && abs(pid)!=211 && abs(pid)!=321)hparticle->Fill(3);
         }


	//"Gabriel's fragmentation"
	/*
     Double_t r = fJetProf->GetRandom()*radius;
	  Double_t alpha = TMath::TwoPi()*gRandom->Uniform();

	  Double_t eta = eta_parton + r*TMath::Sin(alpha);
	  Double_t phi = phi_parton + r*TMath::Cos(alpha);
	*/

	  Double_t eta = partlvPart.Eta();
	  Double_t phi = partlvPart.Phi();
	  Double_t M = 0;

	  new ((*pararr)[goodparpart]) ThrmFourVector(pTpart, eta, phi, M);
	  goodparpart++;
	  conelvPart += partlvPart;

	}//particle loop [in particle level jet]

     jetpTpart = conelvPart.Pt();
     //jetpTpart = conelvPart.P();
   if(!goodparpart) continue;
   found = kTRUE;

	if(charged && doScale==1) scaler = pT/jetpTpart; //rescale jet constituents' pT in order to get desired pT of the (particle level) charged jet
	else scaler=1;

	for(Int_t ipart = 0; ipart < pararr->GetEntries(); ipart++)
     {
        ThrmFourVector *fv = (ThrmFourVector*)pararr->At(ipart);
        TLorentzVector lv = fv->GetTLorentzVector();
		  Double_t pTpart=scaler*lv.Pt();
	     Double_t eta = lv.Eta();
   	  Double_t phi = lv.Phi();
     	  Double_t M = 0;
        fv->SetPtEtaPhiM(pTpart, eta, phi, M); //save the particle in (particle level jet) array again with scaled pT

			//tracking efficency
		   if(effi)
			{
			   //Double_t epsilon=efficiency11(pTpart, kCentral, ratio, cutset, eff_increment); //efficiency function 
			   Double_t epsilon=efficiencyAlex(pTpart, kCentral,ratio, eff_increment); //efficiency function 
			   //cout<<"pT: "<<pTpart<<" epsilon: "<<epsilon<<endl;
				if(fillEffHisto) hTrackEff->Fill(-pTpart); //we fill all particles to the negative side of the histogram and accepted particles to the positive 
	      	Double_t rand = gRandom->Uniform(0,1);
	      	if(rand>epsilon) continue;//drop the particle
				if(fillEffHisto) hTrackEff->Fill(pTpart);
	  	   }//apply tracking efficiency

			//TOF and/or BEMC matching efficency
			if(TOF_effi)
			{
				Double_t epsilon_TOF=efficiencyTOFvBEMC(pTpart,tof_increment);
				//cout<<"pT:"<<pTpart<<" eff TOF:"<<epsilon_TOF<<endl;
				Double_t rand_TOF = gRandom->Uniform(0,1);
				if(rand_TOF>epsilon_TOF) continue;
				//cout<<"PASSED"<<endl;
			}

		  //pTsmearing
	     if(pTsmear)
         {
            //Double_t sigma = (mom_smear+smear_err)*pTpart*pTpart;
            Double_t sigma = mom_res(pTpart,mom_smear);
				//cout<<"pT: "<<pTpart<<" sigma: "<<sigma<<"  sigma/pT^2: "<<sigma/(pTpart*pTpart)<<endl;
            Double_t ptn = gRandom->Gaus(pTpart,sigma);
            //cout<<"old "<<pTpart<<" new "<<ptn<<endl;
				if(ptn<0.2)continue;
            pTpart=ptn;
         }
 
        if(pTpart>pTlead_dete)pTlead_dete=pTpart;
		  partlvDete.SetPtEtaPhiM(pTpart, eta, phi, M);
		  conelvDete += partlvDete;
	  
		  new ((*detearr)[goodparticle]) ThrmFourVector(pTpart, eta, phi, M);
		  goodparticle++;

     }//particle loop [in particle and detector  level jet]

      jetpTdete = conelvDete.Pt();
	}//particle level jet exist 	

//if(pT>25)cout<<"jet pT at detector level:"<<jetpTdete<<endl;
//if(pT>25)cout<<"--------------------------------"<<endl;
  delete simarr;
  //delete fJetProf;
  //return jetpTdete;
  return jetpTpart;
}
//---------------------------------------------
Double_t MakeCone2(TPythia6 *pythia, TClonesArray *detearr, Double_t pT, Double_t jetEta, Double_t jetPhi, Bool_t charged, Bool_t effi, Bool_t pTsmear, TString frag,bool fillEffHisto, TH1D* hTrackEff)
{
	TClonesArray *pararr=new TClonesArray("ThrmFourVector", 1000);
	Double_t jetpTpart=MakeCone(pythia, detearr, pararr, pT, jetEta, jetPhi, charged, effi, pTsmear, frag, fillEffHisto, hTrackEff);
	pararr->Delete();
	delete pararr;
	return jetpTpart;
}

