#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TBranch.h"
#include "TRandom.h"
#include "TClonesArray.h"

#include "Rtypes.h"
#include "Riostream.h"

#include "ThrmJet.h"
#include "ThrmEmbedding.h"
#include "ThrmFourVector.h"

#include "JetReconstruction.h"

using namespace fastjet;
using namespace std;

//==========================================================================
JetReconstruction::JetReconstruction(TString wrkdir, Double_t rparam, Double_t pTcut) :
  fktjetsarr("ThrmJet", 100),
  faktjetsarr("ThrmJet", 100),
  fembeddingarr("ThrmEmbedding", 10)
{
  // FASTJET
  //akt_data_emb = new FJWrapper();

  frparam = rparam;
  frparam_BG = rparam; //we use same R for BG estimation as for jet reconstruction
  fpTcut = pTcut;
  fetaMax=1.0; //maximum pseudo rapidity acceptance for tracks
	TString outId=gSystem->Getenv("JOBID");

  TString str;
  
  // INPUT DATA
  str = Form("%s/toymodel_events_%s.root", wrkdir.Data(),outId.Data());
  finput = new TFile(str.Data(), "OPEN");
  ftreetoy = (TTree*)finput->Get("ToyModel");
  ftreetoy->SetBranchAddress("particles", &fpartarr);
  gROOT->cd();

  // OUTPUT DATA
  str = Form("%s/jets_R%.1lf_pTcut%.1lf_%s.root", wrkdir.Data(), frparam, fpTcut,outId.Data());
  foutput = new TFile(str.Data(), "RECREATE");
  ftreejets = new TTree("RecoJets","Reconstructed Jets");
  TBranch *br_embedding = ftreejets->Branch("embedding", &fembeddingarr);
  TBranch *br_ktjets = ftreejets->Branch("kt_jets", &fktjetsarr);
  TBranch *br_aktjets = ftreejets->Branch("akt_jets", &faktjetsarr);
  TBranch *br_rho = ftreejets->Branch("rho", frho);
  TBranch *br_sigma = ftreejets->Branch("sigma", fsigma);
}

//==========================================================================
JetReconstruction::~JetReconstruction()
{
  finput->Close();
  foutput->Close();

  delete finput;
  delete foutput;
}

//==========================================================================
void JetReconstruction::Run()
{
  Int_t nevts = ftreetoy->GetEntries();

  for(Int_t ievt = 0; ievt < nevts; ievt++)
    {
      
      ftreetoy->GetEntry(ievt);

      FillInputVector();

      RecoJets();

      RunEmbedding();
      
      ftreejets->Fill();

      Reset();

      if(ievt%100 == 0)
		cout << Form("Event #%d reconstructed", ievt) << endl;
    }

  foutput->Write();
}

//==========================================================================
void JetReconstruction::FillInputVector()
{
  Int_t nparticles = fpartarr->GetEntries();
  
  for(Int_t ipart = 0; ipart < nparticles; ipart++)
    {
      ThrmFourVector *partfv = (ThrmFourVector*)fpartarr->At(ipart);
      TLorentzVector partlv = partfv->GetTLorentzVector();

      if(partlv.Pt() < fpTcut) continue;
      
      PseudoJet inp_particle(partlv.Px(),
				      partlv.Py(),
				      partlv.Pz(),
				      partlv.Energy());
      inp_particle.set_user_index(ipart);
      input_data.push_back(inp_particle);
    }
}

//==========================================================================
void JetReconstruction::RecoJets()
{
	//FJ3
 	JetDefinition jet_def(antikt_algorithm, frparam);
 	// jet area definition
 	Double_t ghost_maxrap = 1.0; // Fiducial cut for background estimation
 	GhostedAreaSpec area_spec(ghost_maxrap);
 	//AreaDefinition area_def(active_area, area_spec);
 	AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));


 	ClusterSequenceArea clust_seq_hard(input_data, jet_def, area_def);
 	vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(fpTcut));
 	Selector Fiducial_cut_selector = SelectorAbsEtaMax(fetaMax - frparam); // Fiducial cut for jets
 	vector<PseudoJet> jets = Fiducial_cut_selector(jets_all);

   // background estimation
   JetDefinition jet_def_bkgd(kt_algorithm, frparam_BG);
   AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
	for(Int_t i = 0; i < 3; i++)
   {
   	Selector selector = SelectorAbsEtaMax(fetaMax) * (!SelectorNHardest(i));
   	JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
   	bkgd_estimator.set_particles(input_data);
   	frho[i]   = bkgd_estimator.rho();
   	fsigma[i] = bkgd_estimator.sigma();
   }


	// SAVING RECONSTRUCTED JETS
	Int_t Njets_par=jets.size();
	if(Njets_par<1);
	
		for(Int_t pjet=0; pjet<jets.size(); pjet++){
		
      	Double_t phi = jets[pjet].phi();
      	Double_t eta = jets[pjet].eta();
      	Double_t pT = jets[pjet].perp();
      	Double_t M = jets[pjet].m();
      	Double_t area = jets[pjet].area();
      	//Double_t area_par = akt_data.clust_seq->area(jets[pjet]);
      	//vector<PseudoJet> constituents = sorted_by_pt(akt_data.clust_seq->constituents(jets[pjet]));
			vector<PseudoJet> constituents=sorted_by_pt(clust_seq_hard.constituents(jets[pjet]));
      	Double_t pTleading = constituents[0].perp();
			int N=constituents.size();
	
			ThrmJet *jet;
	  		new ((faktjetsarr)[pjet]) ThrmJet();
	  		jet = (ThrmJet*)faktjetsarr.At(pjet);

      	ThrmFourVector jet_fv(pT, eta, phi, M);
      	jet->jet_fv = jet_fv;
      	jet->Nconst = N;
      	jet->area = area;
      	jet->pTleading = pTleading;
      	jet->embeddedJet_idx = constituents[0].user_index();

		}//loop over jets


//FJ2
/*
  akt_data->r = frparam;
  akt_data->maxrap = 1.;
  akt_data->algor = antikt_algorithm;
  akt_data->ghost_area = 0.01;
  akt_data->input_particles = input_data;
  akt_data->Run();

  kt_data->r = frparam;
  kt_data->maxrap = 1.;
  kt_data->algor = kt_algorithm;
  kt_data->ghost_area = 0.01;
  kt_data->input_particles = input_data;
  kt_data->Run();

  for(Int_t i = 0; i < 3; i++)
    {
      if((Int_t)kt_data->inclusive_jets.size() > i+1)
	kt_data->GetMedianAndSigma(frho[i], fsigma[i], i);
    }      
*/
}


//==========================================================================
void JetReconstruction::Reset()
{
  // RELEASING MEMORY
  //delete kt_data;
  //delete akt_data;

  input_data.clear();

  fembeddingarr.Delete();
  fktjetsarr.Delete();
  faktjetsarr.Delete();
}

//==========================================================================
void JetReconstruction::RunEmbedding()
{

	const Int_t nEmb=13;
	Float_t pTemb[] = {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0};
	for(Int_t iemb = 0; iemb < nEmb; iemb++)
	{
	      
      vector<PseudoJet> input_data_emb = input_data;
	          
      TLorentzVector emb_lv;
	
      Double_t eta = gRandom->Uniform((-fetaMax+frparam),(fetaMax-frparam));
      Double_t phi = gRandom->Uniform(-1.*TMath::Pi(), +1.*TMath::Pi());
      Double_t pT = pTemb[iemb];
/*
      // RUN PYTHIA
      // MAKE PARTICLE CONE
		
		if(kPythiaEmb){
      Double_t jetPt;
      Double_t jetEta = gRandom->Uniform((-1+frparam),(1-frparam));
      Double_t jetPhi = gRandom->Uniform(-1.*TMath::Pi(), +1.*TMath::Pi());
      if(kDeteCorr) jetPt = MakeCone(fpythia, simarr, pT, jetEta, jetPhi, kCharged, kEfficorr,kpTsmear,0, effL, effH);
      else jetPt = MakeCone(fpythia, simarr, pT, jetEta, jetPhi, kCharged, 0,kpTsmear);

      // SAVE RECO JET
      Int_t Nparticles = simarr->GetEntries();
      nhadrons += Nparticles;
      Float_t pTleading = 0;
      TLorentzVector total(0, 0, 0, 0);
      for(Int_t ipart = 0; ipart < Nparticles; ipart++)
   {
     ThrmFourVector *fv = (ThrmFourVector*)simarr->At(ipart);
     TLorentzVector lv = fv->GetTLorentzVector();

     Double_t pTpart = lv.Pt();
     Double_t eta = lv.Eta();
     Double_t phi = lv.Phi();
     Double_t M = lv.M();

     total += lv;

     if(pTpart > pTleading) pTleading = pTpart;

     fhpydNdpT->Fill(pTpart);

     new (fpartarr[goodparticle]) ThrmFourVector(pTpart, eta, phi, M);
     goodparticle++;
   }
}

*/
      emb_lv.SetPtEtaPhiM(pT, eta, phi, 0);
	
      PseudoJet emb_particle(emb_lv.Px(),
				      emb_lv.Py(),
				      emb_lv.Pz(),
				      emb_lv.Energy());
	
      emb_particle.set_user_index(99999);
      input_data_emb.push_back(emb_particle);


		//========================================
		//RUN JET RECONSTRUCTION
		//========================================
		//using Fastjet 3
	   JetDefinition jet_def(antikt_algorithm, frparam);
  		// jet area definition
   	Double_t ghost_maxrap = 1.0; // Fiducial cut for background estimation
   	GhostedAreaSpec area_spec(ghost_maxrap);
   	//AreaDefinition area_def(active_area, area_spec);
   	AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));

   	ClusterSequenceArea clust_seq_hard(input_data_emb, jet_def, area_def);
   	vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(fpTcut));
   	Selector Fiducial_cut_selector = SelectorAbsEtaMax(fetaMax - frparam); // Fiducial cut for jets
   	vector<PseudoJet> embedding_jets = Fiducial_cut_selector(jets_all);
		
		//========================================
		//find embedded jet
		//========================================
  		Bool_t jetEmbedded = kFALSE;
		double pTreco,Area,pTleading;
	
		for(Int_t ijet = 0; ijet <embedding_jets.size(); ijet++)
   	{
      	Double_t eta_jet = embedding_jets[ijet].eta();
      	if(TMath::Abs(eta_jet) > fetaMax - frparam) continue;

			//let's look at the jet constituents	      
      	vector<PseudoJet> constituents = sorted_by_pt(clust_seq_hard.constituents(embedding_jets[ijet]));
			int N=constituents.size();
      	for(Int_t iConst = 0; iConst < N; iConst++)
				if(constituents[iConst].user_index() == 99999)
	  			{
	    			pTreco = embedding_jets[ijet].perp();
	    			Area = embedding_jets[ijet].area();
	    			jetEmbedded = kTRUE;
	    			pTleading = constituents[0].perp();
	  			}
      		if(jetEmbedded) break;
  	 	}    

		Double_t etaEmb = eta;
		Double_t phiEmb = phi;
   	Double_t PtEmb = pT;
   	Double_t PtEmbReco = pTreco;
   	Double_t AreaEmb = Area;
   	Double_t deltapT = pTreco - frho[0]*Area - pT;

   	new ((fembeddingarr)[iemb]) ThrmEmbedding(jetEmbedded, PtEmb, etaEmb, phiEmb, AreaEmb, deltapT, PtEmbReco, pTleading);
   }//pTembed
}

//==========================================================================
