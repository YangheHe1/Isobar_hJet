#define StMiniMcTree_cxx
#include "StMiniMcTree.h"
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
/*
  #include "StRoot/StarRoot/TRVector.h"
  #include "StRoot/StarRoot/TRArray.h"
  #include "StRoot/StarRoot/TRMatrix.h"
  #include "StRoot/StarRoot/TRSymMatrix.h"
  #include "StRoot/StarRoot/TRDiagMatrix.h"
  #include "StRoot/StarRoot/TRSymMatrix.h"
  #include "TRVector.h"
  #include "TRMatrix.h"
  #include "TRSymMatrix.h"
*/
#include "TROOT.h"
#include "Riostream.h"
#include <stdio.h>
#include "TSystem.h"
#include "TMath.h"
#include "cuts.h"

void StMiniMcTree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L StMiniMcTree.C
//      Root > StMiniMcTree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


	int doCentral=atoi(gSystem->Getenv("CENTRAL")); //central or peripheral collisions
   TString particle=gSystem->Getenv("PARTICLE"); //pi,K,p
	bool doGlobal=atoi(gSystem->Getenv("GLOBAL")); //global or primary tracks
	int geantID1=-1, geantID2=-1; //geantID of the embedded particle (and antiparticle)
	int outname= atoi(gSystem->Getenv("OUTPUT"));
	int fNjobs = atoi(gSystem->Getenv("NJOBS"));
	if(particle=="P")
	{
		geantID1=14;
		geantID2=15;
	}
	else if(particle=="Pi")
	{
		geantID1=8;
		geantID2=9;
	}
	else if(particle=="K")
	{
		geantID1=11;
		geantID2=12;
	}

   if (fChain == 0) return;
   
   //TFile *myfile = new TFile("results_minimctree_all.root","RECREATE");
   TFile *myfile = new TFile(Form("embeding/%s/MC_%s_cent%i_%i.root",particle.Data(),particle.Data(),doCentral,outname),"RECREATE");
   myfile->cd();
   
	int npTbins=200;
	double pTmax=20;
	TH2D* hcent_mult=new TH2D("hcent_mult", "centrality bin vs refMult", 11, -0.5 , 10.5, 150,0, 600); 
	TH2D* hpTRec_pTMc=new TH2D("hpTRec_pTMc", "pT reco vs pT MC", npTbins, 0, pTmax, npTbins, 0,pTmax); //for pt smearing calculation

	TH1D* hpTMC_MC=new TH1D("hpTMC_MC", "pTMC of MC tracks", npTbins, 0, pTmax); //for efficiency calculation
	TH1D* hpTMC_Rec=new TH1D("hpTMC_Rec", "pTMC of matched tracks", npTbins, 0, pTmax); //for efficiency calculation

	TH3D* h3d_MC= new TH3D("h3d_MC","pt,eta,phi for MC tracks",npTbins, 0, pTmax,20,-1,1,80,-4,4); //for 3d efficiency
	TH3D* h3d_Rec= new TH3D("h3d_Rec","pt,eta,phi for Rec tracks",npTbins, 0, pTmax,20,-1,1,80,-4,4); //for 3d efficiency

	TH1D* heta_MC = new TH1D("heta_MC","eta for MC tracks",20,-1,1);
	TH1D* heta_Rec = new TH1D("heta_Rec","eta for Rec tracks",20,-1,1);
	TH1D* hphi_MC = new TH1D("hphi_MC","phi for MC tracks",80,-4,4);
	TH1D* hphi_Rec = new TH1D("hphi_Rec","phi for Rec tracks",80,-4,4);
	TH1I* hgeant=new TH1I("hgeant","geant id",1000,0.5,1000.5);
	TH1D* hparticle = new TH1D("hparticle","0P 1Pi 2K",4,0,4);

    fChain->SetBranchStatus("*",0);  // disable all branches
    fChain->SetBranchStatus("mVertexZ",1);  // activate branchname
    fChain->SetBranchStatus("mCentrality",1);  // activate branchname
    fChain->SetBranchStatus("mCentralMult",1);  // activate branchname
    fChain->SetBranchStatus("mNMatchedPair",1);  // activate branchname
    fChain->SetBranchStatus("mMatchedPairs.*",1);  // activate branchname
    fChain->SetBranchStatus("mNMcTrack",1);  // activate branchname
    fChain->SetBranchStatus("mMcTracks.*",1);  // activate branchname
    fChain->SetBranchStatus("mNMatGlobPair",1);  // activate branchname
    fChain->SetBranchStatus("mMatGlobPairs.*",1);  // activate branchname

   Long64_t nentries = fChain->GetEntries();
   cout << "nentries = " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
		int enm=fChain->GetEntryNumber(jentry);
		//if(enm%100==0)cout<<"entry #:"<<enm<<endl;
      nb=fChain->GetEntry(jentry);  
		nbytes += nb;

		if(mVertexZ>zcut)continue;
   
		if(mCentralMult<=refMultCutMin[doCentral] || mCentralMult>refMultCutMax[doCentral]) continue;
		//cout<<"centrality:"<<mCentrality<<" multiplicity:"<<mCentralMult<<endl;  
		hcent_mult->Fill(mCentrality,mCentralMult);

		Double_t PMC_trk1=0,PtMC_trk1=0,PR_trk1=0,PtPR_trk1=0,EtaPr_trk1=0,PhiPr_trk1=0,dcaXYGl_trk1=0,dcaGl_trk1=0,dcaXYGlMc_trk1=0,Chi2Pr_trk1=0;
      Int_t TPC_trk1=0,geant_trk1=0,parentGeant_trk1=0,parentKey_trk1,Q_trk1=0, key_trk1=0, Recokey_trk1=0;

      if(!doGlobal && mNMatchedPair<1) continue;
      if(doGlobal && mNMatGlobPair<1) continue;

      if(jentry%1000 == 0){   cout << " jentry = " << jentry << " eventId = " <<  mEventId << "# of mMatGlobPairs = " << mNMatGlobPair << endl;}
		for(Int_t trk =0; trk < mNMcTrack; trk++) //loop over MC tracks
		{
			//track cuts
			if(TMath::Abs(mMcTracks_mEtaMc[trk])>EtaMax)continue;
			if(mMcTracks_mNHitMc[trk] <TpcMcCut) continue;
			int geant_mctrk  = mMcTracks_mGeantId[trk];

			if(geant_mctrk==14||geant_mctrk==15) hparticle->Fill(0);
			if(geant_mctrk==8||geant_mctrk==9) hparticle->Fill(1);
			if(geant_mctrk==11||geant_mctrk==12) hparticle->Fill(2);


			if((geantID1>0 && geant_mctrk!=geantID1) && (geantID2>0 && geant_mctrk!=geantID2)) continue; //geant id doesn't match the embedded particle/antiparticle geant id
			if(mMcTracks_mPtMc[trk]<pTmin || mMcTracks_mPtMc[trk]>pTmax) continue;

			//fill histograms
			hpTMC_MC->Fill(mMcTracks_mPtMc[trk]);
			heta_MC->Fill(mMcTracks_mEtaMc[trk]);
			hphi_MC->Fill(mMcTracks_mPhiMc[trk]);
			h3d_MC->Fill(mMcTracks_mPtMc[trk],mMcTracks_mEtaMc[trk],mMcTracks_mPhiMc[trk]);
		}

		if(!doGlobal) //primary tracks
		{
	      for(Int_t trk =0; trk < mNMatchedPair; trk++)
			{ 
			  //track cuts
			  //if(TMath::Abs(mMatchedPairs_mEtaMc[trk]) > EtaMax) continue;
			  if(TMath::Abs(mMatchedPairs_mEtaPr[trk]) > EtaMax) continue;
			  if(mMatchedPairs_mPtPr[trk]<pTmin || mMatchedPairs_mPtPr[trk]>pTmax) continue;
			  if(mMatchedPairs_mFitPts[trk] < TpcCut) continue;
			  if(((float)mMatchedPairs_mFitPts[trk]/(float)mMatchedPairs_mNPossible[trk]) < RatioHit) continue;
			  if(mMatchedPairs_mFlag[trk]>FlagMax || mMatchedPairs_mFlag[trk]<=0) continue;

			  dcaGl_trk1      = mMatchedPairs_mDcaGl[trk];
			  Chi2Pr_trk1 = mMatchedPairs_mChi2Pr[trk];
			  if(dcaGl_trk1>DcaCutMax) continue;
			  if(Chi2Pr_trk1>Chi2Cut) continue;

			  geant_trk1       = mMatchedPairs_mGeantId[trk];
			  if((geantID1>0 && geant_trk1!=geantID1) && (geantID2>0 && geant_trk1!=geantID2)) continue; //geant id doesn't match the embedded particle/antiparticle geant id
				
			  PtPR_trk1        = mMatchedPairs_mPtPr[trk];
			  PtMC_trk1        = mMatchedPairs_mPtMc[trk];
			  
			  //fill histograms
			  hpTRec_pTMc->Fill(PtPR_trk1,PtMC_trk1); 
			  hpTMC_Rec->Fill(PtMC_trk1);
			  heta_Rec->Fill(mMatchedPairs_mEtaMc[trk]);
			  hphi_Rec->Fill(mMatchedPairs_mPhiMc[trk]);
			  h3d_Rec->Fill(PtMC_trk1,mMatchedPairs_mEtaMc[trk],mMatchedPairs_mPhiMc[trk]);
			  hgeant->Fill(geant_trk1);

			}//track loop
		}//primary tracks
	else //global tracks
	{

		//*******************************************

      for(Int_t trk =0; trk < mNMatGlobPair; trk++)
		{ 
	
		   if(mMatGlobPairs_mPtGl[trk]<pTmin || mMatGlobPairs_mPtGl[trk]>pTmax) continue;
		   if(mMatGlobPairs_mFitPts[trk] < TpcCut) continue;
		   if(((float)mMatGlobPairs_mFitPts[trk]/(float)mMatGlobPairs_mNPossible[trk]) < RatioHit) continue;
		   if(TMath::Abs(mMatGlobPairs_mEtaPr[trk]) > EtaMax) continue;
		   //if(mMatGlobPairs_mFlag[trk]>FlagMax || mMatGlobPairs_mFlag[trk]<=0) continue;

		   dcaGl_trk1      = mMatGlobPairs_mDcaGl[trk];
		   //Chi2Pr_trk1 = mMatGlobPairs_mChi2Pr[trk];
		   //if(dcaGl_trk1>DcaCutMax) continue;
		   //if(Chi2Pr_trk1>Chi2Cut) continue;

	      geant_trk1       = mMatGlobPairs_mGeantId[trk];
		   if((geantID1>0 && geant_trk1!=geantID1) && (geantID2>0 && geant_trk1!=geantID2)) continue; //geant id doesn't match the embedded particle/antiparticle geant id

		   PtPR_trk1        = mMatGlobPairs_mPtGl[trk];
		   PtMC_trk1        = mMatGlobPairs_mPtMc[trk];

	  		//fill histograms
	  		hpTRec_pTMc->Fill(PtPR_trk1,PtMC_trk1); 
			hpTMC_Rec->Fill(PtMC_trk1);
			heta_Rec->Fill(mMatchedPairs_mEtaMc[trk]);
			hphi_Rec->Fill(mMatchedPairs_mPhiMc[trk]);
			h3d_Rec->Fill(PtMC_trk1,mMatchedPairs_mEtaMc[trk],mMatchedPairs_mPhiMc[trk]);
	  		hgeant->Fill(geant_trk1);
		}//track loop
	}//global tracks
	}//event loop

	for(int j=1; j<npTbins+1; j++){
		double ntot_rm=hpTMC_Rec->GetBinContent(j);
	//loop over detector level jet pT bins
	for(int i=1; i<npTbins+1; i++){
		double binc = hpTRec_pTMc->GetBinContent(i,j);
		double rat;
		if(ntot_rm>0){
		  rat=binc/ntot_rm;
		  hpTRec_pTMc->SetBinContent(i,j,rat);
		}

	}//i - detector level
	}//j - particle level

float scaler=(float)1.0/fNjobs;
hpTRec_pTMc->Scale(scaler);

myfile->Write();
hpTMC_Rec->Divide(hpTMC_MC);
hpTMC_Rec->Scale(scaler);
hpTMC_Rec->Write("heffi");

heta_Rec->Divide(heta_MC);
heta_Rec->Scale(scaler);
heta_Rec->Write("heta_effi");

hphi_Rec->Divide(hphi_MC);
hphi_Rec->Scale(scaler);
hphi_Rec->Write("hphi_effi");

h3d_Rec->Divide(h3d_MC);
h3d_Rec->Scale(scaler);
h3d_Rec->Write("h3d_effi");


myfile->Close();

}
