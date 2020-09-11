#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TDatime.h"
#include "TBranch.h"
#include "TRandom.h"
#include "TSystem.h"
#include "Riostream.h"
#include "TClonesArray.h"
#include "dNdpT.h"
#include "ThrmFourVector.h"
#include "ThrmRecPythia.h"
#include "ThrmMemStat.h"
#include "ThrmEmbedding.h"
#include "ThrmPythiaEmb.h"

//FJ3
#include <fastjet/config.h>            
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/ClusterSequenceActiveArea.hh>
#include <fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/Subtractor.hh> 
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>

using namespace fastjet;
using namespace std;

//=================================================

ThrmPythiaEmb::ThrmPythiaEmb()
:fembeddingarr("ThrmEmbedding", 20)
{

//load parameters----------------------------------
fR = atof(gSystem->Getenv("RPARAM"));
fmax_rap = atof(gSystem->Getenv("MAXRAP"));
fAreaCut = atof(gSystem->Getenv("AREACUT"));
kCharged = atoi(gSystem->Getenv("CHARGED"));
fNjobs = atoi(gSystem->Getenv("NJOBS"));
kEfficorr=atoi(gSystem->Getenv("EFFICORR")); //apply efficiency correction cuts (for detector level jets)
kpTsmear=atoi(gSystem->Getenv("PTSMEAR")); //apply track pT smearing (for detector level)
fpTcut=0.2;
frag = gSystem->Getenv("JETFRAG"); //jet fragmentation: u | g
//TString outDir = gSystem->Getenv("OUTPUTDIR");
TString outDir = gSystem->Getenv("SCRATCH");
TString outId=gSystem->Getenv("JOBID");
kFixedSeed = atoi(gSystem->Getenv("FIXEDSEED")); //1: use fixed number (fSeed) as a seed for RN generator 0: use computer time as a seed
kSeed = atoi(gSystem->Getenv("SEED")); //Seed for random number generatr
//TString eff_path = gSystem->Getenv("EFFPATH"); //path to the single hadron efficiency functions
//TString eff_type = gSystem->Getenv("EFFTYPE"); //single hadron efficiency: AuAu | pp - like scenario

pTmin=0.5;
pTmax=100;

//OUTPUT
TString str = Form("%s/pythia_emb_R%.1lf_%s.root", outDir.Data(), fR, outId.Data());
foutput = new TFile(str.Data(), "RECREATE");
cout<<"opening output file "<<str.Data()<<endl;
//ftreejets = new TTree("PythiaEmbJets","Pythia Jets");
//TBranch *br_pythia = ftreejets->Branch("pythia", &fembeddingarr);

//set acceptance
etaMinCut=-(fmax_rap-fR);
etaMaxCut=(fmax_rap-fR);
//etaMinCut=-1.;
//etaMaxCut=+1.;


//**************************
//histogram definitions
//**************************
nptbins=800;
Float_t ptminbin=-pTmax;
Float_t ptmaxbin=pTmax;

//tracking efficiency (just to be sure it looks OK)
hTrackEff=new TH1D("hTrackEff","histo for calculation of tracking efficiency",120,-30,30);
hTrackEpsilon=new TH1D("hTrackEpsilon","tracking efficiency",60,0,30);

hparticle = new TH1D("hparticle","0P 1Pi 2K",4,0,4);

for(Int_t pTl=0; pTl<npTlead; pTl++)
{
  TString name="";
  //response matrix
  name=Form("hresponse_pTl%i",pTl);
  hresponse[pTl]=new TH2D(name, "ResponseMatrix;p_{T,jet}^{det,ch}(GeV/c);p_{T,jet}^{part,ch}(GeV/c);entries", nptbins, ptminbin, ptmaxbin, nptbins, ptminbin, ptmaxbin);
  name=Form("hreponse_Y_pTl%i",pTl);
  hRM_projectionY[pTl]=new TH1D(name, "RM projection to Y axis",nptbins, ptminbin, ptmaxbin);
  //jet reconstruction efficiency histogram
  name=Form("hJetEff1D_pTl%i",pTl);
  hJetEff1D[pTl]=new TH1D(name, "jet reco. efficiency after detector corrections",nptbins, ptminbin, ptmaxbin);
  name=Form("hJetEff2D_pTl%i",pTl);
  hJetEff2D[pTl]=new TH2D(name, "hist. for calc. of jet reco. efficiency after dete. corr.",nptbins, ptminbin, ptmaxbin, nptbins, ptminbin, ptmaxbin);
  //old efficiency histograms
  //name=Form("heffi_par_pTl%i",pTl);
  //heffi_par[pTl]=new TH1D(name, "jet reco. efficiency after detector corrections",nptbins/2, 0, ptmaxbin); //not used anymore
  //name=Form("heffi_dete_pTl%i",pTl);
  //heffi_dete[pTl]=new TH1D(name, "jet reco. efficiency after detector corrections",nptbins/2, 0, ptmaxbin); //not used anymore
  //histograms for saving the number of generated jets
  name=Form("hntrue_all_pTl%i",pTl);
  hntrue_all[pTl]=new TH1D(name, "n true generated (for efficiency);p_{T,jet}^{part}(GeV/c);entries", nptbins, ptminbin, ptmaxbin);
  name=Form("hntrue_rm_pTl%i",pTl);
  hntrue_rm[pTl]=new TH1D(name, "n true generated (for RM);p_{T,jet}^{part}(GeV/c);entries", nptbins, ptminbin, ptmaxbin);
  name=Form("hndete_rm_pTl%i",pTl);
  hndete_rm[pTl]=new TH1D(name, "n dete generated (for RM);p_{T,jet}^{det}(GeV/c);entries", nptbins,0, 200);
name=Form("hnpart_rm_pTl%i",pTl);
  hnpart_rm[pTl]=new TH1D(name, "n part generated (for RM);p_{T,jet}^{part}(GeV/c);entries", nptbins, 0,200);

  //Jet Energy Resolution histogram
  name=Form("hJER_pTl%i",pTl);
  hJER[pTl]=new TH2D(name, "Jet Energy Resolution;p_{T}^{part};(p_{T}^{dete}-p_{T}^{part})/p_{T}^{part};entries", 50, 0, 50, 80,-0.5,0.3);

  //jet pt vs eta
  name=Form("hntruerm_vs_eta_pTl%i",pTl);
  hntruerm_vs_eta[pTl]=new TH2D(name,"true jet pt(for RM) vs eta;#eta;p_{T,jet}^{part}(GeV/c)",40,-2,2,nptbins, ptminbin, ptmaxbin);
  name=Form("hndeterm_vs_eta_pTl%i",pTl);
  hndeterm_vs_eta[pTl]=new TH2D(name,"dete jet pt(for RM) vs eta;#eta;p_{T,jet}^{det}(GeV/c)",40,-2,2,nptbins, ptminbin, ptmaxbin);

}

// SETTING RANDOM SEED
  TDatime dt;
  UInt_t curtime = dt.Get();
  UInt_t procid = gSystem->GetPid();
  UInt_t seed;
  if(kFixedSeed)
	seed= kSeed + procid;
  else
	seed = curtime - procid;
  gRandom->SetSeed(seed);

// PYTHIA
fpythia = new TPythia6();
fpythia->SetMRPY(1, seed);



}

//=================================================
ThrmPythiaEmb::~ThrmPythiaEmb()
{
delete effL;
delete effH;
delete fpythia; 

foutput->Close();
delete ftreejets;
delete foutput;
}

//=================================================
double ThrmPythiaEmb::getDistance(Double_t eta1, Double_t eta2,Double_t phi1, Double_t phi2)
{
		Double_t delta_eta=TMath::Abs(eta1-eta2);
		Double_t delta_phi=TMath::Abs(phi1-phi2);
		if(delta_phi>TMath::Pi()) delta_phi=2*TMath::Pi()-delta_phi;
		double dist=TMath::Sqrt(delta_eta*delta_eta+delta_phi*delta_phi);
		return dist;
}

//=================================================
void ThrmPythiaEmb::Run(Int_t nevents)
{

//const Int_t nEmb=18;
//Float_t pTemb[] = {0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0, 30.0, 50.0, 75.0, 100.0};

TClonesArray *detearr = new TClonesArray("ThrmFourVector", 1000);
TClonesArray *pararr = new TClonesArray("ThrmFourVector", 1000);

memstat.Start();

//EVENT LOOP 
cout<<"Run set for "<<nevents<<" events"<<endl;
for(Int_t ievt = 0; ievt < nevents; ievt++)
{
if(ievt%1000==0){
cout<<"filling event "<<ievt<<endl;
cout<<"Memory used: "<<memstat.Used()<<endl;
}
  //for(Int_t iemb = 0; iemb < nEmb-1; iemb++)
	//{
	detearr->Delete();
	pararr->Delete();

	Double_t jetEta = gRandom->Uniform(-fmax_rap, +fmax_rap);
	Double_t jetPhi = gRandom->Uniform(-1., +1.)*TMath::Pi();
	Double_t pT = gRandom->Uniform(pTmin,pTmax); 

	//randomly select u-quark or gluon fragmentation if we want a mixture
	TString usefrag=frag;
	if (frag=="2u1g")
	{
		float randn=gRandom->Uniform(0,3);
		if(randn<2) usefrag="u";
		else usefrag="g";
	}
	
	//create pythia jet
	Double_t PtJet=MakeCone(fpythia, detearr, pararr, pT, jetEta, jetPhi, kCharged, kEfficorr,kpTsmear,usefrag, true, hTrackEff,hparticle); 

//cout<<"req: "<<pT<<" gen: "<<PtJet<<endl;

	//loop over particles in particle-jet
	Int_t Nparticles = pararr->GetEntries();
   vector<PseudoJet> input_vector;
	//Float_t pTleading = 0;
	//TLorentzVector total(0, 0, 0, 0);
	for(Int_t ipart = 0; ipart < Nparticles; ipart++)
  	{
     ThrmFourVector *fv = (ThrmFourVector*)pararr->At(ipart);
     TLorentzVector lv = fv->GetTLorentzVector();
     
     Double_t pTpart = lv.Pt();
     //total += lv;
     //if(pTpart > pTleading) pTleading = pTpart;

		if(lv.Pt() < fpTcut) continue;
      if(TMath::Abs(lv.Eta()) > fmax_rap) continue;
      PseudoJet inp_particle(lv.Px(),lv.Py(),lv.Pz(),lv.Energy()); 
		input_vector.push_back(inp_particle);
	}

	//FASTJET
/*
//FJ2
	FJWrapper akt_data;
   akt_data.r = fR;
   akt_data.maxrap = fmax_rap;
   akt_data.algor = antikt_algorithm;
   akt_data.input_particles = input_vector;
   akt_data.Run();

   vector<PseudoJet> jets = sorted_by_pt(akt_data.inclusive_jets);
//FJ2
*/

//FJ3
 JetDefinition jet_def(antikt_algorithm, fR);
 // jet area definition
 Double_t ghost_maxrap = 1.0; // Fiducial cut for background estimation
 GhostedAreaSpec area_spec(ghost_maxrap);
 //AreaDefinition area_def(active_area, area_spec);
 AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));

 ClusterSequenceArea clust_seq_hard(input_vector, jet_def, area_def);
 vector<PseudoJet> jets_all = sorted_by_pt(clust_seq_hard.inclusive_jets(fpTcut));
 Selector Fiducial_cut_selector = SelectorAbsEtaMax(1.0 - fR); // Fiducial cut for jets
 vector<PseudoJet> jets = Fiducial_cut_selector(jets_all);
//FJ3

	Int_t Njets_par=jets.size();
	if(Njets_par<1)continue;

	//[A] PARTICLE LEVEL: one-jet events only
	/*
   if(Njets_par>1)continue;
   Int_t pjet=0;
	*/

	//[B] PARTICLE LEVEL: random jet
	/*
	Double_t rnd= gRandom->Uniform(1,Njets_par);
	Int_t pjet=(Int_t) rnd;
	*/

	//[C] PARTICLE LEVEL: highest pT jet
	//Int_t pjet=0;

	//[D] PARTICLE LEVEL: all jets
		
	//cout<<"rnd:"<<rnd<<" pjet: "<<pjet<<" njets: "<<Njets_par<<endl;
	//cout<<" pjet: "<<pjet<<" njets: "<<Njets_par<<endl;

	for(Int_t pjet=0; pjet<jets.size(); pjet++){
	
      Double_t phi_par = jets[pjet].phi();
      Double_t eta_par = jets[pjet].eta();
      Double_t pT_par = jets[pjet].perp();
      Double_t M_par = jets[pjet].m();
      Double_t area_par = jets[pjet].area();
      //Double_t area_par = akt_data.clust_seq->area(jets[pjet]);
      //vector<PseudoJet> constituents = sorted_by_pt(akt_data.clust_seq->constituents(jets[pjet]));
		vector<PseudoJet> constituents = sorted_by_pt(clust_seq_hard.constituents(jets[pjet]));
      Double_t pTleading_par = constituents[0].perp();

		if(eta_par<etaMinCut || eta_par>etaMaxCut) continue; // fiducial acceptance 

	   //FILLING EFFICIENCY HISTOGRAMS
		for(Int_t pTl=0; pTl<npTlead; pTl++){
			if(pTleading_par>pTl /*&& area_par>fAreaCut*/)
			{
				//heffi_par[pTl]->Fill(pT_par); //not used anymore
				hntrue_all[pTl]->Fill(pT_par);
			}
		}//loop over pTlead cuts


	//loop over particles in detector-jet 
	Nparticles = detearr->GetEntries();
	input_vector.clear();
	//pTleading = 0;
	//total.SetXYZM(0,0,0,0);
	for(Int_t ipart = 0; ipart < Nparticles; ipart++)
  	{
     ThrmFourVector *fv = (ThrmFourVector*)detearr->At(ipart);
     TLorentzVector lv = fv->GetTLorentzVector();
     //Double_t pTpart = lv.Pt();
     //total += lv;
     //if(pTpart > pTleading) pTleading = pTpart;
		if(lv.Pt() < fpTcut) continue;
      if(TMath::Abs(lv.Eta()) > fmax_rap) continue;
      PseudoJet inp_particle(lv.Px(),lv.Py(),lv.Pz(),lv.Energy()); 
		input_vector.push_back(inp_particle);
		//delete fv;
		}

	/*
   Double_t eta_dete = total.Eta(); 
	Double_t phi_dete = total.Phi();
	Double_t pT_dete = total.Pt();
	Double_t pTleading_dete = pTleading;
	Double_t deltapT=pT_dete - pT_par;
   */

	//FASTJET
/*
//FJ2
  	FJWrapper akt_data2;
   akt_data2.r = fR;
   akt_data2.maxrap = fmax_rap;
   akt_data2.algor = antikt_algorithm;
   akt_data2.input_particles = input_vector;
	akt_data2.Run();

   //jets.clear();
   vector<PseudoJet> jets2 = sorted_by_pt(akt_data2.inclusive_jets);
//FJ2
*/
//FJ3
 ClusterSequenceArea clust_seq_hard2(input_vector, jet_def, area_def);
 vector<PseudoJet> jets_all2 = sorted_by_pt(clust_seq_hard2.inclusive_jets(fpTcut));
 vector<PseudoJet> jets2 = Fiducial_cut_selector(jets_all2);
//FJ3 

	Int_t Njets_dete=jets2.size();
	if(jets2.size()<1)continue;

   //variables for jet matching	
	Float_t mindist1=100;
   bool matched1=false;
	Int_t closest_dete=-1;
	Int_t closest_par=-1;
	
   Double_t phi_dete0;
   Double_t eta_dete0;
   Double_t pT_dete0=-1;
   Double_t M_dete0;
   Double_t area_dete0;
	Double_t pTleading_dete0;

	//DETECTOR LEVEL JET LOOP
	for(Int_t djet=0; djet<jets2.size(); djet++){
      Double_t phi_dete = jets2[djet].phi();
      Double_t eta_dete = jets2[djet].eta();
      Double_t pT_dete = jets2[djet].perp();
      Double_t M_dete = jets2[djet].m();
      Double_t area_dete = jets2[djet].area();
      //Double_t area_dete = akt_data2.clust_seq->area(jets2[djet]);

		if((eta_dete<etaMinCut || eta_dete>etaMaxCut) || area_dete<fAreaCut) continue;// fiducial acceptance + area cut
     //Int_t N = akt_data2.clust_seq->constituents(jets2[0]).size();

      //vector<PseudoJet> constituents = sorted_by_pt(akt_data2.clust_seq->constituents(jets2[djet]));
		vector<PseudoJet> constituents2 = sorted_by_pt(clust_seq_hard2.constituents(jets2[djet]));
	   Double_t pTleading_dete = constituents2[0].perp();

		//FILLING EFFICIENCY HISTOGRAMS - for all jets
		if(pjet==0){	//make sure we fill the histogram for each detector level jet only once
		for(Int_t pTl=0; pTl<npTlead; pTl++){
		   if(pTleading_dete>pTl)	
			{
				//not used anymore
				//heffi_dete[pTl]->Fill(pT_dete); //for efficiency callculation - fill pT of detected detector level jet
			}
		}//pTleading loop
		}

  	   //JET MATCHING I - find closest detector level jet && largest energy fraction 
		//Double_t jet_energy_fraction=pT_dete/pT_par;
		Double_t dist=getDistance(eta_par,eta_dete,phi_par,phi_dete);
		if(dist<mindist1 && dist<fR){
		mindist1=dist;
		matched1=true;	
		closest_dete=djet;
		phi_dete0=phi_dete;
		eta_dete0=eta_dete;
		pT_dete0=pT_dete;
		M_dete0=M_dete;
		area_dete0=area_dete;
		pTleading_dete0=pTleading_dete;
		}
	}//detector jets loop

	if(!matched1)continue;

	//JET MATCHING II - find closest particle level jet
	Float_t mindist2=100;
	bool matched2=false;
	for(Int_t ijet=0; ijet<jets.size(); ijet++){
     	Double_t phi_p = jets[ijet].phi();
	   Double_t eta_p = jets[ijet].eta();
		Double_t dist2=getDistance(eta_p,eta_dete0,phi_p,phi_dete0);
		if(eta_p<etaMinCut || eta_p>etaMaxCut) continue; // fiducial acceptance 
		if(dist2<mindist2 && dist2<fR)
			{
				closest_par=ijet;
				mindist2=dist2;
			}
		}//2nd particle level jets loop
	
	if(closest_par==pjet)
		{
		matched2=true;
		
		}

	//-------------------------------
	//     FILLING RESPONSE MATRIX HISTOGRAMS - for matched jets
	//-------------------------------

	if(matched1 && matched2){
	//cout<<"Matched jets: "<<closest_par<<" "<<closest_dete<<endl;
	for(Int_t pTl=0; pTl<npTlead; pTl++)
		{
	   if(pTleading_dete0>pTl)
		{
			hJetEff2D[pTl]->Fill(pT_dete0, pT_par);//FILLING EFFICIENCY HISTOGRAM
		  	if(pTleading_par>pTl)
		  	{
			  double jetEratio = pT_dete0/pT_par;
			  if(jetEratio>0.15){
			  hresponse[pTl]->Fill(pT_dete0, pT_par);
			  hntrue_rm[pTl]->Fill(pT_par); //for normalization
			  hnpart_rm[pTl]->Fill(pT_par);
			  hndete_rm[pTl]->Fill(pT_dete0);
			  hndeterm_vs_eta[pTl]->Fill(eta_dete0,pT_dete0);
			  hntruerm_vs_eta[pTl]->Fill(eta_par,pT_par);
			  hJER[pTl]->Fill(pT_par, (pT_dete0-pT_par)/pT_par);
			  }
			}//part>pTleadmin
		}//dete>pTleadmin
  
		}//pTleading loop
	}//matched
	
	}//particle level jets loop
//cout<<"-----------"<<endl;

//}//end of embedding loop

}//end of event loop

	//calculate tracking efficiency
	for(int bn=1;bn<=hTrackEpsilon->GetNbinsX();bn++)
	{
		double track_pT = hTrackEpsilon->GetBinCenter(bn);
		int bin_true=hTrackEff->FindBin(-track_pT);
		int bin_dete=hTrackEff->FindBin(track_pT);
		double ntrue=hTrackEff->GetBinContent(bin_true);
		double ndete=hTrackEff->GetBinContent(bin_dete);
		if(ntrue>0)hTrackEpsilon->SetBinContent(bn,ndete/ntrue);
	}


//normalize Response matrix and efficiency histogram along x-axis (integral over measured pT for given generated pT is equal to 1 and total efficiency respectively)
float scaler=(float)1.0/fNjobs;
for(int pTl=0; pTl<npTlead; pTl++){

	//loop over particle level jet pT bins
	for(int j=1; j<nptbins+1; j++){
		double ntot_all=hntrue_all[pTl]->GetBinContent(j);
		double ntot_rm=hntrue_rm[pTl]->GetBinContent(j);
	//loop over detector level jet pT bins
	for(int i=1; i<nptbins+1; i++){
		double binc = hresponse[pTl]->GetBinContent(i,j);
		double rat;
		if(ntot_rm>0){
		  rat=binc/ntot_rm;
		  hresponse[pTl]->SetBinContent(i,j,rat);
		}

		binc = hJetEff2D[pTl]->GetBinContent(i,j);
		if(ntot_all>0) //we have a reasonable statistics
		{
			rat=binc/ntot_all;
		  	hJetEff2D[pTl]->SetBinContent(i,j,rat);
		}
		else 
		{
			//double trEff=hTrackEpsilon->GetBinContent(hTrackEpsilon->FindBin(pTl)); //tracking efficiency
			hJetEff2D[pTl]->SetBinContent(i,j,0); //SetBinContent(i,j,0) is not an optimal solution - it would decrease the efficiency when summing up several histograms. However when runing the simulation for large statistics  this condition should apply only on bins with pT<pTlead which we don't use anyway
		}
		
	}//i - detector level
	}//j - particle level
	hresponse[pTl]->Scale(scaler); //we will add more histograms together, so we have to divide the value by the number of histograms
	hJetEff2D[pTl]->Scale(scaler); //we will add more histograms together, so we have to divide the value by the number of histograms
	

	//make 1D projection of the efficiency histogram
	hJetEff1D[pTl]=(TH1D*) hJetEff2D[pTl]->ProjectionY(Form("hepsilon_pTl%i",pTl), 1, hJetEff2D[pTl]->GetXaxis()->GetNbins());
	hJetEff1D[pTl]->SetTitle("jet reconstruction efficiency");
	hRM_projectionY[pTl]=(TH1D*) hresponse[pTl]->ProjectionY(Form("hresponseY_pTl%i",pTl), 1, hresponse[pTl]->GetXaxis()->GetNbins());
	}//pTlead
	hTrackEpsilon->Scale(scaler);

  foutput->cd();
  if(foutput->IsOpen())cout<<"output file is OPEN"<<endl;
  else cout<<"WARNING - output file cannot be accessed"<<endl;
  foutput->Write();
  //ftreejets->Write();
  foutput->Close();

	delete detearr;
	delete pararr;
memstat.Stop();
return;
}
