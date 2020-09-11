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
//#include "fjwrapper.h"
//#include <fastjet/PseudoJet.hh>
#include "ThrmPythiaTest.h"

#include <fastjet/config.h>             // will allow a test for FJ3
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/Subtractor.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>

using namespace fastjet;
using namespace std;

//=================================================

ThrmPythiaTest::ThrmPythiaTest()
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
TString outDir = gSystem->Getenv("OUTPUTDIR");
//TString eff_path = gSystem->Getenv("EFFPATH"); //path to the single hadron efficiency functions
//TString eff_type = gSystem->Getenv("EFFTYPE"); //single hadron efficiency: AuAu | pp - like scenario
fweight = atoi(gSystem->Getenv("WEIGHT"));

pTmin=3.0;
pTmax=50;

TString chname="full";
if(kCharged) chname="charged";

//OUTPUT
TString str = Form("%s/pythia_%s_test_flat_%.1lf_%s.root", outDir.Data(), chname.Data(), fR, frag.Data());
if(fweight==1)str = Form("%s/pythia_%s_test_pT5_%.1lf_%s.root", outDir.Data(), chname.Data(),fR, frag.Data());
foutput = new TFile(str.Data(), "RECREATE");
//ftreejets = new TTree("PythiaEmbJets","Pythia Jets");
//TBranch *br_pythia = ftreejets->Branch("pythia", &fembeddingarr);

//set acceptance
etaMinCut=-(fmax_rap-fR);
etaMaxCut=(fmax_rap-fR);
//etaMinCut=-1.;
//etaMaxCut=+1.;


//histogram definitions
nptbins=250;
Float_t ptminbin=0;
Float_t ptmaxbin=pTmax;

for(Int_t pTl=0; pTl<npTlead; pTl++)
{
  TString name=Form("hpTptlJet_pTl%i",pTl);
  hpTptlJet[pTl]=new TH1D(name, "PLJ n generated;p_{T}^{true};entries", nptbins, ptminbin, ptmaxbin);
  name=Form("hpTptlJet0_pTl%i",pTl);
  hpTptlJet0[pTl]=new TH1D(name, "PLJ n generated;p_{T}^{true};entries", nptbins, ptminbin, ptmaxbin);
  name=Form("hpTptl4v_pTl%i",pTl);
  hpTptl4v[pTl]=new TH1D(name, "PLJ n generated;p_{T}^{true};entries", nptbins, ptminbin, ptmaxbin);
}
hpTjet=new TH1D("hpTjet", "pT of PLJ reconstructed;p_{T}^{true};entries", nptbins, ptminbin, ptmaxbin);
hpt4v1=new TH1D("hpt4v1", "pT of PLJ 4vector;p_{T}^{true};entries", nptbins, ptminbin, ptmaxbin);
hpt4v2=new TH1D("hpt4v2", "pT of PLJ 4vector 2;p_{T}^{true};entries", nptbins, ptminbin, ptmaxbin);

hdpt_nparts=new TH2D("hdpt_nparts","dpt vs n parts(in 4vector)",25,-20,5,10,0.5,10.5);
hdpt_pt=new TH2D("hdpt_pt","dpt vs pT(of 4vector)",25,-20,5,30, ptminbin, ptmaxbin);
hdpt_z=new TH2D("hdpt_z","dpt vs z of pTleading",25,-20,5,21, 0, 1.05);
hpT_pTdete_pTl5=new TH2D("hpT_pTdete_pTl5","pt4v vs pTjet with pTleading>5",nptbins, ptminbin, ptmaxbin,nptbins, ptminbin, ptmaxbin);

hz_pT=new TH2D("hz_pT","FF vs jet pT",51,0,1.02,nptbins, ptminbin, ptmaxbin);
hnJets=new TH1D("hnJets","# of jets with given pT",nptbins, ptminbin, ptmaxbin);

/*
if(kEfficorr)
{
//Tracking efficiency
TFile* efffile = new TFile(Form("%s/eff_%s.root",eff_path.Data(),eff_type.Data()),"OPEN");
effL = (TF1*)efffile->Get("effhL");
effH = (TF1*)efffile->Get("effhH");
efffile->Close();
}
*/


// SETTING RANDOM SEED
  TDatime dt;
  UInt_t curtime = dt.Get();
  UInt_t procid = gSystem->GetPid();
  UInt_t seed = curtime - procid;
  gRandom->SetSeed(seed);

// PYTHIA
fpythia = new TPythia6();
fpythia->SetMRPY(1, seed);

}

//=================================================
ThrmPythiaTest::~ThrmPythiaTest()
{
delete effL;
delete effH;
delete fpythia; 

foutput->Close();
delete ftreejets;
delete foutput;
}

//=================================================
double ThrmPythiaTest::getDistance(Double_t eta1, Double_t eta2,Double_t phi1, Double_t phi2)
{
		Double_t delta_eta=TMath::Abs(eta1-eta2);
		Double_t delta_phi=TMath::Abs(phi1-phi2);
		if(delta_phi>TMath::Pi()) delta_phi=2*TMath::Pi()-delta_phi;
		double dist=TMath::Sqrt(delta_eta*delta_eta+delta_phi*delta_phi);
		return dist;
}

//=================================================
void ThrmPythiaTest::Run(Int_t nevents)
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
	detearr->Delete();
	pararr->Delete();

	//Double_t jetEta = gRandom->Uniform(-fmax_rap, +fmax_rap);
	Double_t jetEta = 0;
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
	Double_t PtJet=MakeCone(fpythia, detearr, pararr, pT, jetEta, jetPhi, kCharged, kEfficorr,kpTsmear,usefrag); 

//cout<<"req: "<<pT<<" gen: "<<PtJet<<endl;
		float weight=1.0;
		if(fweight==1) weight=TMath::Power(pT,-5);

	float pTlead4v=0;
	//loop over particles in particle-jet
	Int_t Nparticles = pararr->GetEntries();
	//Float_t pTleading = 0;
	TLorentzVector total(0, 0, 0, 0);
	   vector<PseudoJet> input_vector;
		hnJets->Fill(PtJet);
	for(Int_t ipart = 0; ipart < Nparticles; ipart++)
  	{
     ThrmFourVector *fv = (ThrmFourVector*)pararr->At(ipart);
     TLorentzVector lv = fv->GetTLorentzVector();
     
		//if(lv.Pt() < fpTcut) continue;
      //if(TMath::Abs(lv.Eta()) > fmax_rap) continue;

     Double_t pTpart = lv.Pt();
     Double_t ppart = lv.P();
     total += lv;
     if(pTpart > pTlead4v) pTlead4v = pTpart;
		hz_pT->Fill(ppart/PtJet,PtJet,ppart/PtJet);
      PseudoJet inp_particle(lv.Px(),lv.Py(),lv.Pz(),lv.Energy()); 
		inp_particle.set_user_index(ipart);
		input_vector.push_back(inp_particle);
	
	}

//cout<<"===================================================="<<endl;
//cout<<"4vector pT:"<<PtJet<<endl;
/*
if(pTlead4v-total.Pt()>0.2 && pTlead4v>4)
{
	cout<<"WARNING !!!!!!!!! pTlead4v>pT4v."<<endl;
	cout<<"pT="<<total.Pt()<<" pTlead="<<pTlead4v<<" eta:"<<total.Eta()<<" phi:"<<total.Phi()<<endl;
	for(Int_t ipart = 0; ipart < Nparticles; ipart++)
   {
     ThrmFourVector *fv = (ThrmFourVector*)pararr->At(ipart);
     TLorentzVector lv = fv->GetTLorentzVector();

      if(lv.Pt() < fpTcut) continue;
      if(TMath::Abs(lv.Eta()) > fmax_rap) continue;

     Double_t pTpart = lv.Pt();
		Double_t eta = lv.Eta();
		Double_t phi=lv.Phi();
		
	cout<<"part "<<ipart<<" pT="<<lv.Pt()<<" eta:"<<lv.Eta()<<" phi:"<<lv.Phi()<<endl;

   }

}*/
float pT4v=total.Pt();
		hpt4v2->Fill(pT4v);
		for(Int_t pTl=0; pTl<npTlead; pTl++){
			if(pTlead4v>pTl ){
				hpTptl4v[pTl]->Fill(pT4v,weight);
			}
		}
//cout<<"4vector pT II:"<<pT4v<<" pTlead_4v:"<<pTlead4v<<endl;

//FASTJET

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

   for(Int_t pjet=0; pjet<jets.size(); pjet++){
   
      Double_t phi_par = jets[pjet].phi();
      Double_t eta_par = jets[pjet].eta();
      Double_t pT_par = jets[pjet].perp();
      Double_t M_par = jets[pjet].m();
      vector<PseudoJet> constituents = sorted_by_pt(clust_seq_hard.constituents(jets[pjet]));
		Double_t pTleading_par = constituents[0].perp();

      //cout<<pjet<<"reco:"<<pT_par<<" pTlead:"<<pTleading_par<<"  eta:"<<eta_par<<" phi:"<<phi_par<<endl;


//FILLING HISTOGRAMS
		hpTjet->Fill(pT);
		hpt4v1->Fill(PtJet);
		for(Int_t pTl=0; pTl<npTlead; pTl++){
			if(pTleading_par>pTl )
			{
				hpTptlJet[pTl]->Fill(pT_par,weight);
				if(pjet==0)hpTptlJet0[pTl]->Fill(pT_par,weight);
			}
		}//loop over pTlead cuts

float dpT=pT_par-pT4v;
float z=pTlead4v/pT4v;
hdpt_nparts->Fill(dpT,Nparticles,weight);
hdpt_pt->Fill(dpT,pT4v,weight);
hdpt_z->Fill(dpT,z,weight);
if(pTleading_par>5) hpT_pTdete_pTl5->Fill(pT4v,pT_par,weight);

}//jets



/*
// background estimation
                            //cout << "Define JetDefinition" << endl;
                            //JetDefinition jet_def_bkgd(kt_algorithm, 0.4);
                            JetDefinition jet_def_bkgd(kt_algorithm, jet_R_background); // <--
                            //JetDefinition jet_def_bkgd(antikt_algorithm, jet_R); // test
                            //cout << "Define AreaDefinition" << endl;
                            AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
                            //AreaDefinition area_def_bkgd(active_area,GhostedAreaSpec(ghost_maxrap,1,0.005));
                            //cout << "Define selector" << endl;
                            //Selector selector = SelectorAbsRapMax(1.0) * (!SelectorNHardest(Remove_N_hardest)); // 2
                            Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(Remove_N_hardest)); // <--
                            //Selector selector = SelectorAbsEtaMax(1.0 - jet_R); // test


                            //cout << "Define JetMedianBackgroundEstimator" << endl;
                            JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd); // <--
                            //JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def, area_def); // test
                            //cout << "Define Subtractor" << endl;
                            Subtractor subtractor(&bkgd_estimator);
                            //cout << "Define bkgd_estimator" << endl;
                            bkgd_estimator.set_particles(particles[i_orig_smear]);

                            //cout << "Calculate jet_rho and jet_sigma" << endl;
                            Double_t jet_rho   = bkgd_estimator.rho();
*/
/*
	FJWrapper kt_data;
   kt_data.r = fR;
   kt_data.maxrap = fmax_rap;
   kt_data.algor = fastjet::kt_algorithm;
   kt_data.input_particles = input_vector;
   kt_data.Run();

   std::vector<fastjet::PseudoJet> ktjets = kt_data.inclusive_jets;
	double  rho, sigma;
	kt_data.GetMedianAndSigma(rho, sigma,0);

	FJWrapper akt_data;
   akt_data.r = fR;
   akt_data.maxrap = fmax_rap;
   akt_data.algor = fastjet::antikt_algorithm;
   akt_data.input_particles = input_vector;
   akt_data.Run();

   //std::vector<fastjet::PseudoJet> jets = sorted_by_pt(akt_data.inclusive_jets);
   std::vector<fastjet::PseudoJet> jets = akt_data.inclusive_jets;

	Int_t Njets_par=jets.size();
	//if(Njets_par<1)continue;

	
	//cout<<"rnd:"<<rnd<<" pjet: "<<pjet<<" njets: "<<Njets_par<<endl;
	//cout<<" pjet: "<<pjet<<" njets: "<<Njets_par<<endl;

	for(Int_t pjet=0; pjet<jets.size(); pjet++){
	
      Double_t phi_par = jets[pjet].phi();
      Double_t eta_par = jets[pjet].eta();
      Double_t pT_par = jets[pjet].perp();
      Double_t M_par = jets[pjet].m();
      Double_t area_par = akt_data.clust_seq->area(jets[pjet]);
      std::vector<fastjet::PseudoJet> constituents = sorted_by_pt(akt_data.clust_seq->constituents(jets[pjet]));
      Double_t pTleading_par = constituents[0].perp();

		cout<<pjet<<"reco:"<<pT_par<<" pTlead:"<<pTleading_par<<"  eta:"<<eta_par<<" phi:"<<phi_par<<endl;

		if(eta_par<etaMinCut || eta_par>etaMaxCut) continue; // fiducial acceptance 
	   //FILLING EFFICIENCY HISTOGRAMS
		hpTjet->Fill(pT);
		hpt4v1->Fill(PtJet);
		for(Int_t pTl=0; pTl<npTlead; pTl++){
			if(pTleading_par>pTl )
			{
				//heffi_par[pTl]->Fill(pT_par);
				hpTptlJet[pTl]->Fill(pT_par,weight);

				//heffi_par[pTl]->Fill(pTleading_par);
			}
		}//loop over pTlead cuts
}//jet loop
*/

}//event loop
  foutput->cd();
  foutput->Write();
  //ftreejets->Write();
  foutput->Close();

	delete detearr;
	delete pararr;
memstat.Stop();
return;
}
