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

#include "pythia_fit_pars.h"
#include "dNdpT.h"
#include "ThrmSimPythiaJet.h"
#include "ThrmFourVector.h"
#include "ThrmRecPythia.h"
#include "ThrmMemStat.h"

//==========================================================================
ThrmSimPythiaJet::ThrmSimPythiaJet(TString wrkdir, bool jetonly, bool charge, bool efficorr, bool pTsmear) :
  fpartarr("ThrmFourVector", 1000)
{

  frag = gSystem->Getenv("JETFRAG"); //jet fragmentation: u | g | sp
  //eff_type = gSystem->Getenv("EFFTYPE"); //single hadron efficiency: AuAu | pp - like scenario
  //eff_path = gSystem->Getenv("EFFPATH"); //path to the single hadron efficiency functions
  hjettype=atoi(gSystem->Getenv("HARDJET_TYPE"));
  increment=atof(gSystem->Getenv("EFF_INCREMENT")); //vary efficiency for systematic studies
  kCentral=atoi(gSystem->Getenv("CENTRAL")); //central/peripheral collisions
  fR=atof(gSystem->Getenv("RADIUS")); //jet size R parameter
  RAA=atof(gSystem->Getenv("RAA")); //HARD JET SPECTRUM RAA
  RAApTdep=atoi(gSystem->Getenv("RAAPTDEP"));//make pT dependent RAA 0: no | 1: yes
  pTMAX=atof(gSystem->Getenv("PTMAX")); //track max pT
  pTMIN=atof(gSystem->Getenv("PTCUT")); //track min pT
  pTMINhard=atof(gSystem->Getenv("PTMINHARD")); //start of the hard distribution
  meanpT=atof(gSystem->Getenv("MEANPT")); //Background <pT>
	meanpT=meanpT/1000.0; //MeV->GeV
	TString outId=gSystem->Getenv("JOBID");
	kFixedSeed = atoi(gSystem->Getenv("FIXEDSEED")); //1: use fixed number (fSeed) as a seed for RN generator 0: use computer time as a seed
	kSeed = atoi(gSystem->Getenv("SEED")); //Seed for random number generatr

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


  // VALUES
  fNevents = 0;
  fNbinaryColl = 0;
  fMultiplicity = 0;
  fSigmaNbinaryColl = 0;
  fSigmaMultiplicity = 0;
  kJetOnly = jetonly;
  kThrmOnly = kFALSE;
  kCharged = charge;
  kEfficorr = efficorr;	
  kpTsmear = pTsmear;

  // FUNCTIONS
  Double_t pTmin = 0.2;
  Double_t pTmax = 200.;
//  Double_t pTjetMin = 1.;

  fbkgd = new TF1("fbkgd", "[0]*TMath::Power([1], 2)*x*TMath::Exp(-[1]*x)", pTmin, pTmax);
  fbkgd->SetParNames("Amplitude", "b (GeV/c)^{-1}");
/*
if(hjettype==0){
  fjet = new TF1("fjet", "[0]*TMath::Power(x, -[1])", pTjetMin, pTmax);
  fjet->SetParNames("Amplitude", "power");
}
else{
  fjet = new TF1("fjet", "[0]* TMath::Exp(-[1]*x)*TMath::Power(x,-[2])",pTjetMin,pTmax);
  fjet->SetParNames("Amplitude", "exponent", "power");
}*/
  Int_t npar = 2; //power law distribution
  if(hjettype==1) npar = 12; //double Levy fit to PYTHIA full jets * RAA
  else if(hjettype==2) npar = 4; //Tsalis fit to PYTHIA full jets * RAA
  dNdpT *fdNdpTjet = new dNdpT(npar,hjettype,RAApTdep); //hard jet distribution
  fjet = new TF1("fjet", fdNdpTjet, pTMINhard, pTmax, npar, "dNdpT");

  // PYTHIA
  fpythia = new TPythia6();
  fpythia->SetMRPY(1, seed);

  // OUTPUT
  foutput_path =  Form("%s/toymodel_events_%s.root", wrkdir.Data(),outId.Data());
  foutput_histo_path = Form("%s/control_histos_%s.root", wrkdir.Data(),outId.Data());
  ftree = new TTree("ToyModel","Toy Model Particle Production");
  TBranch *br_event = ftree->Branch("particles", &fpartarr);

  //tracking efficiency parameters
/*
  if(kEfficorr){
  efffile = new TFile(Form("%s/eff_%s.root",eff_path.Data(),eff_type.Data()),"OPEN");
  effL = (TF1*)efffile->Get("effhL");
  effH = (TF1*)efffile->Get("effhH");
  efffile->Close();
  }
*/

  // dNdpT histograms
  Int_t nbins = 400;
  Double_t xmin = -200;
  Double_t xmax = +200;
  TString name;
  TString title;

  name = Form("hfulljetpT");
  title = Form("soft+hard jet pT;p_{T} (GeV/c);Entries");
  fhalldNdpT = new TH1D(name.Data(), title.Data(), nbins, xmin, xmax);
  fhalldNdpT->Sumw2();

  name = Form("hpartpT");
  title = Form("jet constituent pT;p_{T} (GeV/c);Entries");
  fhjetpartdNdpT = new TH1D(name.Data(), title.Data(), nbins, xmin, xmax);
  fhjetpartdNdpT->Sumw2();

  name = Form("hjetgenpT");
  title = Form("generated hard jet pT;p_{T} (GeV/c);Entries");
  fhjetdNdpT = new TH1D(name.Data(), title.Data(), nbins, xmin, xmax);
  fhjetdNdpT->Sumw2();

  name = Form("hjetreqpT");
  title = Form("required hard jet pT;p_{T} (GeV/c);Entries");
  fhjetreqpT = new TH1D(name.Data(), title.Data(), nbins, xmin, xmax);
  fhjetreqpT->Sumw2();

  name = Form("hboltzdNdpT");
  title = Form("soft (BG) jet pT;p_{T} (GeV/c);Entries");
  fhboltzdNdpT = new TH1D(name.Data(), title.Data(), nbins, xmin, xmax);
  fhboltzdNdpT->Sumw2();

  name = Form("hpTpTleadReq");
  title = Form("requred pTleading vs jet pT;p_{T}^{leading} (GeV/c);p_{T}^{jet} (GeV/c);Entries");
  fhpTpTleading = new TH2D(name.Data(), title.Data(), nbins, 0, xmax, nbins, 0, xmax);
  fhpTpTleading->Sumw2();

  name = Form("hpTpTleadGen");
  title = Form("generated pTleading vs jet pT;p_{T}^{leading} (GeV/c);p_{T}^{jet} (GeV/c);Entries");
  fhpTpTleadingGen = new TH2D(name.Data(), title.Data(), nbins, 0, xmax, nbins, 0, xmax);
  fhpTpTleadingGen->Sumw2();
}

//==========================================================================
ThrmSimPythiaJet::~ThrmSimPythiaJet()
{
  //delete fdNdpTsp;
  //delete fdNdpTjet;

  delete fjet;
  delete fbkgd;


  delete fhalldNdpT;
  delete fhjetpartdNdpT;
  delete fhjetdNdpT;
  delete fhboltzdNdpT;
  delete fhjetreqpT;

  delete ftree;
  delete foutput;
  delete foutput_histo;

/*
 	if(kEfficorr){
	delete efffile;
	delete effL;
	delete effH;}
*/
}

//==========================================================================
void ThrmSimPythiaJet::Run()
{
ThrmMemStat memstat;
  memstat.Start();


  cout << "Running toy model with:" << endl;
  cout << Form("Nevts = %d", fNevents) << endl;
  
  cout << Form("Multiplicity = %d", fMultiplicity) << endl;
  cout << Form("Sigma Multiplicity = %d", fSigmaMultiplicity) << endl;

  cout << Form("Nbin = %d", fNbinaryColl) << endl;
  cout << Form("Sigma Nbin = %d", fSigmaNbinaryColl) << endl;  

	if(kThrmOnly)cout<<"Producing only thermal background."<<endl;
	if(kJetOnly)cout<<"Producing only hard jet distribution."<<endl;

  Int_t ievent = 0;
  while(ievent < fNevents)
    {
      if(ievent%100 == 0){
	cout << Form("Event #%d processed", ievent) << endl;
cout<<"Memory used 0: "<<memstat.Used()<<endl;}
      fpartarr.Delete();

      Int_t multiplicity;
      
		//cout<<"configuring event"<<endl;
      ConfigEvent(multiplicity);

		//cout<<"producing particles"<<endl;
      ProduceParticles(multiplicity);

		//cout<<"filling tree"<<endl;
      ftree->Fill();
      ievent++;
    }

   foutput_histo= new TFile(foutput_histo_path, "RECREATE");
   //foutput_histo->cd();
   fhalldNdpT->Write();
   fhjetpartdNdpT->Write();
   fhjetdNdpT->Write();
   fhboltzdNdpT->Write();
   fhjetreqpT->Write();
	fhpTpTleadingGen->Write();
	fhpTpTleading->Write();
   foutput_histo->Close();

   foutput= new TFile(foutput_path, "RECREATE");
  cout<<"creating output file: "<<foutput_path.Data()<<endl;
   //foutput->cd();
   ftree->Write();
   foutput->Close();

memstat.Stop();
}

//==========================================================================
void ThrmSimPythiaJet::ConfigEvent(Int_t &multiplicity)
{
  Double_t b = 0.0;
  Double_t power = 0.0;
  Double_t jetsignal = 0.0;
  Double_t TAA=0;

  Double_t pT0jet = 10.;
  Double_t pTjetmax = 200.;

  if(fkinematics == "RHIC" )
    {
      b = 2./meanpT; // <pT> = 500 MeV/c
      power = 6.0;
      jetsignal = 2.67e-06;
		TAA=22.5; //TAA for 0-10% most central collisions
		if(!kCentral)TAA=0.5; //60-80%
		//RAA=0.3;
    }
  
if(fkinematics == "RHIC_CHARGED")
   {
      b = 2./meanpT; // <pT> = 500 MeV/c
      power = 6.0;
      jetsignal = 2.67e-06;
		TAA=22.5;
		//RAA=0.3;
		if(!kCentral){
			TAA=0.5;
		}
   }
if(fkinematics == "LHC")
    {
      b = 2./meanpT; // <pT> = 700 MeV/c
      power = 5.0;
      jetsignal = 5.85e-04;
		TAA=20;
     }

  multiplicity = gRandom->Gaus(fMultiplicity, fSigmaMultiplicity);
  Int_t nbin = gRandom->Gaus(fNbinaryColl, fSigmaNbinaryColl);

	//Int_t ridx=fR*10-2; //0:R=0.2 1:R=0.3 2:R=0.4

  //Set parameters 
  Double_t signal = nbin*jetsignal;
  Double_t jet_norm = signal; 
  if(hjettype==0)//power law spectrum
		fjet->SetParameters(jet_norm,power);
  else if(hjettype==1) //double Levy fit to PYTHIA full jets
	{
  		 //TF1->SetParameters works only with max. 11 parameters... :-(
		fjet->SetParameter(0,TAA);
		fjet->SetParameter(1,B1);
		fjet->SetParameter(2,T1);
		fjet->SetParameter(3,n1);
		fjet->SetParameter(4,m1);
		fjet->SetParameter(5,mu1);
		fjet->SetParameter(6,B2);
		fjet->SetParameter(7,T2);
		fjet->SetParameter(8,n2);
		fjet->SetParameter(9,m2);
		fjet->SetParameter(10,mu2);
		fjet->SetParameter(11,RAA); 
 	}
  else if(hjettype==2) //Tsalis fit to PYTHIA full jets
	   fjet->SetParameters(TAA*A_tsal,n_tsal,T_tsal,RAA);

  //Double_t bkgd_norm = (multiplicity - fjet->Integral(pT0jet, pTjetmax)) / fbkgd->Integral(0.2, pTjetmax);
  Double_t bkgd_norm = multiplicity;
  fbkgd->SetParameters(bkgd_norm,b);

/*
TCanvas *cspectra = new TCanvas("cspectra","cspectra",10,10,800,600);
cspectra->SetLogy();
cspectra->SetGrid();
fjet->SetLineColor(kBlue);
fjet->SetNpx(1000);
fjet->Draw("same");*/
return;
}


//==========================================================================
void ThrmSimPythiaJet::ProduceParticles(Int_t multiplicity)
{
  Double_t pTthresh = pTMINhard;

  TClonesArray *simarr = new TClonesArray("ThrmFourVector", 1000);

  Int_t triggered = 0;
  Int_t goodparticle = 0;
  Int_t nhadrons = 0;
 
  // HARD JETS 
  
  if(!kThrmOnly){
  Double_t jetprob = fjet->Integral(pTthresh, fjet->GetXmax())*2*2*TMath::Pi();//jet probbability * eta range * phi range
	//cout<<"integral="<<fjet->Integral(pTthresh, fjet->GetXmax())<<endl;
	//cout<<"jetprob:"<<jetprob<<endl;
  Int_t jetmult = (Int_t) jetprob; 
//cout<<"Jet mult: "<<jetmult<<endl;
  //Int_t jetmult = 1;
  if(gRandom->Uniform() < jetprob - jetmult) jetmult++;
  
  while(triggered < jetmult)
    {
      simarr->Delete();      
      
      Double_t pT = fjet->GetRandom();
		//cout<<"jet pT: "<<pT<<endl;
      //Double_t pT = gRandom->Uniform(25,45);
		//if(pT < pTthresh) continue; //not necessary

      // RUN PYTHIA
      // MAKE PARTICLE CONE
      Double_t pTjetParton=pT;
		Double_t jetEta = gRandom->Uniform(-1., +1.);
		Double_t jetPhi = gRandom->Uniform(-1., +1.)*TMath::Pi();
	   if(frag!="sp"){	
			double pTjetPL = MakeCone2(fpythia, simarr, pT, jetEta, jetPhi, kCharged, kEfficorr,kpTsmear,frag);
		}

      
      // SAVE RECO JET
      Int_t Nparticles = 0;
      if(frag!="sp") Nparticles = simarr->GetEntries();
		else Nparticles=1;
		if(Nparticles==0)continue;
		triggered++;

      nhadrons += Nparticles;
      Float_t pTleading = 0;
      TLorentzVector total(0, 0, 0, 0);
   for(Int_t ipart = 0; ipart < Nparticles; ipart++)
	{
		if(frag!="sp")
		{
	  		ThrmFourVector *fv = (ThrmFourVector*)simarr->At(ipart);
			TLorentzVector lv = fv->GetTLorentzVector();
	  
		  	Double_t pTpart = lv.Pt();
			if(pTpart>pTMAX || pTpart<pTMIN)
			{
					  nhadrons=nhadrons-1;
					  continue;
			}
		  	Double_t eta = lv.Eta();
	 	 	Double_t phi = lv.Phi();
	  		Double_t M = lv.M();

	  		total += lv;

	  		if(pTpart > pTleading) pTleading = pTpart;
	  		fhjetpartdNdpT->Fill(pTpart);
	
	  		new (fpartarr[goodparticle]) ThrmFourVector(pTpart, eta, phi, M);
	  		goodparticle++;
		}
		else
		{
			pTleading = pT;
			fhjetpartdNdpT->Fill(pT);
			new (fpartarr[goodparticle]) ThrmFourVector(pT, jetEta, jetPhi, 0);
	  		goodparticle++;
		}
	}

    fhpTpTleading->Fill(pTleading, pTjetParton);
	 fhjetreqpT->Fill(pTjetParton);
    if(frag!="sp"){
		fhpTpTleadingGen->Fill(pTleading, total.Pt());
    	fhjetdNdpT->Fill(total.Pt());
	   fhalldNdpT->Fill(total.Pt());
	}
  }//triggered
  delete simarr;
 } //HARD JETS

  //BOLTZMAN BKG
 if (!kJetOnly)
 {
//Double_t npart=0;
   triggered = 0;
	 while(triggered < (multiplicity - nhadrons))
	{
	  Double_t pT = fbkgd->GetRandom();
		if(pT!=pT)continue;
	  if(pT < pTMIN) continue;
//cout<<triggered<<" pT: "<<pT<<endl;

//smearing and tracking eddiciency doesn't need to be applied to the background
/*	  
		if(kpTsmear)
      {
            Double_t sigma = 0.01*pT*pT;
            Double_t pTn = gRandom->Gaus(pT,sigma);
            //cout<<"old "<<ptp<<" new "<<ptn<<endl;
				if(pTn<0.2) continue;
            pT=pTn;
      }
     if(kEfficorr){
        Double_t epsilon=efficiency11(pT, effL, effH, increment, kCentral); //efficiency function from ThrmRecPythia.cxx
        Double_t rand = gRandom->Uniform(0,1);
//cout<<"epsilon "<<epsilon<<endl;
        if(rand>epsilon)continue;
     }
*/

	  Double_t eta = gRandom->Uniform(-1., +1.);
	  Double_t phi = gRandom->Uniform(-1., +1.)*TMath::Pi();
	  Double_t M = 0.0;
	  
	  new (fpartarr[triggered+nhadrons]) ThrmFourVector(pT, eta, phi, M);
	  
	  triggered++;
    fhboltzdNdpT->Fill(pT);
    fhalldNdpT->Fill(pT);
    //npart=triggered;
	}
     // cout<<"particles:"<<npart<<endl;
 }//BOLTZMAN BKG
return;
}


