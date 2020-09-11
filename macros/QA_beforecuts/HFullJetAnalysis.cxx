
//_____________________________________________
// full jet reconstruction reading TPC and BEMC
//
// Yang He  Jun11 2019



//___________________________________________________
//
//
#include "TNtuple.h"
#include "TF1.h"
#include "HFullJetAnalysis.h"


static const double Mass_pi = 0.140;
Double_t Jet_R=0.2;


ClassImp(HFullJetAnalysis)
HFullJetAnalysis::HFullJetAnalysis(
	const char*     Out,
	StPicoDstMaker* picoDstMaker,
	const char*     _bad_run_list,
	const char*     _run_number_list
	)
:StMaker ()
{

	picoReader = picoDstMaker;


	std::ifstream nfile ;
	nfile.open(_bad_run_list);
//	std::vector<int> badrunlist;
        std::string nline;
        while (std::getline(nfile, nline))
        {
                std::istringstream iss(nline);
                int nf;
                vector<int> p__;
                while (iss >> nf)
                {
                        if(nf>0)        p__.push_back(nf);
                }


                if(p__.size()<=0) continue;
                std::cout <<p__.at(0)<<"\n";
                badrunlist.push_back(p__.at(0));
                vector<int>().swap(p__);
        }
        nfile.close();
        for(size_t l_=0; l_<badrunlist.size(); l_++) std::cout <<badrunlist.at(l_)<<endl;


        std::ifstream infile ;
	infile.open(_run_number_list);
//	std::vector<int> runidlist;
        std::string nid;
        while (std::getline(infile, nid))
        {
                std::istringstream is(nid);
                int n1;
                vector<int>id_;
                while(is>>n1)
                {
                        if(n1>0)    id_.push_back(n1);
                }
                if(id_.size()<=0) continue;
                runidlist.push_back(id_.at(0));
                vector<int>().swap(id_);
        }
        infile.close();
	for(size_t l_=0; l_<runidlist.size(); l_++) std::cout <<runidlist.at(l_)<<endl;       
        
	//const Char_t *Out = "Picotest";
        infile2_.append(Out);
        infile2_.append(".root");
        cout<<"Will write to the file: "<<infile2_.c_str()<<endl;

};


HFullJetAnalysis::~HFullJetAnalysis()
{};


Int_t HFullJetAnalysis::Init() {
    fout = new TFile(infile2_.c_str(),"RECREATE");

	hVtxXvsY[0] = new TH2D("hVtxXvsY_0",
                        "hVtxXvsY before cuts",
                        400,-20.,20.,400,-20.,20.);
	hVtxXvsY[0]->GetXaxis()->SetTitle("Vx");
    hVtxXvsY[0]->GetYaxis()->SetTitle("Vy");

    hVtxXvsY[1] = new TH2D("hVtxXvsY_1",
                        "hVtxXvsY after cuts",
                        400,-20.,20.,400,-20.,20.);
	hVtxXvsY[1]->GetXaxis()->SetTitle("Vx");
    hVtxXvsY[1]->GetYaxis()->SetTitle("Vy");


	HVz[0] = new TH1D("Vz_0","Vz before cuts",1200,-300,300);
    HVz[0]->GetXaxis()->SetTitle("Vz");
    HVz[0]->GetYaxis()->SetTitle("Counts");

	HVz[1] = new TH1D("Vz_1","Vz after cuts",800,-200,200);
    HVz[1]->GetXaxis()->SetTitle("Vz");
    HVz[1]->GetYaxis()->SetTitle("Counts");



    Hdelta_Vz[0] = new TH1D("Hdelta_Vz_0","Hdelta_Vz before cuts",1000,0,100);
    Hdelta_Vz[0]->GetXaxis()->SetTitle("|V_{z,TPC}-V_{z,VPD}|");
    Hdelta_Vz[0] ->GetYaxis()->SetTitle("events");

    Hdelta_Vz[1] = new TH1D("Hdelta_Vz_1","Hdelta_Vz after cuts",100,0,10);
    Hdelta_Vz[1]->GetXaxis()->SetTitle("|V_{z,TPC}-V_{z,VPD}|");
    Hdelta_Vz[1] ->GetYaxis()->SetTitle("events");



    H2d_Vz[0] = new TH2D("H2d_Vz_0","H2d_Vz before cuts",800,-220,220,800,-220,220);
    H2d_Vz[0]->GetXaxis()->SetTitle("V_{z,TPC}");
    H2d_Vz[0]->GetYaxis()->SetTitle("V_{z,VPD}");

    H2d_Vz[1] = new TH2D("H2d_Vz_1","H2d_Vz after cuts",400,-120,120,400,-120,120);
    H2d_Vz[1]->GetXaxis()->SetTitle("V_{z,TPC}");
    H2d_Vz[1]->GetYaxis()->SetTitle("V_{z,VPD}");


    Htof_ref[0] = new TH2D("Htof_ref_0","Htof_ref before cuts",1000,0,1000,1000,0,1000);
    Htof_ref[0]->GetXaxis()->SetTitle("TOFMatch");
    Htof_ref[0]->GetYaxis()->SetTitle("refMult");

    Htof_ref[1] = new TH2D("Htof_ref_1","Htof_ref after cuts",1000,0,1000,1000,0,1000);
    Htof_ref[1]->GetXaxis()->SetTitle("TOFMatch");
    Htof_ref[1]->GetYaxis()->SetTitle("refMult");

    Htof_ntrk[0] = new TH2D("Htof_ntrk_0","Htof_ntrk before cuts",1000,0,1000,1000,0,1000);
    Htof_ntrk[0]->GetXaxis()->SetTitle("TOFMatch");
    Htof_ntrk[0]->GetYaxis()->SetTitle("charged track");

    Htof_ntrk[1] = new TH2D("Htof_ntrk_1","Htof_ntrk after cuts",1000,0,1000,1000,0,1000);
    Htof_ntrk[1]->GetXaxis()->SetTitle("TOFMatch");
    Htof_ntrk[1]->GetYaxis()->SetTitle("charged track");
    
	
	refMultHist[0] = new TH1D("refMult_0","refMult before cuts",800,-0.5,799.5);
    refMultHist[0] ->GetXaxis()->SetTitle("refMult");
    refMultHist[0] ->GetYaxis()->SetTitle("Counts");

	refMultHist[1] = new TH1D("refMult_1","refMult after cuts",800,-0.5,799.5);
    refMultHist[1] ->GetXaxis()->SetTitle("refMult");
    refMultHist[1] ->GetYaxis()->SetTitle("Counts");

        refMultHist[2] = new TH1D("refMult_2","refMult after all cuts",800,-0.5,799.5);
    refMultHist[2] ->GetXaxis()->SetTitle("refMult");
    refMultHist[2] ->GetYaxis()->SetTitle("Counts");

	vz_refmult[0]=new TH2D("vz_refmult_0","vz_refmult before cuts",1200,-300,300,800,0,800);
    vz_refmult[0]->GetXaxis()->SetTitle("Vz");
    vz_refmult[0]->GetYaxis()->SetTitle("refmult");

	vz_refmult[1]=new TH2D("vz_refmult_1","vz_refmult after cuts",800,-200,200,800,0,800);
    vz_refmult[1]->GetXaxis()->SetTitle("Vz");
    vz_refmult[1]->GetYaxis()->SetTitle("refmult");


	etahist[0] = new TH1D("etahist_0","etahist before cuts",200,-6,6);
    etahist[0]->GetXaxis()->SetTitle("#eta");
    etahist[0]->GetYaxis()->SetTitle("Counts");
    etahist[0]->Sumw2();

	etahist[1] = new TH1D("etahist_1","etahist after cuts",200,-6,6);
    etahist[1]->GetXaxis()->SetTitle("#eta");
    etahist[1]->GetYaxis()->SetTitle("Counts");
    etahist[1]->Sumw2();



    phihist[0] = new TH1D("phihist_0","phihist before cuts",300, -4, 4);
    phihist[0]->GetXaxis()->SetTitle("#phi");
    phihist[0]->GetYaxis()->SetTitle("Counts");
    phihist[0]->Sumw2();

	phihist[1] = new TH1D("phihist_1","phihist after cuts",300, -4, 4);
    phihist[1]->GetXaxis()->SetTitle("#phi");
    phihist[1]->GetYaxis()->SetTitle("Counts");
    phihist[1]->Sumw2();



	Pt_all_hist[0]= new TH1D("Pt_all_hist_0","pt distribution before cuts",200,0,100);
    Pt_all_hist[0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Pt_all_hist[0]->GetYaxis()->SetTitle("counts");
    Pt_all_hist[0]->Sumw2();

	Pt_all_hist[1]= new TH1D("Pt_all_hist_1","pt distribution after cuts",200,0,100);
    Pt_all_hist[1]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Pt_all_hist[1]->GetYaxis()->SetTitle("counts");
    Pt_all_hist[1]->Sumw2();


    Pthist_C= new TH1D("Pthist_C","after cut pt distribution 0-10%",200,0,100);
    Pthist_C->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Pthist_C->GetYaxis()->SetTitle("counts");
    Pthist_C->Sumw2();

	Pthist_P= new TH1D("Pthist_P","after cut pt distribution 60-80%",200,0,100);
    Pthist_P->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Pthist_P->GetYaxis()->SetTitle("counts");
    Pthist_P->Sumw2();

    Pthist_Cut0= new TH1D("Pthist_Cut0","before cut pt distribution",200,0,100);
    Pthist_Cut0->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Pthist_Cut0->GetYaxis()->SetTitle("counts");
    Pthist_Cut0->Sumw2();

    Pthist_Cut1= new TH1D("Pthist_Cut1","dca<3 pt distribution",200,0,100);
    Pthist_Cut1->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Pthist_Cut1->GetYaxis()->SetTitle("counts");
    Pthist_Cut1->Sumw2();

    Pthist_Cut2= new TH1D("Pthist_Cut2","dca<1 pt distribution",200,0,100);
    Pthist_Cut2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Pthist_Cut2->GetYaxis()->SetTitle("counts");
    Pthist_Cut2->Sumw2();

    Pthist_Cut3= new TH1D("Pthist_Cut3","dca<1 nhitsratio>0.52 pt distribution",200,0,100);
    Pthist_Cut3->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Pthist_Cut3->GetYaxis()->SetTitle("counts");
    Pthist_Cut3->Sumw2();
	
	nutrigger =new TH1D("nutrigger","trigger number distribution",8,0,8);
    nutrigger->GetXaxis()->SetTitle("trigger number");
    nutrigger->GetYaxis()->SetTitle("events");

	triggerdist_C[0] = new TH1D("triggerdist_C_0","0-10% 7-30 trigger pt distribution",80,0,40);
    triggerdist_C[0]->GetXaxis()->SetTitle("p_{T}^{trig} (GeV/c)");
    triggerdist_C[0]->GetYaxis()->SetTitle("counts");
	triggerdist_C[0]->Sumw2();

	triggerdist_P[0] = new TH1D("triggerdist_P_0","60-80% 7-30 trigger pt distribution",80,0,40);
    triggerdist_P[0]->GetXaxis()->SetTitle("p_{T}^{trig} (GeV/c)");
    triggerdist_P[0]->GetYaxis()->SetTitle("counts");
	triggerdist_P[0]->Sumw2();

    triggerdist_C[1] = new TH1D("triggerdist_C_1","0-10% 9-30 trigger pt distribution",80,0,40);
    triggerdist_C[1]->GetXaxis()->SetTitle("p_{T}^{trig} (GeV/c)");
    triggerdist_C[1]->GetYaxis()->SetTitle("counts");
	triggerdist_C[1]->Sumw2();

	triggerdist_P[1] = new TH1D("triggerdist_P_1","60-80% 9-30 trigger pt distribution",80,0,40);
    triggerdist_P[1]->GetXaxis()->SetTitle("p_{T}^{trig} (GeV/c)");
    triggerdist_P[1]->GetYaxis()->SetTitle("counts");
	triggerdist_P[1]->Sumw2();

	eventC =new TH1D("eventC","event number after cuts 0-10%",2,0,2);
    eventP =new TH1D("eventP","event number after cuts60-80%",2,0,2);

	HNCharge_C[0] = new TH1D("HNCharge_C_0","number of charged partiles before pileup in 0-10%",800,0,800);
	HNCharge_C[0]->GetXaxis()->SetTitle("N_{ch}");
	HNCharge_C[0]->GetYaxis()->SetTitle("events");

	HNCharge_P[0] = new TH1D("HNCharge_P_0","number of charged partiles before pileup in 60-80%",800,0,800);
	HNCharge_P[0]->GetXaxis()->SetTitle("N_{ch}");
	HNCharge_P[0]->GetYaxis()->SetTitle("events");

    HNCharge_C[1] = new TH1D("HNCharge_C_1","number of charged partiles after pileup in 0-10%",800,0,800);
	HNCharge_C[1]->GetXaxis()->SetTitle("N_{ch}");
	HNCharge_C[1]->GetYaxis()->SetTitle("events");

	HNCharge_P[1] = new TH1D("HNCharge_P_1","number of charged partiles after pileup in 60-80%",800,0,800);
	HNCharge_P[1]->GetXaxis()->SetTitle("N_{ch}");
	HNCharge_P[1]->GetYaxis()->SetTitle("events");


	return kStOK;
};

Int_t HFullJetAnalysis::Finish() {
    
    fout->Write();
    fout->Close();
    return kStOK;
};


void HFullJetAnalysis::Clear(Option_t *opt) { };

Int_t HFullJetAnalysis::centrality( float _refMult ){
    
    float   CentralityBins  [NCENT] = {213,149,102,68,44,27,15,8,4} ;  //
    Int_t   MiddleBinID     [NCENT] = { 0,1,2,3,4,5,6,7,8};  // ID Number
    
    Int_t   myCentrality=-1;
    
    for(int i=0;i!=NCENT;i++){
        if( _refMult > CentralityBins[i] ){
            myCentrality = MiddleBinID[i];
            break;
        }
        else{
            myCentrality = -1;
        }
    }
    
    return myCentrality ;
    
}

bool HFullJetAnalysis::refmult_check (int __nBTOFMatch, int __refMult){

//We do the Y-projection in a nBTOFMatch window, fit by double negative binomial distribution then get these parameters.
///Recommand use max as 3, set min as 4.

double a0=6.44225, a1=1.37398, a2=-0.00446098, a3=2.60578e-05, a4= -6.84022e-08, a5=6.18087e-11;
//double b0=2.52126730672253, b1=0.128066911940844, b2=-0.000538959206681944, b3=1.21531743671716e-06, b4=-1.01886685404478e-09;
//double c0=4.79427731664144, c1=0.187601372159186, c2=-0.000849856673886957, c3=1.9359155975421e-06, c4=-1.61214724626684e-09;



double refmultcutmax=a0+a1*(__nBTOFMatch)+a2*pow(__nBTOFMatch,2)+a3*pow(__nBTOFMatch,3)+a4*pow(__nBTOFMatch,4)+a5*pow(__nBTOFMatch,5);


        //cout<<"refmult,maxcut,mincut= "<<__refMult<<" "<<refmultcutmax<<" "<<refmultcutmin<<endl;

        if( __refMult<refmultcutmax){
                return true;
        }else{
                return false;
        }

}

bool HFullJetAnalysis::ntrk_check(int __nBTOFMatch, int __ntrk){

    double b0=-13.8407, b1=1.00334, b2=0.000421093, b3=-1.88309e-06, b4=2.27559e-09;
    double c0=10.4218, c1=1.84005, c2=-0.00289939, c3=1.01996e-05, c4=-1.472e-08;

    double ntrkcutmin=b0+b1*(__nBTOFMatch)+b2*pow(__nBTOFMatch,2)+b3*pow(__nBTOFMatch,3)+b4*pow(__nBTOFMatch,4);
    double ntrkcutmax=c0+c1*(__nBTOFMatch)+c2*pow(__nBTOFMatch,2)+c3*pow(__nBTOFMatch,3)+c4*pow(__nBTOFMatch,4);

    if( __nBTOFMatch<450){
    if( __ntrk<ntrkcutmax && __ntrk>ntrkcutmin){
                return true;
        }else{
                return false;
        }
    }

    if( __nBTOFMatch>450){
    if( __ntrk>ntrkcutmin){
                return true;
        }else{
                return false;
        }
    }



}


Int_t HFullJetAnalysis::Make() {

	if ( !picoReader) return kStWarn;
	StPicoDst *pico = picoReader->picoDst();
	StPicoEvent *picoEvent = pico->event();
    if( !picoEvent ) {return kStWarn;}


    TVector3 pRcVx = picoEvent->primaryVertex();


    Double_t zTpc = pRcVx.z();
	

    Double_t zVpd = picoEvent->vzVpd();

                
	Double_t ZDC=picoEvent->ZDCx();
	Double_t BBC=picoEvent->BBCx();
	

	Int_t NoGlobalTracks = pico->numberOfTracks();
		
	Int_t gRefMult=picoEvent->grefMult();
/*		StRefMultCorr* mRefMultCorr = CentralityMaker::instance()->getgRefMultCorr() ;
		mRefMultCorr->init((Int_t)picoEvent->runId());
                mRefMultCorr->initEvent(gRefMult,zTpc,ZDC);
                float refMultCorr = mRefMultCorr->getRefMultCorr() ;
*/
	Int_t eventNumber=picoEvent->eventId();
	float refMultCorr =picoEvent->refMult();
//		vector<PseudoJet> particles;
//		particles.clear();

    Int_t CentId = centrality(refMultCorr);

//	if(mRefMultCorr->isBadRun(picoEvent->runId())){cout<<"RefMultCorr says bad run"<<endl;return kStOk;}
	bool IsTrigger=false;
        IsTrigger=(picoEvent->isTrigger(600001)==1 ||
                   picoEvent->isTrigger(600011)==1 ||
		   picoEvent->isTrigger(600021)==1 ||
		   picoEvent->isTrigger(600031)==1 );


	if(!IsTrigger){return kStOk;}



	
    bool CHECKRUN=1;
	for(size_t __i=0;__i<badrunlist.size(); __i++) if(picoEvent->runId()==badrunlist.at(__i) && CHECKRUN==1){CHECKRUN*=0; break;}
	if(!CHECKRUN){return kStOk;}
/*
	Int_t runidx;
        for(size_t __j=0;__j<runidlist.size();__j++) if(picoEvent->runId()==runidlist.at(__j)){runidx=__j;break;}
	//cout<<"runid"<<picoEvent->runId()<<endl;
*/

    //event plots before cuts
    hVtxXvsY[0]->Fill(pRcVx.x(),pRcVx.y());
    HVz[0]->Fill(zTpc);
    refMultHist[0]->Fill(refMultCorr);
    vz_refmult[0]->Fill(zTpc,refMultCorr);
    Htof_ref[0]->Fill(picoEvent->nBTOFMatch(),refMultCorr);
    H2d_Vz[0]->Fill(zTpc,zVpd);
    Hdelta_Vz[0]->Fill(fabs(zTpc-zVpd));
	//event cut_____________________________________________
	//
	if(zTpc>25.||zTpc<-35){return kStOk;}


    if(sqrt(pRcVx.x()*pRcVx.x()+pRcVx.y()*pRcVx.y()) >2.0){return kStOk;}

	

	if(fabs(zTpc-zVpd)>5){return kStOk;}

    //Event QA plots
    hVtxXvsY[1]->Fill(pRcVx.x(),pRcVx.y());
    HVz[1]->Fill(zTpc);
    //refMultHist[1]->Fill(refMultCorr);
    vz_refmult[1]->Fill(zTpc,refMultCorr);
    H2d_Vz[1]->Fill(zTpc,zVpd);
    Hdelta_Vz[1]->Fill(fabs(zTpc-zVpd));
	//pile-up rejection_________________________
/*	double NBTOFMultfit =0;
	if( (picoEvent->nBTOFMatch())<14 ){
			NBTOFMultfit = 0.59207 + 2.1317*(picoEvent->nBTOFMatch());
		}else{
			NBTOFMultfit = 18.770 + 1.0699*(picoEvent->nBTOFMatch());
		}
	if( (picoEvent->refMult() ) >= NBTOFMultfit ) return kStOK;
*/
    if(!refmult_check(picoEvent->nBTOFMatch(),refMultCorr)) return kStOk;

    
    Htof_ref[1]->Fill(picoEvent->nBTOFMatch(),refMultCorr);
    //refMultHist[1]->Fill(refMultCorr);
    //tof vs ntrk pile up
    int Nch=0;
    for(Int_t iTrk=0; iTrk<NoGlobalTracks; iTrk++) {

		StPicoTrack *gTrack = pico->track(iTrk);

        if(!gTrack) continue;
	
		if(!gTrack->isPrimary()) continue;

        Double_t pt=gTrack->pMom().Perp();
        Double_t eta=gTrack->pMom().PseudoRapidity();
		Double_t dca=TMath::Abs(gTrack->gDCA(pRcVx).Mag());
		
		//track cut____________________________________
		//
		if(gTrack->nHitsFit()<=15){ continue;}
		
        if(dca>=1) continue;
		if(TMath::Abs(eta)>1.) continue;
	    if(pt<0.2||pt>30) continue;
		if(gTrack->nHitsFit()*1./gTrack->nHitsMax()<0.52) continue;

        Nch++;
    }

    Htof_ntrk[0]->Fill(picoEvent->nBTOFMatch(),Nch);
    if(CentId==0) HNCharge_C[0]->Fill(Nch);
    if(CentId==6||CentId==7) HNCharge_P[0]->Fill(Nch);

    if(ntrk_check(picoEvent->nBTOFMatch(),Nch)) refMultHist[2]->Fill(refMultCorr);
    if(!ntrk_check(picoEvent->nBTOFMatch(),Nch)) return kStOk;
  
    //QA plots
    refMultHist[1]->Fill(refMultCorr);
    Htof_ntrk[1]->Fill(picoEvent->nBTOFMatch(),Nch);
    if(CentId==0){
        HNCharge_C[1]->Fill(Nch);
        eventC->Fill(1); }
    if(CentId==6||CentId==7){
        HNCharge_P[1]->Fill(Nch);
        eventP->Fill(1); }




	double number_of_charge = 0;
	int NPt =0;
    TVector3 Trigger[8];
    TVector3 Trigger1[8];
	
    Int_t trigger=0;
	Int_t trigger1=0;


	//track loop
	for(Int_t iTrk=0; iTrk<NoGlobalTracks; iTrk++) {

		StPicoTrack *pTrack = pico->track(iTrk);

        if(!pTrack) continue;
	
		if(!pTrack->isPrimary()) continue;


		//tofmatch_______________________________________________________________
	/*	if( pTrack->isTofTrack() ) {
                StPicoBTofPidTraits *trait = pico->btofPidTraits( pTrack->bTofPidTraitsIndex() );
                if(TMath::Abs(pTrack->gDCA(pRcVx).Mag())<2 && fabs(pTrack->pMom().PseudoRapidity())<0.5 &&
                             pTrack->nHitsFit()>10 && trait->btofBeta()>0){TOFMatch ++;}
                }
		*/

		Double_t pt=pTrack->pMom().Perp();
        Double_t eta=pTrack->pMom().PseudoRapidity();
        Double_t phi=pTrack->pMom().Phi();
        Double_t charge=pTrack->charge();
        Double_t px=pTrack->pMom().x();
        Double_t py=pTrack->pMom().y();
        Double_t pz=pTrack->pMom().z();
		Double_t dca=TMath::Abs(pTrack->gDCA(pRcVx).Mag());
		
        //before cut plots
        Pt_all_hist[0]->Fill(pt);
        etahist[0]->Fill(eta);
        phihist[0]->Fill(phi);
        Pthist_Cut0->Fill(pt);

		//track cut____________________________________
		//
		if(pTrack->nHitsFit()<=15){ continue;}
		
		if(TMath::Abs(eta)>1.) continue;
	    if(pt<0.2) continue;

        if(dca<3) Pthist_Cut1->Fill(pt); //pt cut test

        if(dca>=1) continue;

        Pthist_Cut2->Fill(pt); //pt cut test


		if(pTrack->nHitsFit()*1./pTrack->nHitsMax()<0.52) continue;

        //_______________________________________
        Pt_all_hist[1]->Fill(pt);
        etahist[1]->Fill(eta);
        phihist[1]->Fill(phi);
		if(CentId==0) Pthist_C->Fill(pt);
        if(CentId==6||CentId==7) Pthist_P->Fill(pt);
        Pthist_Cut3->Fill(pt);


		NPt++;

        if(pt<30&&pt>=7){
		    trigger +=1;
		    Trigger[trigger-1].SetPtEtaPhi(pt,eta,phi);
		}

	    if(pt<30&&pt>=9){
            trigger1 +=1;
            Trigger1[trigger1-1].SetPtEtaPhi(pt,eta,phi);
        }


	}
	//track loop end
	
    nutrigger->Fill(trigger);

    double triggerpt=0.;
    if(trigger>0){
        if(trigger==1){
            triggerpt=Trigger[0].Perp();
			if(CentId==0) triggerdist_C[0]->Fill(triggerpt);
            if(CentId==6||CentId==7) triggerdist_P[0]->Fill(triggerpt);
	    }
	    else{
            Int_t z=gRandom->Uniform(0,trigger);
            triggerpt=Trigger[z].Perp();
			if(CentId==0) triggerdist_C[0]->Fill(triggerpt);
            if(CentId==6||CentId==7) triggerdist_P[0]->Fill(triggerpt);
		}
	}

    triggerpt=0.;
    if(trigger1>0){
        if(trigger1==1){
            triggerpt=Trigger1[0].Perp();
			if(CentId==0) triggerdist_C[1]->Fill(triggerpt);
            if(CentId==6||CentId==7) triggerdist_P[1]->Fill(triggerpt);
	    }
	    else{
            Int_t z=gRandom->Uniform(0,trigger1);
            triggerpt=Trigger1[z].Perp();
			if(CentId==0) triggerdist_C[1]->Fill(triggerpt);
            if(CentId==6||CentId==7) triggerdist_P[1]->Fill(triggerpt);
		}
	}




	return kStOk;



};	



































