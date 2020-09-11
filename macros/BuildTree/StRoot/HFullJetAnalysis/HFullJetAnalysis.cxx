
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


	outTree = new TTree("EventTree","EventTree");

	outTree->Branch("RunId",&RunId,"RunId/I");
	outTree->Branch("EventId",&EventId,"EventId/I");
  	outTree->Branch("refmult",&refmult,"refmult/S");
  	outTree->Branch("Vz",&Vz,"Vz/F");
	outTree->Branch("Vx",&Vx,"Vx/F");
	outTree->Branch("Vy",&Vy,"Vy/F");
	outTree->Branch("Vz_VPD",&Vz_VPD,"Vz_VPD/F");
	outTree->Branch("zdc",&zdc,"zdc/F");
	outTree->Branch("bbc",&bbc,"bbc/F");
	//outTree->Branch("NCharge",&NCharge,"NCharge/D");
	outTree->Branch("NBTOFMultfit",&NBTOFMultfit,"NBTOFMultfit/F");
	
	//TPC
	outTree->Branch("numTrk",&numTrk,"numTrk/I");

	
	outTree->Branch("Pt",&Pt,"Pt[numTrk]/F");
	//outTree->Branch("Px",&Px,"Px[numTrk]/D");
	//outTree->Branch("Py",&Py,"Py[numTrk]/D");
	outTree->Branch("Pz",&Pt,"Pz[numTrk]/F");
	outTree->Branch("Eta",&Eta,"Eta[numTrk]/F");
	outTree->Branch("Phi",&Phi,"Phi[numTrk]/F");
	//outTree->Branch("Charge",&Charge,"Charge[numTrk]/F");
	outTree->Branch("Dca",&Dca,"Dca[numTrk]/F");
	outTree->Branch("nHits",&nHits,"nHits[numTrk]/F");

	//EPD
	//outTree->Branch("EPDmaxHit",&EPDmaxHit, "EPDmaxHit/I");
    outTree->Branch("NumEPDHit",&NumEPDHit, "NumEPDHit/S");
    
    outTree->Branch("NumMip", NumMip, "NumMip[NumEPDHit]/F");
    //outTree->Branch("TileID", TileID, "TileID[NumEPDHit]/I");
    
    outTree->Branch("EPDX", EPDX, "EPDX[NumEPDHit]/F");
    outTree->Branch("EPDY", EPDY, "EPDY[NumEPDHit]/F");
    outTree->Branch("EPDZ", EPDZ, "EPDZ[NumEPDHit]/F");
	//outTree->Branch("EPD_eta", EPD_eta, "EPD_eta[NumEPDHit]/D");
	//outTree->Branch("EPD_phi", EPD_phi, "EPD_phi[NumEPDHit]/D");
	outTree->Branch("EPD_PP", EPD_PP, "EPD_PP[NumEPDHit]/F");
	outTree->Branch("EPD_TT", EPD_TT, "EPD_TT[NumEPDHit]/F");

	//BBC
	outTree->Branch("mBbcQ", mBbcQ, "mBbcQ[48]/I");
  	
	//Event plane
	//outTree->Branch("EPD_Psi",&EPD_Psi, "EPD_Psi/D");
	//outTree->Branch("TPC_Psi",&TPC_Psi, "TPC_Psi/D");


	return kStOK;
};

Int_t HFullJetAnalysis::Finish() {
    outTree->Print();
    fout->Write();
    fout->Close();
    return kStOK;
};


void HFullJetAnalysis::Clear(Option_t *opt) { };



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
	Int_t runnumber = picoEvent->runId();
	float refMultCorr =picoEvent->refMult();
//		vector<PseudoJet> particles;
//		particles.clear();


//	if(mRefMultCorr->isBadRun(picoEvent->runId())){cout<<"RefMultCorr says bad run"<<endl;return kStOk;}
/*	bool IsTrigger=false;
        IsTrigger=(picoEvent->isTrigger(450050)==1 ||
                   picoEvent->isTrigger(450060)==1 ||
		   picoEvent->isTrigger(450025)==1 ||
		   picoEvent->isTrigger(450015)==1 );


	if(!IsTrigger){return kStOk;}
*/


/*	
        bool CHECKRUN=1;
	for(size_t __i=0;__i<badrunlist.size(); __i++) if(picoEvent->runId()==badrunlist.at(__i) && CHECKRUN==1){CHECKRUN*=0; break;}
	if(!CHECKRUN){return kStOk;}
*/
	Int_t runidx;
        for(size_t __j=0;__j<runidlist.size();__j++) if(picoEvent->runId()==runidlist.at(__j)){runidx=__j;break;}
	//cout<<"runid"<<picoEvent->runId()<<endl;


	//event cut_____________________________________________
	//
	if(fabs(zTpc)>35.){return kStOk;}


    if(sqrt(pRcVx.x()*pRcVx.x()+pRcVx.y()*pRcVx.y()) >2.0){return kStOk;}

	

	if(fabs(zTpc-zVpd)>5){return kStOk;}

	//pile-up rejection_________________________
/*	double NBTOFMultfit =0;
	if( (picoEvent->nBTOFMatch())<14 ){
			NBTOFMultfit = 0.59207 + 2.1317*(picoEvent->nBTOFMatch());
		}else{
			NBTOFMultfit = 18.770 + 1.0699*(picoEvent->nBTOFMatch());
		}
	if( (picoEvent->refMult() ) >= NBTOFMultfit ) return kStOK;
*/
	//initialize_____________________________
	memset(Pt, 0, sizeof(Pt));
	memset(Px, 0, sizeof(Px));
	memset(Py, 0, sizeof(Py));
	memset(Pz, 0, sizeof(Pz));
	memset(Eta, 0, sizeof(Eta));
	memset(Phi, 0, sizeof(Phi));
	memset(Charge, 0, sizeof(Charge));
	memset(Dca, 0, sizeof(Dca));
	memset(nHits, 0, sizeof(nHits));



	double number_of_charge = 0;
	int NPt =0;


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
		
		//track cut____________________________________
		//
		if(pTrack->nHitsFit()<=15){ continue;}
		
        if(TMath::Abs(pTrack->gDCA(pRcVx).Mag())>=3) continue;
		if(TMath::Abs(eta)>1.5) continue;
	    if(pt<0.1) continue;
		if(pTrack->nHitsFit()*1./pTrack->nHitsMax()<0.52) continue;

        //_______________________________________

		if(charge!=0) number_of_charge++;
		

		Pt[NPt] = pt*charge;
		Px[NPt] = px;
		Py[NPt] = py;
		Pz[NPt] = pz;
		Eta[NPt] = eta;
		Phi[NPt] = phi;
		Charge[NPt] = charge;
		Dca[NPt] = dca;
		nHits[NPt] = pTrack->nHitsFit();

		NPt++;

	}
	//track loop end
	RunId 			= runnumber;
	EventId			= eventNumber;
	refmult         = refMultCorr;
	Vz              = pRcVx.z();
	Vx				= pRcVx.x();
	Vy				= pRcVx.y();
	zdc				= ZDC;
	bbc				= BBC;
	NCharge			= number_of_charge;
	numTrk			= NPt;
	Vz_VPD 			= zVpd;
	NBTOFMultfit	= picoEvent->nBTOFMatch();


	//EPD part___________________________________
	//_________________________________________________________________
	memset(NumMip, 0, sizeof(NumMip));
    memset(TileID, 0, sizeof(TileID));
    
    memset(EPDX, 0, sizeof(EPDX));
    memset(EPDY, 0, sizeof(EPDY));
    memset(EPDZ, 0, sizeof(EPDZ));
	memset(EPD_PP, 0, sizeof(EPD_PP));
	memset(EPD_TT, 0, sizeof(EPD_TT));
	memset(EPD_eta, 0, sizeof(EPD_eta));
	memset(EPD_phi, 0, sizeof(EPD_phi));

	for(int i=0;i<2;i++){
		for(int j=0;j<3;j++){
			QxPos[i][j]=0;
			QyPos[i][j]=0;
			QwPos[i][j]=0;
			QxNeg[i][j]=0;
			QyNeg[i][j]=0;
			QwNeg[i][j]=0;
		}
	}

	int NEPD = pico->numberOfEpdHits();
    cout<<"EPD hits "<<NEPD<<endl;

    EPDmaxHit 		= NEPD;

    NumEPDHit = 0;
    StEpdGeom* mEpdGeom = new StEpdGeom();

	for (int ihit=0; ihit<NEPD; ihit++){
		
		int tileId,ring,TT,PP,ADC;
        float nMip;
        
        StPicoEpdHit *epdHit = pico->epdHit(ihit);
		
		tileId = epdHit->id();
        //EW = (tileId<0)?0:1;
        ring = epdHit->row();
        TT = epdHit->tile();
        PP = epdHit->position();
        ADC = epdHit->adc();
        nMip = (TT<10)?(double)ADC/160.0:(double)ADC/115.0;

		double TileWeight = (nMip<2.0)?nMip:2.0;

		TVector3 StraightLine = mEpdGeom->TileCenter(tileId) - pRcVx;
        TVector3 StraightLine1 = mEpdGeom->TileCenter(tileId);
        
        double phi = StraightLine.Phi();
        double eta = StraightLine.Eta();
        
        float _x = StraightLine1.X();
        float _y = StraightLine1.Y();
        float _z = StraightLine1.Z();


		if(fabs(eta)<2.1 || fabs(eta)>5.1) continue;

		NumMip[NumEPDHit] = nMip;
        TileID[NumEPDHit] = tileId;
        
        EPDX[NumEPDHit] = _x;
        EPDY[NumEPDHit] = _y;
        EPDZ[NumEPDHit] = _z;
		EPD_PP[NumEPDHit] = PP;
		EPD_TT[NumEPDHit] = TT;
		EPD_eta[NumEPDHit] = eta;
		EPD_phi[NumEPDHit] = phi;
        
        NumEPDHit++;
		
		TVector3 XYZ(_x, _y, _z);
        
        float q_phi = XYZ.Phi();
        float q_eta = XYZ.Eta();
		bool IsPos;
        
        if(q_eta > 2.1 && q_eta < 5.1)   IsPos = true;
        if(q_eta > -5.1 && q_eta < -2.1) IsPos = false;
		/*
		float EpdWeight = getEpdWeight(_zbin, _cent, eta, phi, IsPos);
        TileWeight = TileWeight*EpdWeight;
        
        cout<<"EpdWeight = "<<EpdWeight<<" TileWeight= "<<TileWeight<<endl;
		*/
		/*
		for(int ih = 0; ih != 2; ih++){
            
            if(q_eta > -5.1 && q_eta < -3.6)
            {
                
                QwNeg[ih][0]+=TileWeight;
                QxNeg[ih][0] += TileWeight*cos((ih+2)*q_phi);
                QyNeg[ih][0] += TileWeight*sin((ih+2)*q_phi);
                
            }
            
            if(q_eta > -3.6 && q_eta < -2.1)
            {
                
                QwNeg[ih][1]+=TileWeight;
                QxNeg[ih][1] += TileWeight*cos((ih+2)*q_phi);
                QyNeg[ih][1] += TileWeight*sin((ih+2)*q_phi);
                
            }
            
            QxNeg[ih][2] = QxNeg[ih][0] + QxNeg[ih][1];
            QyNeg[ih][2] = QyNeg[ih][0] + QyNeg[ih][1];
            QwNeg[ih][2] = QwNeg[ih][0] + QwNeg[ih][1];
            
            if(q_eta > 2.1 && q_eta < 3.6)
            {
                
                QwPos[ih][0]+=TileWeight;
                QxPos[ih][0] += TileWeight*cos((ih+2)*q_phi);
                QyPos[ih][0] += TileWeight*sin((ih+2)*q_phi);
                
            }
            
            if(q_eta > 3.6 && q_eta < 5.1)
            {
                
                QwPos[ih][1]+=TileWeight;
                QxPos[ih][1] += TileWeight*cos((ih+2)*q_phi);
                QyPos[ih][1] += TileWeight*sin((ih+2)*q_phi);
                
            }
            
            QxPos[ih][2] = QxPos[ih][0] + QxPos[ih][1];
            QyPos[ih][2] = QyPos[ih][0] + QyPos[ih][1];
            QwPos[ih][2] = QwPos[ih][0] + QwPos[ih][1];
            
        }//ih
	*/

	}	//EPD loop end
/*
	Double_t Psi2 = TMath::ATan2(QyPos[0][2]+QyNeg[0][2],QxPos[0][2]+QxNeg[0][2]);
    Psi2 /= 2.0;

	EPD_Psi = Psi2;

  */  

	//BBC
	memset(mBbcQ,0,sizeof(mBbcQ));
  	for(int ch=0;ch<24;ch++) {
    	mBbcQ[ch]    = picoEvent->bbcAdcEast(ch);
    	mBbcQ[ch+24] = picoEvent->bbcAdcWest(ch);

  	}


	outTree->Fill();
	return kStOk;



};	



































