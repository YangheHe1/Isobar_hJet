#include "TROOT.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TMath.h"
#include "TFile.h"
#include "TH2.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TComplex.h"

#include <algorithm>
#include <TRandom.h>
#include "readtree.C"

using namespace std;


const Int_t maxTRACK=10000;
const Int_t NCENT = 9; //centrality

const int NHAR = 2; //Psin
const double PI=acos(-1.0);
const int NK = 12;  //flattening kmax

//for class bin
const Int_t refMult_bin=4;
const Int_t Psi_bin=2;
const Int_t Vz_bin = 2;

Float_t mult_start_stop_delta[3] = {181.0,616.0,1.0};
Float_t Psi_start_stop_delta[3] = {-TMath::Pi()/2.0,TMath::Pi()/2.0,1.0};
Float_t vertex_z_start_stop_delta[3] = {-40.0,40.0,1.0};


Double_t QxRaw[2];
Double_t QyRaw[2];
Double_t QwRaw[2];


//switch mode
char *FileInput=0;
char *FileInput1=0;
char *FileOutput=0;

Int_t Cen_mode;     //0 all; 1 central; 2 peripheral 
Int_t Psi_mode;     //0 no psi bin;1 epd psi; 2 tpc psi


ifstream inData;
ofstream outData;


//_______class file___________________________
TFile*          fout[refMult_bin][Vz_bin][Psi_bin];

TTree          *outEventTree[refMult_bin][Vz_bin][Psi_bin];

Double_t    M_refmult;
Double_t    M_Vx;
Double_t    M_Vy;
Double_t    M_Vz;
Double_t    M_Psi2;
Int_t       M_numTrk;
Double_t    M_TofMatch;
	
Double_t		M_Pt[maxTRACK];
Double_t		M_Eta[maxTRACK];
Double_t		M_Phi[maxTRACK];
Double_t  	    M_Charge[maxTRACK];
Double_t		M_Dca[maxTRACK];
//____________________________________________________________
TH1D *href[refMult_bin][Vz_bin][Psi_bin];

int centrality( int _refMult ){
    
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

bool refmult_check (int __nBTOFMatch, int __refMult){

//We do the Y-projection in a nBTOFMatch window, fit by double negative binomial distribution then get these parameters.

double a0=6.44225, a1=1.37398, a2=-0.00446098, a3=2.60578e-05, a4= -6.84022e-08, a5=6.18087e-11;

double refmultcutmax=a0+a1*(__nBTOFMatch)+a2*pow(__nBTOFMatch,2)+a3*pow(__nBTOFMatch,3)+a4*pow(__nBTOFMatch,4)+a5*pow(__nBTOFMatch,5);


        //cout<<"refmult,maxcut,mincut= "<<__refMult<<" "<<refmultcutmax<<" "<<refmultcutmin<<endl;

        if( __refMult<refmultcutmax){
                return true;
        }else{
                return false;
        }

}

bool ntrk_check(int __nBTOFMatch, int __ntrk){

    double b0=-14.7884, b1=1.07711, b2=-0.000186869, b3=1.20443e-06, b4=-1.745e-09;
    double c0=11.1393, c1=1.88776, c2=-0.0024491, c3=8.43354e-06, c4=-1.29845e-08;

    double ntrkcutmin=b0+b1*(__nBTOFMatch)+b2*pow(__nBTOFMatch,2)+b3*pow(__nBTOFMatch,3)+b4*pow(__nBTOFMatch,4);
    double ntrkcutmax=c0+c1*(__nBTOFMatch)+c2*pow(__nBTOFMatch,2)+c3*pow(__nBTOFMatch,3)+c4*pow(__nBTOFMatch,4);

    if( __ntrk<ntrkcutmax && __ntrk>ntrkcutmin){
                return true;
        }else{
                return false;
        }


}


int MultBin_Cen( int _NPt){
    
    float   MBins  [refMult_bin] = {470,390,320,231} ;
    Int_t   BinID     [refMult_bin] = { 3,2,1,0};
    
    Int_t   myBin =-1;

    for(int i=0;i!=refMult_bin;i++){
        if( _NPt > MBins[i] ){
            myBin = BinID[i];
            break;
        }
        else{
            myBin = -1;
        }
    }
    
    return myBin ;

}

int MultBin_Per( int _NPt){
    
    float   MBins  [refMult_bin] = {62,34,17,1} ;
    Int_t   BinID     [refMult_bin] = { 3,2,1,0};
    
    Int_t   myBin =-1;

    for(int i=0;i!=refMult_bin;i++){
        if( _NPt > MBins[i] ){
            myBin = BinID[i];
            break;
        }
        else{
            myBin = -1;
        }
    }
    
    return myBin ;

}

void InitHist(){

    //_____________mult bin for peripheral________________
    if(Cen_mode==2){
        mult_start_stop_delta[0] = 2.0;
        mult_start_stop_delta[1] = 110.0;
        mult_start_stop_delta[2] = 1.0;
    }
    if(Cen_mode==0){
        mult_start_stop_delta[0] = 0.0;
        mult_start_stop_delta[1] = 616.0;
        mult_start_stop_delta[2] = 1.0;
    }

 //__________information for class_______________________


    int stop_refM=refMult_bin;
	int stop_vz=Vz_bin;
    //int stop_psi=Psi_bin;
	int stop_psi;

    if(Psi_mode==1||Psi_mode==2) stop_psi=Psi_bin;
    if(Psi_mode==0) stop_psi=1;
    
	for(int i_refM=0;i_refM<stop_refM;i_refM++){
		for(int i_vz=0;i_vz<stop_vz;i_vz++){
			for(int i_psi=0;i_psi<stop_psi;i_psi++){
                //TString mixed_event_outfile_name = "/gpfs/mnt/gpfs01/star/pwg/yanghe/";
				TString mixed_event_outfile_name = "/star/data05/scratch/yanghe/ntrk_mixclass/";
                if(Cen_mode==1) mixed_event_outfile_name +="central/F_mixedevents_mult_";
                if(Cen_mode==2) mixed_event_outfile_name +="peripheral/F_mixedevents_mult_";

				mixed_event_outfile_name += i_refM;
				mixed_event_outfile_name += "_vz_";
				mixed_event_outfile_name += i_vz;

                mixed_event_outfile_name += "_Psi_";
				mixed_event_outfile_name += i_psi;

                mixed_event_outfile_name += "/";
                mixed_event_outfile_name += FileOutput;
				mixed_event_outfile_name += ".root";

				fout[i_refM][i_vz][i_psi] = new TFile(mixed_event_outfile_name.Data(),"RECREATE");

				cout << "ME output file" << mixed_event_outfile_name.Data() << " created" << endl;

				outEventTree[i_refM][i_vz][i_psi] = new TTree("MCTree","MCTree");
				outEventTree[i_refM][i_vz][i_psi]->Branch("M_refmult",&M_refmult,"M_refmult/D");
				outEventTree[i_refM][i_vz][i_psi]->Branch("M_Vz",&M_Vz,"M_Vz/D");
				outEventTree[i_refM][i_vz][i_psi]->Branch("M_Vx",&M_Vx,"M_Vx/D");
				outEventTree[i_refM][i_vz][i_psi]->Branch("M_Vy",&M_Vy,"M_Vy/D");
                outEventTree[i_refM][i_vz][i_psi]->Branch("M_Psi2",&M_Psi2,"M_Psi2/D");
				//outEventTree[i_refM][i_vz][i_psi]->Branch("M_TofMatch",&M_TofMatch,"M_TofMatch/D");
				//outTrackNtuple[i_refM][i_vz][i_psi] = new TNtuple("TrackTree","TrackTree","Px:Py:Pz:Pt:Eta:Phi:Charge:Dca");
                outEventTree[i_refM][i_vz][i_psi]->Branch("M_numTrk",&M_numTrk,"M_numTrk/I");

				outEventTree[i_refM][i_vz][i_psi]->Branch("M_Pt",&M_Pt,"M_Pt[M_numTrk]/D");
				outEventTree[i_refM][i_vz][i_psi]->Branch("M_Eta",&M_Eta,"M_Eta[M_numTrk]/D");
				outEventTree[i_refM][i_vz][i_psi]->Branch("M_Phi",&M_Phi,"M_Phi[M_numTrk]/D");
				outEventTree[i_refM][i_vz][i_psi]->Branch("M_Charge",&M_Charge,"M_Charge[M_numTrk]/D");
				outEventTree[i_refM][i_vz][i_psi]->Branch("M_Dca",&M_Dca,"M_Dca[M_numTrk]/D");

                //href[i_refM][i_vz][i_psi] = new TH1D("href","href",1000,0,1000);

			}
		}
	}




	
}


int main(int argc, char **argv){
 
    FileInput = argv[1];
    FileOutput = argv[2];
    Cen_mode =  atoi(argv[3]);
    Psi_mode = atoi(argv[4]);
    
    cout<<"We read in "<<FileInput<<endl;
    
    TChain *mychain = new TChain("EventTree");
    readtree pptree(mychain);

    //==================== read in list====================================
    int fileNumber = 0;
    char FileList[512];
    ifstream* inputStream = new ifstream;
    inputStream->open(FileInput);
    if (!(inputStream))
    {
        printf("can not open list file\n");
        return 0;
    }
    for (;inputStream->good();)
    {
        inputStream->getline(FileList,512);
        if  ( inputStream->good() )
        {
            TFile *ftmp = new TFile(FileList);
            if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys()))
            {
                printf(" file %s error in opening!!!\n",FileList);
            }
            else
            {
                printf(" read in file %s\n",FileList);
                mychain->Add(FileList);
                fileNumber++;
            }
            delete ftmp;
        }
    }
    printf(" files read in %d\n",fileNumber);

    //==============================================================================

    Long64_t nentries = mychain->GetEntries();
    cout<<"We have "<<nentries<<" events"<<endl;

    InitHist();

    Long64_t nbytes = 0, nb = 0;

    //nentries = 10000;
    for (Long64_t jentry=0; jentry<nentries;jentry++){
        nb = mychain->GetEntry(jentry);   nbytes += nb;
        if(jentry%1000==0) cout<<"working on "<<jentry<<endl;



        double zTpc = pptree.Vz;
        double refMultCorr = pptree.refmult;
        double j_Vx = pptree.Vx;
        double j_Vy = pptree.Vy;
        int _NTracks = pptree.numTrk;
        

        int CentID = centrality(refMultCorr);
        if(CentID<0 || CentID > NCENT ) continue;

        if(Cen_mode==1){ if(CentID!=0) continue;}
        if(Cen_mode==2){ if( !(CentID==6||CentID==7) ) continue;}
	    //if(refMultCorr < mult_start_stop_delta[0] || refMultCorr > mult_start_stop_delta[1]) continue;
        /*
        if(Cen_mode==1){
                if(refMultCorr < 255 || refMultCorr > 400) continue; }
        if(Cen_mode==2){
                if(refMultCorr < 15 || refMultCorr > 44) continue; }
        */
        //cout<<"mode "<<Cen_mode<<endl;
        //if(Psi_mode == 1) rawEPDPsi();

	if(zTpc>25) continue;
	//pile up rejection________________________________________
	
    int Ntofmatch = pptree.NBTOFMultfit;
	/*
    double NBTOFMult_fit;
	if(Ntofmatch < 14){ NBTOFMult_fit = 0.59207 + 2.1317*Ntofmatch; }
	else{ NBTOFMult_fit = 18.770 + 1.0699*Ntofmatch; }
	if(refMultCorr >= NBTOFMult_fit) continue;
    */
    if(!refmult_check(Ntofmatch,refMultCorr)) continue;

        //cal TPC psi________________________________________________________________________________
        if(Psi_mode == 2) {


            memset(QxRaw, 0, sizeof(QxRaw));
            memset(QyRaw, 0, sizeof(QyRaw));
            memset(QwRaw, 0, sizeof(QwRaw));

            for(int j=0; j <_NTracks; j++){
                double _pt = pptree.Pt[j];
                double _eta = pptree.Eta[j];
                double _phi = pptree.Phi[j];
                double _dca = pptree.Dca[j];
                double _charge = pptree.Charge[j];

                if(_pt>5.0) continue;
                if(TMath::Abs(_eta)>1.0) continue;
                if(TMath::Abs(_dca)>1.0) continue;
                double weight = (_pt<2.0)?_pt:2.0;

                for(int ih = 0; ih != 2; ih++){
            
            
                
                    QwRaw[ih]+= weight;
                    QxRaw[ih] += weight*cos((ih+2)*_phi);
                    QyRaw[ih] += weight*sin((ih+2)*_phi);
                
                }//ih

            }   // track loop

        
        }
        //tpc psi end________________________________

    
        double Psi2 = TMath::ATan2(QyRaw[0],QxRaw[0]);
        Psi2 /= 2.0;
        //cout<<"psi "<<Psi2<<endl;
        double evt_psi2 = 0;

        if(Psi_mode !=0){
            evt_psi2 = Psi2;
            if(Psi_mode == 1) evt_psi2 = pptree.EPD_Psi;

        }


        //cout<<"zbin "<<z_bin<<" mbin "<<mult_bin<<" psibin "<<psi_bin<<endl;   

        //initialize_____________________________
	    
        memset(M_Pt, 0, sizeof(M_Pt));
        memset(M_Eta, 0, sizeof(M_Eta));
        memset(M_Phi, 0, sizeof(M_Phi));
        memset(M_Charge, 0, sizeof(M_Charge));
        memset(M_Dca, 0, sizeof(M_Dca));

        //________________________________________________________

	    int NPt =0;


        for(int itrack=0; itrack<_NTracks; itrack++){
        
            double _pt = pptree.Pt[itrack];
            double _eta = pptree.Eta[itrack];
            double _phi = pptree.Phi[itrack];
            double _charge = pptree.Charge[itrack];
            double _dca = pptree.Dca[itrack];

            
	    if(TMath::Abs(_eta)>1) continue;
            if(TMath::Abs(_dca)>1) continue;
            if(_pt<0.2||_pt>30)  continue;
            

            M_Pt[NPt] = _pt;
		    M_Eta[NPt] = _eta;
		    M_Phi[NPt] = _phi;
		    M_Charge[NPt] = _charge;
		    M_Dca[NPt] = _dca;
            //cout<<"pt "<<M_Pt[NPt]<<endl;
		    NPt++;



        } //track loop
	

    /*
	if(Cen_mode==2){
		double CorNch = 0.45*NPt-4;
		if(Ntofmatch<CorNch) continue;
		double CorN = 0.75*NPt-21;
                if(Ntofmatch<CorN) continue;
		//if(NPt>110) continue;
	}       
	*/
    //________________tof vs nch pile up__________________

    if(!ntrk_check(Ntofmatch,NPt)) continue;

    //________________________________________________

        double totaltrack= NPt;

        //____select event bins____________
	    Float_t delta_z = (vertex_z_start_stop_delta[1] - vertex_z_start_stop_delta[0])/((Float_t)Vz_bin);
        vertex_z_start_stop_delta[2] = delta_z;

	    Float_t delta_mult = (mult_start_stop_delta[1] - mult_start_stop_delta[0])/((Float_t)refMult_bin);
        mult_start_stop_delta[2] = delta_mult;

        Float_t delta_Psi = (Psi_start_stop_delta[1] - Psi_start_stop_delta[0])/((Float_t)Psi_bin);
        Psi_start_stop_delta[2] = delta_Psi;

	    Int_t z_bin = -1;
        if(zTpc > vertex_z_start_stop_delta[0] && zTpc < vertex_z_start_stop_delta[1])
        {
            z_bin = (Int_t)((zTpc-vertex_z_start_stop_delta[0])/vertex_z_start_stop_delta[2]);
        }
        /*
        Int_t mult_bin = -1;
        if(totaltrack > mult_start_stop_delta[0] && totaltrack < mult_start_stop_delta[1])
        {
            mult_bin = (Int_t)((totaltrack-mult_start_stop_delta[0])/mult_start_stop_delta[2]);
        }
        */
        Int_t psi_bin = -1;
        if(Psi_mode==1||Psi_mode==2 ){
            if(evt_psi2 > Psi_start_stop_delta[0] && evt_psi2 < Psi_start_stop_delta[1])
            {
                psi_bin = (Int_t)((evt_psi2-Psi_start_stop_delta[0])/Psi_start_stop_delta[2]);
            }
        }
        if(Psi_mode==0) psi_bin = 0;

        Int_t mult_bin = -1;
        if(Cen_mode==1) mult_bin= MultBin_Cen(NPt);
        if(Cen_mode==2) mult_bin= MultBin_Per(NPt);

	
	    if(z_bin==-1||mult_bin==-1||psi_bin==-1) continue;
 
        //M_refmult         = 1;
        M_refmult         = refMultCorr;
	    M_Vz              = zTpc;
	    M_Vx			  = j_Vx;
	    M_Vy			  = j_Vy;
        M_numTrk          = NPt;
        M_Psi2            = evt_psi2;
        M_TofMatch	  = Ntofmatch;
        //href[mult_bin][z_bin][psi_bin]->Fill(refMultCorr);
        //cout<<"mpsi"<<M_Psi2<<endl;
        outEventTree[mult_bin][z_bin][psi_bin]->Fill();
     


    } //event loop

    int stop_refM=refMult_bin;
	int stop_vz=Vz_bin;
	//int stop_psi=Psi_bin;
    int stop_psi;

    if(Psi_mode==1||Psi_mode==2) stop_psi=Psi_bin;
    if(Psi_mode==0) stop_psi=1;

	for(int i_refM=0;i_refM<stop_refM;i_refM++){
		for(int i_vz=0;i_vz<stop_vz;i_vz++){
			for(int i_psi=0;i_psi<stop_psi;i_psi++){
				fout[i_refM][i_vz][i_psi]->cd();
				outEventTree[i_refM][i_vz][i_psi]->Write();
				//outTrackNtuple[i_refM][i_vz][i_psi]->Write();
                fout[i_refM][i_vz][i_psi]->Write();
				fout[i_refM][i_vz][i_psi]->Close();
			}
		}
	}



}
