#include "TROOT.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TMath.h"
#include "TFile.h"
#include "TH2.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TComplex.h"
#include "TClonesArray.h"
#include <iostream>
#include <stdlib.h>

#include <algorithm>
#include <TRandom.h>
#include "readtree.C"

#include "mcTrackParticle.h"
#include "mcTrackEvent.h"


using std::vector;
using std::string;
using namespace std;


typedef vector< vector<Int_t> >  Two_dim_Int_vector;
static vector<mcTrackEvent*> MCTrackEvent_ME;
static mcTrackParticle  *MCTrackParticle_ME;
static TRandom ran;

char filename[200];
char name[200];

char *FileInput=0;
char *FileInput1=0;
char *FileOutput=0;

ifstream inData;
ofstream outData;

const Int_t maxTRACK=10000;
//_________ME tree________________
Double_t        ME_refmult;
Int_t           ME_numTrk;
Double_t        ME_Vz;
Double_t        ME_Vx;
Double_t        ME_Vy;
Double_t        ME_Psi2;
Double_t  ME_Pt[maxTRACK];
Double_t  ME_Eta[maxTRACK];
Double_t  ME_Phi[maxTRACK];
Double_t  ME_Charge[maxTRACK];
Double_t  ME_Dca[maxTRACK];


TFile  *fout;
TTree *outTree;
TH2D *h_tracks_vs_z_vertex;
TH2D *h_tracks_vs_z_vertex_Fill;
TH2D *h_ref_vs_z_vertex;
TH2D *h_ref_vs_z_vertex_Fill;
TH1D *Hme_number;
TH1D *Hnum_tracks;


//_________buffer information settings____________________________
Long64_t stop_event_use = 830000;

Int_t N_max_events; //max refmult

//__________test for min_ref__________
const Int_t N_max_events1 = 100;
Double_t Buffer_nTrk[N_max_events1];


Int_t Cen_mode;


void InitHist(){

    if(Cen_mode==1) sprintf(name,"/star/data05/scratch/yanghe/MixEvent/ntrk_M4_central/%s.root",FileOutput);
    if(Cen_mode==2) sprintf(name,"/star/data05/scratch/yanghe/MixEvent/ntrk_M4_peripheral/%s.root",FileOutput);
    fout = new TFile(name,"recreate");

    outTree = new TTree("METree","METree");
    outTree->Branch("ME_refmult",&ME_refmult,"ME_refmult/D");
    outTree->Branch("ME_numTrk",&ME_numTrk,"ME_numTrk/I");
    outTree->Branch("ME_Vz",&ME_Vz,"ME_Vz/D");
    outTree->Branch("ME_Vx",&ME_Vx,"ME_Vx/D");
    outTree->Branch("ME_Vy",&ME_Vy,"ME_Vy/D");
    outTree->Branch("ME_Psi2",&ME_Vy,"ME_Psi2/D");

    outTree->Branch("ME_Pt",&ME_Pt,"ME_Pt[ME_numTrk]/D");
    outTree->Branch("ME_Eta",&ME_Eta,"ME_Eta[ME_numTrk]/D");
    outTree->Branch("ME_Phi",&ME_Phi,"ME_Phi[ME_numTrk]/D");
    outTree->Branch("ME_Charge",&ME_Charge,"ME_Charge[ME_numTrk]/D");
    outTree->Branch("ME_Dca",&ME_Dca,"ME_Dca[ME_numTrk]/D");
/*
    TemTree = new TTree("temTree","temTree");
    TemTree->Branch("tem_refmult",&tem_refmult,"tem_refmult/D");
    TemTree->Branch("tem_Vz",&tem_Vz,"tem_Vz/D");
    TemTree->Branch("tem_Vx",&tem_Vx,"tem_Vx/D");
    TemTree->Branch("tem_Vy",&tem_Vy,"tem_Vy/D");
    TemTree->Branch("tem_Pt",&tem_Pt);
    TemTree->Branch("tem_Eta",&tem_Eta);
    TemTree->Branch("tem_Phi",&tem_Phi);
    TemTree->Branch("tem_Charge",&tem_Charge);
    TemTree->Branch("tem_Dca",$tem_Dca);
*/
    h_tracks_vs_z_vertex = new TH2D("h_tracks_vs_z_vertex","h_tracks_vs_z_vertex",100,-50,50,800,0,800);
    h_ref_vs_z_vertex = new TH2D("h_ref_vs_z_vertex","h_ref_vs_z_vertex",100,-50,50,800,0,800);

    h_tracks_vs_z_vertex_Fill = new TH2D("h_tracks_vs_z_vertex_Fill","h_tracks_vs_z_vertex_Fill",100,-50,50,800,0,800);
    h_ref_vs_z_vertex_Fill = new TH2D("h_ref_vs_z_vertex_Fill","h_ref_vs_z_vertex",100,-50,50,800,0,800);

    Hme_number = new TH1D("Hme_number","Hme_number",600,-0.5,599.5);

    Hnum_tracks = new TH1D("Hnum_tracks","Hnum_tracks",1000,0,1000);

}    

int main(int argc, char **argv){

    FileInput = argv[1];
    FileOutput = argv[2];
    Cen_mode = atoi(argv[3]);
    
    cout<<"We read in "<<FileInput<<endl;
    
    TChain *mychain = new TChain("MCTree");
    readtree mtree(mychain);

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

    

    //_____________set end number_________________
    Long64_t stop_event_use_loop = stop_event_use;
    if(stop_event_use_loop > nentries) stop_event_use_loop = nentries;

    Int_t ME_event_counter = 0;

    //_________set refmult pool____________
    Long64_t i_nbytes = 0, i_nb = 0;
    for (Long64_t i=0; i<nentries;i++){
        i_nb = mychain->GetEntry(i);    i_nbytes += i_nb;
        double i_RefMult = mtree.M_refmult;
        double i_vz = mtree.M_Vz;

        int i_NTracks =mtree.M_numTrk;

        //cout<<"NTracks"<<i_NTracks<<endl;
        Hnum_tracks->Fill(mtree.M_numTrk);
        h_tracks_vs_z_vertex->Fill(i_vz,i_NTracks);
        h_ref_vs_z_vertex->Fill(i_vz,i_RefMult);
    }

    //_________________initialize buffer_______________
    int NLastBin = Hnum_tracks->FindLastBinAbove(0,1);
    N_max_events=(Int_t)Hnum_tracks->GetBinCenter(NLastBin);
    N_max_events = N_max_events+1;
    //N_max_events=(Int_t)Hnum_tracks->FindLastBinAbove(0,1);

    Two_dim_Int_vector ME_track_number_vector; // stores the track ids for every event which is mixed
    ME_track_number_vector.resize(N_max_events);

    MCTrackEvent_ME.resize(N_max_events);

   for(Int_t mix_loop = 0; mix_loop < N_max_events; mix_loop++)
    {
        MCTrackEvent_ME[mix_loop] = new mcTrackEvent();
        //mcTrackParticle MCTrackEvent_ME[mix_loop];
    }

    int loopNUM=0;
    int meNUM=0;


    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        nb = mychain->GetEntry(jentry);   nbytes += nb;
        if(jentry%1000==0) cout<<"working on "<<jentry<<endl;

        double RefMult = mtree.M_refmult;
        double vz = mtree.M_Vz;
        double vx = mtree.M_Vx;
        double vy = mtree.M_Vy;
        double _psi2 = mtree.M_Psi2;

        int _NTracks = mtree.M_numTrk;

        //________set track id per event________________
        ME_track_number_vector[ME_event_counter].resize(_NTracks);

        for(Int_t i_track = 0; i_track < _NTracks; i_track++)
        {
            ME_track_number_vector[ME_event_counter][i_track] = i_track;
        }
        

        //store buffer event information_________________
        MCTrackEvent_ME[ME_event_counter]->clearParticleList();
        MCTrackEvent_ME[ME_event_counter]->setx(vx);
        MCTrackEvent_ME[ME_event_counter]->sety(vy);
        MCTrackEvent_ME[ME_event_counter]->setz(vz);
        MCTrackEvent_ME[ME_event_counter]->setmult(RefMult);
        MCTrackEvent_ME[ME_event_counter]->setpsi(_psi2);
        MCTrackEvent_ME[ME_event_counter]->setnumtrack(_NTracks);


        for(int itrack=0; itrack<_NTracks; itrack++){
            double pt = mtree.M_Pt[itrack];
            double eta = mtree.M_Eta[itrack];
            double phi = mtree.M_Phi[itrack];
            double charge = mtree.M_Charge[itrack];
            double dca = mtree.M_Dca[itrack];

            

            //store track information
            MCTrackParticle_ME = MCTrackEvent_ME[ME_event_counter]->createParticle();
            MCTrackParticle_ME ->set_dca_to_prim(dca);
            MCTrackParticle_ME ->set_Particle_charge(charge);
            MCTrackParticle_ME ->set_Particle_pt(pt);
            MCTrackParticle_ME ->set_Particle_eta(eta);
            MCTrackParticle_ME ->set_Particle_phi(phi);

        }

        /*
        double r1=MCTrackEvent_ME[ME_event_counter]->getmult();
        cout<<"Ref "<< RefMult<<endl;
        cout<<"refmult "<<r1<<endl;
        */

        ME_event_counter++;

        //test
        //for(int i=0;i<N_max_events1;i++){Buffer_nTrk[i]=9999;}


        if(ME_event_counter == N_max_events) // event buffer is full, start to create mixed events
        {
            for(Int_t mix_loop = 0; mix_loop < N_max_events; mix_loop++)
            {

                //get sample refmult/vz/Ntrack
                Double_t N_tracks_sample_d, N_z_vertex_sample, N_ref_sample_d;
                //h_tracks_vs_z_vertex ->GetRandom2(N_z_vertex_sample,N_tracks_sample_d);
                //Int_t N_tracks_sample = (Int_t)N_tracks_sample_d;
		N_tracks_sample_d = Hnum_tracks ->GetRandom();
		Int_t N_tracks_sample = (Int_t)N_tracks_sample_d;

                //cout<<"N_tracks_sample "<<N_tracks_sample<<endl;

                if(N_tracks_sample > N_max_events)
                {
                    cout << "Error: N_tracks_sample > N_max_events" << endl;
                    break;
                }

                Int_t Fill_flag = 1; // if all events in the loop have some entries left over then fill the event later
                //Double_t Et_total = 0.0;

                memset(ME_Pt, 0, sizeof(ME_Pt));
                memset(ME_Eta, 0, sizeof(ME_Eta));
                memset(ME_Phi, 0, sizeof(ME_Phi));
                memset(ME_Charge, 0, sizeof(ME_Charge));
                memset(ME_Dca, 0, sizeof(ME_Dca));


                for(Int_t ME_loop_counter = 0; ME_loop_counter < N_tracks_sample; ME_loop_counter++)
                {
                    Double_t prim_vertex_x   = MCTrackEvent_ME[ME_loop_counter]->getx();
                    Double_t prim_vertex_y   = MCTrackEvent_ME[ME_loop_counter]->gety();
                    Double_t prim_vertex_z   = MCTrackEvent_ME[ME_loop_counter]->getz();
                    Double_t refMult         = MCTrackEvent_ME[ME_loop_counter]->getmult();
                    Double_t Psi             = MCTrackEvent_ME[ME_loop_counter]->getpsi();

                    Int_t   N_Particles     = MCTrackEvent_ME[ME_loop_counter]->getNumParticle();
                    //Double_t nn     = MCTrackEvent_ME[ME_loop_counter]->getNumParticle();

                    if((N_Particles-mix_loop) < 0) // at least one event in the loop is out of tracks, stop and don't fill the tree
                    {
                        Fill_flag = 0;
                        break;
                    }

                    //cout<<"refmult "<<refMult<<endl;
                    //test buffer efficiency
                    //if(mix_loop==0)Buffer_nTrk[ME_loop_counter] = nn;

                    //set event level information for one ME
                    if(ME_loop_counter == 0){
                        ME_refmult  = refMult;
                        ME_Vx       = prim_vertex_x;
                        ME_Vy       = prim_vertex_y;
                        ME_Vz       = prim_vertex_z;
                        ME_numTrk   = N_tracks_sample;
                        ME_Psi2     = Psi;
                        

                    }
                    
                    //select random particles
                    Int_t track_num_random    = (Int_t)ran.Integer(N_Particles-mix_loop);
                    Int_t track_id_num_random = ME_track_number_vector[ME_loop_counter][track_num_random];
                    ME_track_number_vector[ME_loop_counter][track_num_random] = ME_track_number_vector[ME_loop_counter][N_Particles-mix_loop-1];

                    //pick each track information
                    MCTrackParticle_ME         = MCTrackEvent_ME[ME_loop_counter]->getParticle(track_id_num_random);
                    Double_t mdca                 = MCTrackParticle_ME->get_dca_to_prim();
                    Double_t mpt                  = MCTrackParticle_ME->get_Particle_pt();
                    Double_t meta                  = MCTrackParticle_ME->get_Particle_eta();
                    Double_t mphi                  = MCTrackParticle_ME->get_Particle_phi();
                    Double_t mcharge               = MCTrackParticle_ME->get_Particle_charge();

                    //___Fill ME track tree_________________
                    ME_Pt[ME_loop_counter] = mpt;
                    ME_Eta[ME_loop_counter] = meta;
                    ME_Phi[ME_loop_counter] = mphi;
                    ME_Dca[ME_loop_counter] = mdca;
                    ME_Charge[ME_loop_counter] = mcharge;

                    //cout<<"pt "<<mpt<<endl;

                }   //pick track from different SE

                if(Fill_flag == 1) {
                    outTree->Fill();
                    //h_tracks_vs_z_vertex_Fill ->Fill(N_z_vertex_sample,N_tracks_sample_d);
                    //cout<<"me number  "<<meNUM<<endl;
                    //meNUM++;

                }

            }   //pick ME from SE pool(N=max_events)

            //double real_loopnum_buffer = *min_element(Buffer_nTrk,Buffer_nTrk+N_max_events1);
            //if(real_loopnum_buffer!=9999)Hme_number->Fill(real_loopnum_buffer);

            ME_event_counter = 0;
        } //buffer

        

    } //load event

    cout<<"loop number "<<loopNUM<<"  me number  "<<meNUM<<endl;
    fout -> Write();


}   
