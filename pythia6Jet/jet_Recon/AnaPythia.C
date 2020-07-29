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

#include <algorithm>
#include <TRandom.h>
#include "readtree.C"

#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

using namespace std;
using namespace fastjet;

static const double Mass_pi = 0.140;
float Jet_R=0.2;
float Jet_R_bkg=0.4;

TRandom randomI;
TRandom randomI1;
//___________________________________________
char filename[200];
char name[200];

char *FileInput=0;
char *FileInput1=0;
char *FileOutput=0;

ifstream inData;
ofstream outData;

//_________jet tree_________________________________________

vector<int> JetIndex;
vector<int> JetNCons;
vector<double> JetPt;
vector<double> JetPtCorr;
vector<double> JetEta;
vector<double> JetPhi;
vector<double> JetE;
vector<double> JetArea;
Int_t       EventIndex;
Double_t    refmult;
Int_t       NJets;
Double_t    TrgEta;
Double_t    TrgPhi;
Double_t    TrgEt;
Int_t       TrgId;
Double_t    Rho;
Double_t    Sigma;
Int_t       PrimTrk;

//_________________________________________________

vector<int> TrackID;


TFile  *fout;
TTree *outTree;

TH1D *nutrigger;
TH1D *pthist;
TH1D *triggerhist;
TH1D *pi_trigger;
TH1D *etahist;
TH1D *phihist;
TH1D *gammapt;
TH2D *etavsphi;
TH1D *Hrefmult;
TH1D *NPrimaryTrack;
TH1D *NewEtahist;
TH1D *Nevt20;

void InitHist(){

    TrackID.clear();
    TrackID.push_back(211);
    TrackID.push_back(-211);
    TrackID.push_back(321);
    TrackID.push_back(-321);
    TrackID.push_back(2212);
    TrackID.push_back(-2212);
    //TrackID.push_back(111);

    sprintf(name,"/star/data05/scratch/yanghe/jet_Pythia6/R%.1lf_RootFile/%s.root",Jet_R,FileOutput);
    fout = new TFile(name,"recreate");

    outTree = new TTree("JetTree","JetTree");
    outTree->Branch("eventIndex",&EventIndex,"EventIndex/I");

    outTree->Branch("refmult",&refmult,"refmult/D");
    outTree->Branch("NJets",&NJets,"NJets/I");
    outTree->Branch("PrimTrk",&PrimTrk,"PrimTrk/I");
    outTree->Branch("TrgEta",&TrgEta,"TrgEta/D");
    outTree->Branch("TrgPhi",&TrgPhi,"TrgPhi/D");
    outTree->Branch("TrgEt",&TrgEt,"TrgEt/D");
    outTree->Branch("TrgId",&TrgId,"TrgId/I");
    outTree->Branch("Rho",&Rho,"Rho/D");
    outTree->Branch("Sigma",&Sigma,"Sigma/D");

    outTree->Branch("JetIndex",&JetIndex);
    outTree->Branch("JetNCons",&JetNCons);
    outTree->Branch("JetPt",&JetPt);
    outTree->Branch("JetPtCorr",&JetPtCorr);
    outTree->Branch("JetEta",&JetEta);
    outTree->Branch("JetPhi",&JetPhi);
    outTree->Branch("JetE",&JetE);
    outTree->Branch("JetArea",&JetArea);

    

}

int main(int argc, char **argv){

    FileInput = argv[1];
    FileOutput = argv[2];
    Jet_R = atof(argv[3]);
    Jet_R_bkg = atof(argv[4]);
    
    cout<<"We read in "<<FileInput<<endl;

    TChain *mychain = new TChain("PythiaTree");
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

    //________event level________________

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        nb = mychain->GetEntry(jentry);   nbytes += nb;
        if(jentry%1000==0) cout<<"working on "<<jentry<<endl;

        int RefMult = pptree.mParticles_;
        int eventID = pptree.mEventId;

	    

        //___________for jet___________________
        TVector3 Trigger[5];
        Int_t Trigger_Id[5];
        
        Int_t triggerindex=0;
        

		vector<PseudoJet> particles;
        particles.clear();
        //__________________________

        int _NTracks = pptree.mParticles_;
        int Npt = 0;

        for(int itrack=0; itrack<_NTracks; itrack++){

            double _E = pptree.mParticles_fE[itrack];
            double _px = pptree.mParticles_fPx[itrack];
            double _py = pptree.mParticles_fPy[itrack];
            double _pz = pptree.mParticles_fPz[itrack];
            int _pid =  pptree.mParticles_fPdgCode[itrack];
            
            TLorentzVector iparticle;
            iparticle.SetPxPyPzE(_px,_py,_pz,_E);

            double _eta = iparticle.Eta();
            double _phi = iparticle.Phi();
            double _pt = iparticle.Pt();
            //cout<<"id "<<_pid<<endl;

            if(TMath::Abs(_eta)>1) continue;
            if(_pt<0.2 || _pt>30) continue;

            //_________select charge________________
            bool validTrigger = false;
            for(int i=0;i<TrackID.size();i++){
                if(_pid==TrackID[i]) validTrigger = true;
            }

            if(!validTrigger) continue;

            //________________________________

            
            if(_pt>0.2 && _pt<30){
                //Double_t E = sqrt( (_px*_px)+(_py*_py)+(_pz*_pz)+(Mass_pi*Mass_pi));
                particles.push_back( PseudoJet(_px,_py,_pz,_E));
            }

            
                if(_pt>=7&&_pt<=30){
                    triggerindex +=1;
                    Trigger[triggerindex-1].SetPtEtaPhi(_pt,_eta,_phi);
                    Trigger_Id[triggerindex-1]= _pid;
                }
            
            Npt++;


        }  //track end  

        //__________tigger selection_____________________

        double triggereta=0.;
	    double triggerphi=0.;
	    double triggerpt=0.;
        int    triggerid=0;
	    

        if(triggerindex>0){
            if(triggerindex==1){
			    triggerpt=Trigger[0].Perp();
			    triggereta=Trigger[0].PseudoRapidity();
			    triggerphi=Trigger[0].Phi();
                triggerid = Trigger_Id[0];
                
            }
		    else{
		        Int_t z=randomI.Integer(triggerindex);
                //Int_t z=gRandom->Uniform(0,triggerindex);
			    triggerpt=Trigger[z].Perp();
                triggereta=Trigger[z].PseudoRapidity();
                triggerphi=Trigger[z].Phi();
                triggerid=Trigger_Id[z];
                
            }
	    }

        if(triggerindex==0) continue;

        //__________________jet reconstruction______________________
        JetDefinition jet_def(antikt_algorithm, Jet_R);
        Double_t ghost_maxrap=1.0;    //if particles go up to y=1.0
        AreaDefinition area_def(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));

        vector<PseudoJet> jets;
        // calculates jet areas
        ClusterSequenceArea cs_hard(particles, jet_def, area_def);
        //return a vector of jets sorted into decreasing transverse momentum
        double pt_Min = 0.2;
        vector<PseudoJet> jets_cs = sorted_by_pt(cs_hard.inclusive_jets(pt_Min));
        //select jet with eta<1-R
        Selector Fiducial_cut_selector = SelectorAbsEtaMax(1.0 - Jet_R);
        jets = Fiducial_cut_selector(jets_cs);
        //_____________________________________________________________________________
        //        

        //background substract_____________________________________________
        JetDefinition jet_def_bkgd(kt_algorithm, Jet_R_bkg);
        AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(ghost_maxrap,1,0.01));
        //select eta<1,not the hardest
        int Remove_N_hardest=2;
        Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(Remove_N_hardest));
        //Background estimation
        JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
        bkgd_estimator.set_particles(particles);

        //setting subtractor
        Subtractor subtractor(&bkgd_estimator);
        Double_t jets_rho = bkgd_estimator.rho();         //median of pt/area
        Double_t jets_sigma = bkgd_estimator.sigma();     //fluctuation

        Int_t no0f_jets = jets.size();

	    EventIndex      = eventID;
	    refmult         =RefMult;
	    NJets           = no0f_jets;
	    TrgEta          = triggereta;
        TrgPhi          = triggerphi;
        TrgEt           = triggerpt;
        TrgId           = triggerid;
        Rho             = jets_rho;
        Sigma           = jets_sigma;
        PrimTrk         = Npt;

        JetIndex.clear();
        JetPt.clear();
        JetPtCorr.clear();
        JetEta.clear();
        JetPhi.clear();
        JetE.clear();
        JetArea.clear();

        for (Int_t nj = 0; nj < no0f_jets; nj++){
            //jet information
            double jet_pT= jets[nj].perp();
            double jet_Area = jets[nj].area();
            double jet_pT_sub= jets[nj].perp() - jet_Area*jets_rho;
            double jet_PsudoRap = jets[nj].pseudorapidity();
            double jet_Phi = jets[nj].phi_std();
            double jet_E =jets[nj].e();
            vector<PseudoJet> constituents = jets[nj].constituents();
            Int_t nCon = constituents.size();

            JetIndex.push_back(nj);
            JetPt.push_back(jet_pT);
            JetPtCorr.push_back(jet_pT_sub);
            JetEta.push_back(jet_PsudoRap);
            JetPhi.push_back(jet_Phi);
            JetE.push_back(jet_E);
            JetArea.push_back(jet_Area);
            JetNCons.push_back(nCon);
            
        }  //_________jet tree information____________

        outTree->Fill();

    } //event end

    fout -> Write();
}
