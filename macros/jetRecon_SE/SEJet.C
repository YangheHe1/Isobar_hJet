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
#include <stdlib.h>
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
const Int_t maxTRACK=10000;
//_________________________________________________________
char filename[200];
char name[200];

char *FileInput=0;
char *FileInput1=0;
char *FileOutput=0;

ifstream inData;
ofstream outData;

//_________jet tree____________________________________________

Int_t       EventIndex;
Double_t    refmult;
Double_t    TrgEta;
Double_t    TrgPhi;
Double_t    TrgEt;
Int_t       TrgId;
Double_t    Rho;
Double_t    Sigma;
Int_t       PrimTrk;

Int_t       NJets;
Int_t       JetIndex[maxTRACK];
Int_t       JetNCons[maxTRACK];
Double_t    JetPt[maxTRACK];
Double_t    JetPtCorr[maxTRACK];
Double_t    JetEta[maxTRACK];
Double_t    JetPhi[maxTRACK];
Double_t    JetE[maxTRACK];
Double_t    JetArea[maxTRACK];


//____________________________________________________


TFile  *fout;
TTree *outTree;


void InitHist(){

    
    sprintf(name,"/star/data05/scratch/yanghe/jet_SE_Isobar/M4_R%.1lf_RootFile/%s.root",Jet_R,FileOutput);
    fout = new TFile(name,"recreate");

    outTree = new TTree("JetTree","JetTree");
    outTree->Branch("eventIndex",&EventIndex,"EventIndex/I");

    outTree->Branch("refmult",&refmult,"refmult/D");
    outTree->Branch("PrimTrk",&PrimTrk,"PrimTrk/I");
    outTree->Branch("TrgEta",&TrgEta,"TrgEta/D");
    outTree->Branch("TrgPhi",&TrgPhi,"TrgPhi/D");
    outTree->Branch("TrgEt",&TrgEt,"TrgEt/D");
    outTree->Branch("TrgId",&TrgId,"TrgId/I");
    outTree->Branch("Rho",&Rho,"Rho/D");
    outTree->Branch("Sigma",&Sigma,"Sigma/D");

    outTree->Branch("NJets",&NJets,"NJets/I");
    outTree->Branch("JetIndex",&JetIndex,"JetIndex[NJets]/I");
    outTree->Branch("JetNCons",&JetNCons,"JetNCons[NJets]/I");
    outTree->Branch("JetPt",&JetPt,"JetPt[NJets]/D");
    outTree->Branch("JetPtCorr",&JetPtCorr,"JetPtCorr[NJets]/D");
    outTree->Branch("JetEta",&JetEta,"JetEta[NJets]/D");
    outTree->Branch("JetPhi",&JetPhi,"JetPhi[NJets]/D");
    outTree->Branch("JetE",&JetE,"JetE[NJets]/D");
    outTree->Branch("JetArea",&JetArea,"JetArea[NJets]/D");

    

}


int main(int argc, char **argv){

    FileInput = argv[1];
    FileOutput = argv[2];
    Jet_R = atof(argv[3]);
    Jet_R_bkg = atof(argv[4]);

    
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

    Long64_t nbytes = 0, nb = 0;

    //________event level________________

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        nb = mychain->GetEntry(jentry);   nbytes += nb;
        if(jentry%1000==0) cout<<"working on "<<jentry<<endl;

        double zTpc = mtree.M_Vz;
        double refMultCorr = mtree.M_refmult;
        double j_Vx = mtree.M_Vx;
        double j_Vy = mtree.M_Vy;
        int _NTracks = mtree.M_numTrk;

        //if(refMultCorr<24||refMultCorr>400) continue;
        //if(refMultCorr<235&&refMultCorr>54) continue;
        
        
        //___________for jet___________________
        TVector3 Trigger[5];
        Int_t triggerindex=0;
        Int_t trigger=0;

		vector<PseudoJet> particles;
        particles.clear();
        //__________________________

        int Npt = 0;

        for(int itrack=0; itrack<_NTracks; itrack++){
            
            double _pt = mtree.M_Pt[itrack];
            double _eta = mtree.M_Eta[itrack];
            double _phi = mtree.M_Phi[itrack];
            double _charge = mtree.M_Charge[itrack];
            double _dca = mtree.M_Dca[itrack];

            if(TMath::Abs(_eta)>1) continue;
            if(TMath::Abs(_dca)>1) continue;
            if(_pt<0.2 || _pt>30) continue;

            TVector3 _v;
            _v.SetPtEtaPhi(_pt,_eta,_phi);
            double _px = _v.Px();
            double _py = _v.Py();
            double _pz = _v.Pz();

            //___________trigger_____________
            if(_pt>=7&&_pt<=30){
		        trigger +=1;
                triggerindex +=1;
                Trigger[triggerindex-1].SetPtEtaPhi(_pt,_eta,_phi);
                
            }
            

            if(_pt>0.2 && _pt<30){
                Double_t E = sqrt( (_px*_px)+(_py*_py)+(_pz*_pz)+(Mass_pi*Mass_pi));
                particles.push_back( PseudoJet(_px,_py,_pz,E));
            }

            Npt++;            

        } //____________track loop_________________

        //__________tigger selection_____________________
        double triggereta=0.;
	    double triggerphi=0.;
	    double triggerpt=0.;
        
	    
        if(triggerindex>0){
            if(triggerindex==1){
			    triggerpt=Trigger[0].Perp();
			    triggereta=Trigger[0].PseudoRapidity();
			    triggerphi=Trigger[0].Phi();
            }
		    else{
		        Int_t z=randomI.Integer(triggerindex);
                //Int_t z=gRandom->Uniform(0,triggerindex);
			    triggerpt=Trigger[z].Perp();
                triggereta=Trigger[z].PseudoRapidity();
                triggerphi=Trigger[z].Phi();
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

	    //EventIndex      = eventID;
	    refmult         = refMultCorr;
	    NJets           = no0f_jets;
	    TrgEta          = triggereta;
        TrgPhi          = triggerphi;
        TrgEt           = triggerpt;
        //TrgId           = triggerid;
        Rho             = jets_rho;
        Sigma           = jets_sigma;
        PrimTrk         = Npt;

        memset(JetIndex, 0, sizeof(JetIndex));
        memset(JetPt, 0, sizeof(JetPt));
        memset(JetPtCorr, 0, sizeof(JetPtCorr));
        memset(JetEta, 0, sizeof(JetEta));
        memset(JetPhi, 0, sizeof(JetPhi));
        memset(JetE, 0, sizeof(JetE));
        memset(JetArea, 0, sizeof(JetArea));

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

            JetIndex[nj] = nj;
            JetPt[nj] = jet_pT;
            JetPtCorr[nj] = jet_pT_sub;
            JetEta[nj] = jet_PsudoRap;
            JetPhi[nj] = jet_Phi;
            JetE[nj] = jet_E;
            JetArea[nj] = jet_Area;
            JetNCons[nj] = nCon;
            
        }  //_________jet tree information____________
        
        outTree->Fill();

    } //_________event loop_______________


    fout -> Write();



}

































