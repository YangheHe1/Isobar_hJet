#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include<stdio.h>
#include<algorithm>
#include<vector>
#include<iostream>
#include <string>
#include <stdlib.h> // strtod(), string to double
#include "drawheader_v4.h"
#include "TGraph.h"

using std::vector;
using std::string;

using namespace std;

const int NBINS = 15;
const int NRT   = 4;
const float Xlow  = 1;
const float Xhigh = 500;
const double PI = acos(-1.0);




TFile *fin;
TFile *fin_mix;
TFile *fout;
TLine* line;

// settings____________________________________
double Jet_R=0.2;



//________0 R0.2 1 R0.3 2 R0.4_______________ 
TH1D *hJetPt_SE_trg7_C[3];
TH1D *hNtrigger_SE_trg7_C[3];
TH1D *hJetPt_ME_trg7_C[3];
TH1D *hNtrigger_ME_trg7_C[3];
TH1D *hJetPt_ME_scaled_trg7_C[3];
TH1D *hJet_ME_Sub_SE_trg7_C[3];

TH1D *hRatio_trg7_C[3];
TH1D *hRatio_trg7_scaled_C[3];


TH1D *hJetPt_SE_trg7_P[3];
TH1D *hNtrigger_SE_trg7_P[3];
TH1D *hJetPt_ME_trg7_P[3];
TH1D *hNtrigger_ME_trg7_P[3];
TH1D *hJetPt_ME_scaled_trg7_P[3];
TH1D *hJet_ME_Sub_SE_trg7_P[3];

TH1D *hRatio_trg7_P[3];
TH1D *hRatio_trg7_scaled_P[3];


TH1D *hJetPt_SE_trg9_C[3];
TH1D *hNtrigger_SE_trg9_C[3];
TH1D *hJetPt_ME_trg9_C[3];
TH1D *hNtrigger_ME_trg9_C[3];
TH1D *hJetPt_ME_scaled_trg9_C[3];
TH1D *hJet_ME_Sub_SE_trg9_C[3];

TH1D *hRatio_trg9_C[3];
TH1D *hRatio_trg9_scaled_C[3];

TH1D *hJetPt_SE_trg9_P[3];
TH1D *hNtrigger_SE_trg9_P[3];
TH1D *hJetPt_ME_trg9_P[3];
TH1D *hNtrigger_ME_trg9_P[3];
TH1D *hJetPt_ME_scaled_trg9_P[3];
TH1D *hJet_ME_Sub_SE_trg9_P[3];

TH1D *hRatio_trg9_P[3];
TH1D *hRatio_trg9_scaled_P[3];


double scalingfactor_trg9_C[3];
double scalingfactor_trg9_P[3];
double scalingfactor_trg7_C[3];
double scalingfactor_trg7_P[3];

//__________________________________________________________
TLegend *leg;

TGraphErrors *grtmpErr0[100];
TGraph *gr[100];
char name[200];

TCanvas *can[100];
TPad* pad[20][20];

TH1D *hhtem[100];
TH2D *h2D[100];

double llmrg[] = {0.2, 0.02, 0.2, 0.2};
double Lmrg[]  = {0.2,0.2};

double ratx[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
double raty[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

int mrk[] = {25,24,21,25, 33, 27};
int clr[] = {kBlue,kRed+1, kBlue, 1,  kGreen+1,kGray+2, kCyan+1, kYellow+1,kViolet+1};

const char* type[]={"same event", "mix event","mix event(scaled)","norm.region", "^{4}He+^{14}N 1 PeV", "^{4}He+^{14}N 10 PeV"};

void readin(){
    
    float nraw_trg7_P[3];
    float nraw_trg7_C[3];
    float nraw_trg9_P[3];
    float nraw_trg9_C[3];

    float nmix_trg7_P[3];
    float nmix_trg7_C[3];
    float nmix_trg9_P[3];
    float nmix_trg9_C[3];


    for(int i=0;i<3;i++){

    int j=i+2;
    double jet_radius= 0.1*j;
    //===============================  read in =====================================
    sprintf(name,"R%i/jet_SE_R%i.root",j,j);
    cout<<name<<endl;
    fin = TFile::Open(name);

    //trg9-30
    //central
    //float nraw_trg9_C[3];
    
    sprintf(name,"recoil_jet_trg9_30_cent0");
    hJetPt_SE_trg9_C[i]     = (TH1D*)fin->Get(name);

    sprintf(name,"number_of_trigger9_30_cent0");
    hNtrigger_SE_trg9_C[i]    = (TH1D*)fin->Get(name);

    nraw_trg9_C[i] = hNtrigger_SE_trg9_C[i] -> GetEntries();
    hJetPt_SE_trg9_C[i] -> Scale(1./nraw_trg9_C[i]);

    //peripheral
    //float nraw_trg9_P[3];
    
    sprintf(name,"recoil_jet_trg9_30_cent1");
    hJetPt_SE_trg9_P[i]     = (TH1D*)fin->Get(name);

    sprintf(name,"number_of_trigger9_30_cent1");
    hNtrigger_SE_trg9_P[i]    = (TH1D*)fin->Get(name);

    nraw_trg9_P[i] = hNtrigger_SE_trg9_P[i] -> GetEntries();
    hJetPt_SE_trg9_P[i] -> Scale(1./nraw_trg9_P[i]);
    //___________________________________________

    //trg7-30
    //central
    //float nraw_trg7_C[3];
    
    sprintf(name,"recoil_jet_cent0");
    hJetPt_SE_trg7_C[i]     = (TH1D*)fin->Get(name);

    sprintf(name,"number_of_trigger_cent0");
    hNtrigger_SE_trg7_C[i]    = (TH1D*)fin->Get(name);

    nraw_trg7_C[i] = hNtrigger_SE_trg7_C[i] -> GetEntries();
    hJetPt_SE_trg7_C[i] -> Scale(1./nraw_trg7_C[i]);

    //peripheral
    //float nraw_trg7_P[3];
    
    sprintf(name,"recoil_jet_cent1");
    hJetPt_SE_trg7_P[i]     = (TH1D*)fin->Get(name);

    sprintf(name,"number_of_trigger_cent1");
    hNtrigger_SE_trg7_P[i]    = (TH1D*)fin->Get(name);

    nraw_trg7_P[i] = hNtrigger_SE_trg7_P[i] -> GetEntries();
    hJetPt_SE_trg7_P[i] -> Scale(1./nraw_trg7_P[i]);

    
    //===============================  read in =====================================
    sprintf(name,"R%i/jet_ME_R%i.root",j,j);
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    //trg9-30
    //central
    
    sprintf(name,"recoil_jet_cent0");
    hJetPt_ME_trg9_C[i]     = (TH1D*)fin->Get(name);

    sprintf(name,"number_of_trigger_cent0");
    hNtrigger_ME_trg9_C[i]    = (TH1D*)fin->Get(name);

    nmix_trg9_C[i] = hNtrigger_ME_trg9_C[i] -> GetEntries();
    hJetPt_ME_trg9_C[i] -> Scale(1./nmix_trg9_C[i]);

    //peripheral
    //float nraw_trg9_P[3];
    
    sprintf(name,"recoil_jet_cent1");
    hJetPt_ME_trg9_P[i]     = (TH1D*)fin->Get(name);

    sprintf(name,"number_of_trigger_cent1");
    hNtrigger_ME_trg9_P[i]    = (TH1D*)fin->Get(name);

    nmix_trg9_P[i] = hNtrigger_ME_trg9_P[i] -> GetEntries();
    hJetPt_ME_trg9_P[i] -> Scale(1./nmix_trg9_P[i]);
    //___________________________________________

    sprintf(name,"R%i/jet_ME_R%i.root",j,j);
    cout<<name<<endl;
    fin = TFile::Open(name);
    //trg7-30
    //central
    //float nraw_trg7_C[3];
    
    sprintf(name,"recoil_jet_cent0");
    hJetPt_ME_trg7_C[i]     = (TH1D*)fin->Get(name);

    sprintf(name,"number_of_trigger_cent0");
    hNtrigger_ME_trg7_C[i]    = (TH1D*)fin->Get(name);

    nmix_trg7_C[i] = hNtrigger_ME_trg7_C[i] -> GetEntries();
    hJetPt_ME_trg7_C[i] -> Scale(1./nmix_trg7_C[i]);

    //peripheral
    //float nraw_trg7_P[3];
    
    sprintf(name,"recoil_jet_cent1");
    hJetPt_ME_trg7_P[i]     = (TH1D*)fin->Get(name);

    sprintf(name,"number_of_trigger_cent1");
    hNtrigger_ME_trg7_P[i]    = (TH1D*)fin->Get(name);

    nmix_trg7_P[i] = hNtrigger_ME_trg7_P[i] -> GetEntries();
    hJetPt_ME_trg7_P[i] -> Scale(1./nmix_trg7_P[i]);

    
    //_____________________________________________________________________

    //trg9-30
    //central

            Int_t binmax=hJetPt_SE_trg9_C[i]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_SE_trg9_C[i]->GetBinContent(ibin);
                double deltap = hJetPt_SE_trg9_C[i]->GetBinWidth(ibin);
                double deltaeta = 2-(2*jet_radius);
                double binerr = hJetPt_SE_trg9_C[i]->GetBinError(ibin);
                hJetPt_SE_trg9_C[i]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_SE_trg9_C[i]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }

    //peripheral
            binmax=hJetPt_SE_trg9_P[i]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_SE_trg9_P[i]->GetBinContent(ibin);
                double deltap = hJetPt_SE_trg9_P[i]->GetBinWidth(ibin);
                double deltaeta = 2-(2*jet_radius);
                double binerr = hJetPt_SE_trg9_P[i]->GetBinError(ibin);
                hJetPt_SE_trg9_P[i]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_SE_trg9_P[i]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }

    //trg7-30
    //central

            binmax=hJetPt_SE_trg7_C[i]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_SE_trg7_C[i]->GetBinContent(ibin);
                double deltap = hJetPt_SE_trg7_C[i]->GetBinWidth(ibin);
                double deltaeta = 2-(2*jet_radius);
                double binerr = hJetPt_SE_trg7_C[i]->GetBinError(ibin);
                hJetPt_SE_trg7_C[i]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_SE_trg7_C[i]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }

    //peripheral
            binmax=hJetPt_SE_trg7_P[i]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_SE_trg7_P[i]->GetBinContent(ibin);
                double deltap = hJetPt_SE_trg7_P[i]->GetBinWidth(ibin);
                double deltaeta = 2-(2*jet_radius);
                double binerr = hJetPt_SE_trg7_P[i]->GetBinError(ibin);
                hJetPt_SE_trg7_P[i]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_SE_trg7_P[i]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }
    


    

    //ME——————————————————————————————————————————————
    //trg9-30
    //central

            Int_t binmax=hJetPt_ME_trg9_C[i]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_ME_trg9_C[i]->GetBinContent(ibin);
                double deltap = hJetPt_ME_trg9_C[i]->GetBinWidth(ibin);
                double deltaeta = 2-(2*jet_radius);
                double binerr = hJetPt_ME_trg9_C[i]->GetBinError(ibin);
                hJetPt_ME_trg9_C[i]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_ME_trg9_C[i]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }

    //peripheral
            binmax=hJetPt_ME_trg9_P[i]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_ME_trg9_P[i]->GetBinContent(ibin);
                double deltap = hJetPt_ME_trg9_P[i]->GetBinWidth(ibin);
                double deltaeta = 2-(2*jet_radius);
                double binerr = hJetPt_ME_trg9_P[i]->GetBinError(ibin);
                hJetPt_ME_trg9_P[i]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_ME_trg9_P[i]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }

    //trg7-30
    //central

            binmax=hJetPt_ME_trg7_C[i]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_ME_trg7_C[i]->GetBinContent(ibin);
                double deltap = hJetPt_ME_trg7_C[i]->GetBinWidth(ibin);
                double deltaeta = 2-(2*jet_radius);
                double binerr = hJetPt_ME_trg7_C[i]->GetBinError(ibin);
                hJetPt_ME_trg7_C[i]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_ME_trg7_C[i]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }

    //peripheral
            binmax=hJetPt_ME_trg7_P[i]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_ME_trg7_P[i]->GetBinContent(ibin);
                double deltap = hJetPt_ME_trg7_P[i]->GetBinWidth(ibin);
                double deltaeta = 2-(2*jet_radius);
                double binerr = hJetPt_ME_trg7_P[i]->GetBinError(ibin);
                hJetPt_ME_trg7_P[i]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_ME_trg7_P[i]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }
    

    //ratio calculation_______________________________________
    //trg9
    //centarl
        hRatio_trg9_C[i] = (TH1D*) hJetPt_SE_trg9_C[i]->Clone();
        binmax=hJetPt_SE_trg9_C[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg9_C[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_trg9_C[i]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio_trg9_C[i]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE_trg9_C[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_trg9_C[i]->GetBinError(ibin);
            hRatio_trg9_C[i]->SetBinContent(ibin,y/yref);
            double err = (binerr/yref)-(y*binerref/(yref*yref));
            hRatio_trg9_C[i]->SetBinError(ibin,err);

        }

    //peripheral
        hRatio_trg9_P[i] = (TH1D*) hJetPt_SE_trg9_P[i]->Clone();
        binmax=hJetPt_SE_trg9_P[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg9_P[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_trg9_P[i]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio_trg9_P[i]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE_trg9_P[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_trg9_P[i]->GetBinError(ibin);
            hRatio_trg9_P[i]->SetBinContent(ibin,y/yref);
            double err = (binerr/yref)-(y*binerref/(yref*yref));
            hRatio_trg9_P[i]->SetBinError(ibin,err);

        }

    //trg7
    //centarl
        hRatio_trg7_C[i] = (TH1D*) hJetPt_SE_trg7_C[i]->Clone();
        binmax=hJetPt_SE_trg7_C[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg7_C[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_trg7_C[i]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio_trg7_C[i]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE_trg7_C[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_trg7_C[i]->GetBinError(ibin);
            hRatio_trg7_C[i]->SetBinContent(ibin,y/yref);
            double err = (binerr/yref)-(y*binerref/(yref*yref));
            hRatio_trg7_C[i]->SetBinError(ibin,err);

        }

    //peripheral
        hRatio_trg7_P[i] = (TH1D*) hJetPt_SE_trg7_P[i]->Clone();
        binmax=hJetPt_SE_trg7_P[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg7_P[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_trg7_P[i]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio_trg7_P[i]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE_trg7_P[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_trg7_P[i]->GetBinError(ibin);
            hRatio_trg7_P[i]->SetBinContent(ibin,y/yref);
            double err = (binerr/yref)-(y*binerref/(yref*yref));
            hRatio_trg7_P[i]->SetBinError(ibin,err);

        }



    }

    //read in scaling factor
    ifstream f_T7("scalingfactor_T7.txt");
    for(int i=0;i<3;i++){
        f_T7>>scalingfactor_trg7_C[i]>>scalingfactor_trg7_P[i];
        cout<<scalingfactor_trg7_C[i]<<"  "<<scalingfactor_trg7_P[i]<<endl;
    }  
    f_T7.close();     

    ifstream f_T9("scalingfactor_T9.txt");
    for(int i=0;i<3;i++){
        f_T9>>scalingfactor_trg9_C[i]>>scalingfactor_trg9_P[i];
        cout<<scalingfactor_trg9_C[i]<<"  "<<scalingfactor_trg9_P[i]<<endl;
    }
    f_T9.close();
    

    
    
   
    //_____________ME sclaing_____________________
    for(int i=0;i<3;i++){

        

	    hJetPt_ME_scaled_trg9_C[i] = (TH1D*) hJetPt_ME_trg9_C[i]->Clone(); 
        hJetPt_ME_scaled_trg9_C[i]->Scale(scalingfactor_trg9_C[i]);

        hJetPt_ME_scaled_trg9_P[i] = (TH1D*) hJetPt_ME_trg9_P[i]->Clone(); 
        hJetPt_ME_scaled_trg9_P[i]->Scale(scalingfactor_trg9_P[i]);

        hJetPt_ME_scaled_trg7_C[i] = (TH1D*) hJetPt_ME_trg7_C[i]->Clone(); 
        hJetPt_ME_scaled_trg7_C[i]->Scale(scalingfactor_trg7_C[i]);

        hJetPt_ME_scaled_trg7_P[i] = (TH1D*) hJetPt_ME_trg7_P[i]->Clone(); 
        hJetPt_ME_scaled_trg7_P[i]->Scale(scalingfactor_trg7_P[i]);

        
    //SCALED ratio calculation_______________________________________
    //trg9
    //centarl
        hRatio_trg9_scaled_C[i] = (TH1D*) hJetPt_SE_trg9_C[i]->Clone();
        binmax=hJetPt_SE_trg9_C[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg9_C[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_scaled_trg9_C[i]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio_trg9_scaled_C[i]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE_trg9_C[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_scaled_trg9_C[i]->GetBinError(ibin);
            hRatio_trg9_scaled_C[i]->SetBinContent(ibin,y/yref);
            double err = (binerr/yref)-(y*binerref/(yref*yref));
            hRatio_trg9_scaled_C[i]->SetBinError(ibin,err);

        }

    //peripheral
        hRatio_trg9_scaled_P[i] = (TH1D*) hJetPt_SE_trg9_P[i]->Clone();
        binmax=hJetPt_SE_trg9_P[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg9_P[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_scaled_trg9_P[i]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio_trg9_scaled_P[i]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE_trg9_P[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_scaled_trg9_P[i]->GetBinError(ibin);
            hRatio_trg9_scaled_P[i]->SetBinContent(ibin,y/yref);
            double err = (binerr/yref)-(y*binerref/(yref*yref));
            hRatio_trg9_scaled_P[i]->SetBinError(ibin,err);

        }

    //trg7
    //centarl
        hRatio_trg7_scaled_C[i] = (TH1D*) hJetPt_SE_trg7_C[i]->Clone();
        binmax=hJetPt_SE_trg7_C[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg7_C[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_scaled_trg7_C[i]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio_trg7_scaled_C[i]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE_trg7_C[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_scaled_trg7_C[i]->GetBinError(ibin);
            hRatio_trg7_scaled_C[i]->SetBinContent(ibin,y/yref);
            double err = (binerr/yref)-(y*binerref/(yref*yref));
            hRatio_trg7_scaled_C[i]->SetBinError(ibin,err);

        }

    //peripheral
        hRatio_trg7_scaled_P[i] = (TH1D*) hJetPt_SE_trg7_P[i]->Clone();
        binmax=hJetPt_SE_trg7_P[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg7_P[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_scaled_trg7_P[i]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio_trg7_scaled_P[i]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE_trg7_P[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_scaled_trg7_P[i]->GetBinError(ibin);
            hRatio_trg7_scaled_P[i]->SetBinContent(ibin,y/yref);
            double err = (binerr/yref)-(y*binerref/(yref*yref));
            hRatio_trg7_scaled_P[i]->SetBinError(ibin,err);

        }

        
   
        

    }
    



}

void drawPtCompare(int _form){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "PtCompar_%d",_form);
    
    can[0]= new TCanvas(name,"Graph",10,10,1100,900);    
	pad[0][0]=new TPad("pad1","pad1",0.06,0.4,0.94,0.94);
    pad[0][0]->Draw();
	pad[0][1]=new TPad("pad1","pad1",0.06,0.06,0.94,0.4);
    pad[0][1]->Draw();
    
    pad[0][0]->SetTopMargin(0.08);
	pad[0][0]->SetBottomMargin(0);
    pad[0][0]->cd();
            
            
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);      //remove the entries,mean,RMS in the upper right.
    gStyle->SetOptTitle(0);
            
    hhtem[0] =  (TH1D*)    hJetPt_SE_trg7_C[_form]  ->Clone();
    hhtem[1] =  (TH1D*)    hJetPt_ME_trg7_C[_form]  ->Clone();
         
    
    hhtem[0]->SetMarkerStyle(29);
    hhtem[0]->SetMarkerSize(1);
    hhtem[0]->SetMarkerColor(2);
    hhtem[0]->SetLineColor(2);

    //hhtem[1]->SetMarkerStyle(20);
    //hhtem[1]->SetMarkerSize(0.7);
    //hhtem[1]->SetMarkerColor(4);
    hhtem[1]->SetLineColor(kGray+2);
    hhtem[1]->SetLineWidth(1.5);
    hhtem[1]->SetLineStyle(2);
    hhtem[1]->SetFillColor(kGray+2);
    hhtem[1]->SetFillStyle(3004);

    

            
            //char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="(1/N_{trig})d^{2}N_{jets}/(dp_{T,jet}^{reco,ch}d#eta) (Gev/c)^{-1} ";
            
            //hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            if(_form==0) hhtem[0] -> GetXaxis()->SetRangeUser(-5,30);
            if(_form==1) hhtem[0] -> GetXaxis()->SetRangeUser(-5,30);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            
            hhtem[0] -> GetYaxis()->SetLabelSize(0.05);
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.0);
            hhtem[0] -> GetYaxis()->SetTitleSize(0.05);
            hhtem[0]->GetYaxis()->CenterTitle(true);
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("hist same");
            //gr[1]->DrawClone("same F2");
            //hhtem[2]->DrawClone("same flx");
            gr[0]->DrawClone("same F");
 
            

            //hhtem[1]->DrawClone("same");
            
            tx0=0.50, ty0=0.84;
            if(_form==0) myTextF(tx0,ty0,"Isobar 200GeV, 0-10%",tsize*1.,1,12);
            if(_form==1) myTextF(tx0,ty0,"Isobar 200GeV, 60-80%",tsize*1.,1,12);
            
            tx0=0.50;ty0=0.76;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, R=0.2",tsize*1.,1,12);

	        tx0=0.50, ty0=0.69;
	        myTextF(tx0,ty0,"h^{#pm}+jet,7<P_{T}^{trg}<30GeV/c",tsize*1.,1,12);
            
            float _yy = 0.38;
            leg = mylegF(0.6,_yy,0.8,0.58,0.05);
            
            leg->AddEntry(hhtem[0],type[0],"pl");
	        leg->AddEntry(hhtem[1],"mix events(scaled)","l");
	    //leg->AddEntry(hhtem[2],"Pythia6 STAR tune","l");
	        leg->AddEntry(gr[0],"norm. region","f");
//            leg->AddEntry(gr[1],"Fill region range2","f");
            leg->Draw("same");
            
            
            
    pad[0][1]->SetTopMargin(0.00);
	pad[0][1]->SetBottomMargin(0.3);

    pad[0][1]->cd();

    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);    

    char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
    char *ytile1 ="SE/ME";

    hhtem[3] =  (TH1D*)    hRatio_scaled[_form] ->Clone();
    

    hhtem[3]->SetMarkerStyle(29);
    hhtem[3]->SetMarkerSize(0.7);
    hhtem[3]->SetMarkerColor(2);
    hhtem[3]->SetLineColor(2);

    

    hhtem[3] -> GetXaxis()->SetTitle(xtile);
    hhtem[3] -> GetYaxis()->SetTitle(ytile1);
    hhtem[3] -> GetXaxis()->SetNdivisions(507);
    hhtem[3] -> GetYaxis()->SetNdivisions(507);
    if(_form==0) hhtem[3] -> GetXaxis()->SetRangeUser(-5,30);
    if(_form==1) hhtem[3] -> GetXaxis()->SetRangeUser(-5,30);

    hhtem[3] -> GetYaxis()->SetTitleOffset(0.7); 
    hhtem[3] -> GetYaxis()->SetLabelSize(0.07);
    hhtem[3] -> GetXaxis()->SetLabelSize(0.07);
    hhtem[3]->GetYaxis()->SetTitleSize(0.07);             
    hhtem[3]->GetXaxis()->SetTitleSize(0.08);
    hhtem[3]->GetYaxis()->CenterTitle(true);           
    hhtem[3]->DrawClone("P");
 


    TLine *t1;
    if(_form==0) t1= new TLine(-5,1,30,1);
    if(_form==1) t1= new TLine(-5,1,30,1);
    t1->SetLineStyle(2);
	t1->Draw("same");  



    //______________________zoom in panel
    if(_form==0) pad[0][2]=new TPad("pad1","pad1",0.12,0.54,0.42,0.99);
    if(_form==1) pad[0][2]=new TPad("pad1","pad1",0.58,0.54,0.88,0.99);
    
    pad[0][2]->Draw();

    pad[0][2]->SetBottomMargin(0.35);
    pad[0][2]->SetLeftMargin(0.15);
    pad[0][2]->cd();

    gPad->SetFillStyle(4000); // make it transparent
    
    gPad->SetTickx(1);
    gPad->SetTicky(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);    

    char *xtile ="p_{T,jet}^{reco,ch} [Gev/c]";
    char *ytile1 ="SE/ME";

    hhtem[4] =  (TH1D*)    hRatio_scaled[_form] ->Clone();  

    hhtem[4]->SetMarkerStyle(29);
    hhtem[4]->SetMarkerSize(0.7);
    hhtem[4]->SetMarkerColor(2);
    hhtem[4]->SetLineColor(2);

    

    hhtem[4] -> GetXaxis()->SetTitle(xtile);
    hhtem[4] -> GetYaxis()->SetTitle(ytile1);
    hhtem[4] -> GetXaxis()->SetNdivisions(507);
    hhtem[4] -> GetYaxis()->SetNdivisions(507);
    if(_form==0)hhtem[4] -> GetXaxis()->SetRangeUser(regionC[0],regionC[1]-1);
    if(_form==1)hhtem[4] -> GetXaxis()->SetRangeUser(regionP[0],regionP[1]-1);

    hhtem[4] -> GetYaxis()->SetTitleOffset(0.6); 
    hhtem[4] -> GetXaxis()->SetTitleOffset(1.1);
    hhtem[4]->GetXaxis()->SetTitleSize(0.12);
    hhtem[4]->GetYaxis()->SetTitleSize(0.12);
    hhtem[4] -> GetYaxis()->SetLabelSize(0.11);
    hhtem[4] -> GetXaxis()->SetLabelSize(0.12);
    hhtem[4]->GetYaxis()->CenterTitle(true);  
    hhtem[4]->GetXaxis()->CenterTitle(true);          
    hhtem[4]->DrawClone("P");

    TLine *t2;
    if(_form==0) t2 = new TLine(regionC[0],1,regionC[1],1);
    if(_form==1) {
        t2 = new TLine(regionP[0],1,regionP[1],1);
        if(regionP[0]==0) t2 = new TLine(regionP[0],1,1,1);}
    t2->SetLineStyle(2);
	t2->Draw("same");
            
   
    
}



 void drawMEScale_detail(){
    
    readin();

    //0-0.2 1-0.3 2-0.4
    for(int j=0;j<3;j++){
        drawSESubME_R_Cen(j);
        drawSESubME_R_Per(j);
        //drawtest(j);
    }
   
}
   
