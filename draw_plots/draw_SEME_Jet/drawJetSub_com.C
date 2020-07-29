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

TH1D *hJetPt_SE_trg7_P[3];
TH1D *hNtrigger_SE_trg7_P[3];
TH1D *hJetPt_ME_trg7_P[3];
TH1D *hNtrigger_ME_trg7_P[3];
TH1D *hJetPt_ME_scaled_trg7_P[3];
TH1D *hJet_ME_Sub_SE_trg7_P[3];

TH1D *hJetPt_SE_trg9_C[3];
TH1D *hNtrigger_SE_trg9_C[3];
TH1D *hJetPt_ME_trg9_C[3];
TH1D *hNtrigger_ME_trg9_C[3];
TH1D *hJetPt_ME_scaled_trg9_C[3];
TH1D *hJet_ME_Sub_SE_trg9_C[3];

TH1D *hJetPt_SE_trg9_P[3];
TH1D *hNtrigger_SE_trg9_P[3];
TH1D *hJetPt_ME_trg9_P[3];
TH1D *hNtrigger_ME_trg9_P[3];
TH1D *hJetPt_ME_scaled_trg9_P[3];
TH1D *hJet_ME_Sub_SE_trg9_P[3];

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

        
    //____ME-SE jet___________________________________________
    

    //trg9
    //central

        hJet_ME_Sub_SE_trg9_C[i] = (TH1D*) hJetPt_SE_trg9_C[i]->Clone();
        Int_t binmax=hJetPt_SE_trg9_C[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg9_C[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_scaled_trg9_C[i]->GetBinContent(ibin);
            if(y-yref<0)
			    {hJet_ME_Sub_SE_trg9_C[i]->SetBinContent(ibin,0);
			    continue;}
            double binerr = hJetPt_SE_trg9_C[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_scaled_trg9_C[i]->GetBinError(ibin);
            hJet_ME_Sub_SE_trg9_C[i]->SetBinContent(ibin,y-yref);
            double err = binerr+binerref;
            hJet_ME_Sub_SE_trg9_C[i]->SetBinError(ibin,err);

        }
    
    //peripheral
        hJet_ME_Sub_SE_trg9_P[i] = (TH1D*) hJetPt_SE_trg9_P[i]->Clone();
        binmax=hJetPt_SE_trg9_P[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg9_P[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_scaled_trg9_P[i]->GetBinContent(ibin);
            if(y-yref<0)
			    {hJet_ME_Sub_SE_trg9_P[i]->SetBinContent(ibin,0);
			    continue;}
            double binerr = hJetPt_SE_trg9_P[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_scaled_trg9_P[i]->GetBinError(ibin);
            hJet_ME_Sub_SE_trg9_P[i]->SetBinContent(ibin,y-yref);
            double err = binerr+binerref;
            hJet_ME_Sub_SE_trg9_P[i]->SetBinError(ibin,err);

        }    

    //trg7
    //central

        hJet_ME_Sub_SE_trg7_C[i] = (TH1D*) hJetPt_SE_trg7_C[i]->Clone();
        binmax=hJetPt_SE_trg7_C[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg7_C[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_scaled_trg7_C[i]->GetBinContent(ibin);
            if(y-yref<0)
			    {hJet_ME_Sub_SE_trg7_C[i]->SetBinContent(ibin,0);
			    continue;}
            double binerr = hJetPt_SE_trg7_C[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_scaled_trg7_C[i]->GetBinError(ibin);
            hJet_ME_Sub_SE_trg7_C[i]->SetBinContent(ibin,y-yref);
            double err = binerr+binerref;
            hJet_ME_Sub_SE_trg7_C[i]->SetBinError(ibin,err);

        }
    
    //peripheral
        hJet_ME_Sub_SE_trg7_P[i] = (TH1D*) hJetPt_SE_trg7_P[i]->Clone();
        binmax=hJetPt_SE_trg7_P[i]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE_trg7_P[i]->GetBinContent(ibin);
            double yref = hJetPt_ME_scaled_trg7_P[i]->GetBinContent(ibin);
            if(y-yref<0)
			    {hJet_ME_Sub_SE_trg7_P[i]->SetBinContent(ibin,0);
			    continue;}
            double binerr = hJetPt_SE_trg7_P[i]->GetBinError(ibin);
            double binerref =hJetPt_ME_scaled_trg7_P[i]->GetBinError(ibin);
            hJet_ME_Sub_SE_trg7_P[i]->SetBinContent(ibin,y-yref);
            double err = binerr+binerref;
            hJet_ME_Sub_SE_trg7_P[i]->SetBinError(ibin,err);

        }         
        

    }
    //subtraction



}



void drawSESubME_trg7(int _form){
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "SESubME_trg7_%i",_form);
    
    can[0] = newDivCan2( name, Lmrg,llmrg, ratx,  raty, nx, ny, 400, 400 );
    
    for(int iy=0; iy<ny; iy++){
        for(int ix=0; ix<nx; ix++){
            
            int ipad = ix   + iy* nx;
            
            sprintf(name,"%s_pad_%i_%i",can[0]->GetName(),ix, iy);
            pad[ix][iy] = (TPad*) gROOT->FindObject(name);
            pad[ix][iy]->cd();
            
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            gPad->SetLogy(1);
            gStyle->SetOptStat(0);      //remove the entries,mean,RMS in the upper right.
            gStyle->SetOptTitle(0);
            
            if(_form==0){
            hhtem[0] =  (TH1D*)    hJet_ME_Sub_SE_trg7_C[0]  ->Clone();
            hhtem[1] =  (TH1D*)    hJet_ME_Sub_SE_trg7_C[1]  ->Clone();
            hhtem[2] =  (TH1D*)    hJet_ME_Sub_SE_trg7_C[2]  ->Clone();
            }
            else{
            hhtem[0] =  (TH1D*)    hJet_ME_Sub_SE_trg7_P[0]  ->Clone();
            hhtem[1] =  (TH1D*)    hJet_ME_Sub_SE_trg7_P[1]  ->Clone();
            hhtem[2] =  (TH1D*)    hJet_ME_Sub_SE_trg7_P[2]  ->Clone();

            }
            
            
            
            char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="(1/N_{trig})d^{2}N_{jets}/(dp_{T,jet}^{reco,ch}d#eta) (Gev/c)^{-1} ";
            
    hhtem[0]->SetMarkerStyle(22);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(kGreen+2);
    hhtem[0]->SetLineColor(kGreen+2);

    hhtem[1]->SetMarkerStyle(20);
    hhtem[1]->SetMarkerSize(0.7);
    hhtem[1]->SetMarkerColor(kBlue);
    hhtem[1]->SetLineColor(kBlue);

    hhtem[2]->SetMarkerStyle(25);
    hhtem[2]->SetMarkerSize(0.7);
    hhtem[2]->SetMarkerColor(1);
    hhtem[2]->SetLineColor(1);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetXaxis()->SetRangeUser(0,1100);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,1);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
           
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("Psame");
            hhtem[2]->DrawClone("Psame");
            
            tx0=0.50, ty0=0.84;
            myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.8,1,12);
            tx0=0.50;ty0=0.76;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, 7<P^{trg}_{T}<30 GeV/c",tsize*0.8,1,12);

	        tx0=0.50, ty0=0.69;
	        if(_form==0) myTextF(tx0,ty0,"h^{#pm}+jet,0-10%",tsize*0.8,1,12);
            if(_form==1) myTextF(tx0,ty0,"h^{#pm}+jet,60-80%",tsize*0.8,1,12);


            leg = mylegF(0.50,0.45,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"R=0.2","lp");

            leg->AddEntry(hhtem[1],"R=0.3","lp");

            leg->AddEntry(hhtem[2],"R=0.4","lp");
            
            leg->DrawClone();
            //sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy

}

void drawSESubME_trg9(int _form){
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "SESubME_trg9_%i",_form);
    
    can[0] = newDivCan2( name, Lmrg,llmrg, ratx,  raty, nx, ny, 400, 400 );
    
    for(int iy=0; iy<ny; iy++){
        for(int ix=0; ix<nx; ix++){
            
            int ipad = ix   + iy* nx;
            
            sprintf(name,"%s_pad_%i_%i",can[0]->GetName(),ix, iy);
            pad[ix][iy] = (TPad*) gROOT->FindObject(name);
            pad[ix][iy]->cd();
            
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            gPad->SetLogy(1);
            gStyle->SetOptStat(0);      //remove the entries,mean,RMS in the upper right.
            gStyle->SetOptTitle(0);
            
            if(_form==0){
            hhtem[0] =  (TH1D*)    hJet_ME_Sub_SE_trg9_C[0]  ->Clone();
            hhtem[1] =  (TH1D*)    hJet_ME_Sub_SE_trg9_C[1]  ->Clone();
            hhtem[2] =  (TH1D*)    hJet_ME_Sub_SE_trg9_C[2]  ->Clone();
            }
            else{
            hhtem[0] =  (TH1D*)    hJet_ME_Sub_SE_trg9_P[0]  ->Clone();
            hhtem[1] =  (TH1D*)    hJet_ME_Sub_SE_trg9_P[1]  ->Clone();
            hhtem[2] =  (TH1D*)    hJet_ME_Sub_SE_trg9_P[2]  ->Clone();

            }
            
            
            
            char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="(1/N_{trig})d^{2}N_{jets}/(dp_{T,jet}^{reco,ch}d#eta) (Gev/c)^{-1} ";
            
    hhtem[0]->SetMarkerStyle(20);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(kGreen+2);
    hhtem[0]->SetLineColor(kGreen+2);

    hhtem[1]->SetMarkerStyle(22);
    hhtem[1]->SetMarkerSize(0.7);
    hhtem[1]->SetMarkerColor(kBlue);
    hhtem[1]->SetLineColor(kBlue);

    hhtem[2]->SetMarkerStyle(25);
    hhtem[2]->SetMarkerSize(0.7);
    hhtem[2]->SetMarkerColor(1);
    hhtem[2]->SetLineColor(1);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetXaxis()->SetRangeUser(0,1100);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,1);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
           
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("Psame");
            hhtem[2]->DrawClone("Psame");
            
            tx0=0.50, ty0=0.84;
            myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.8,1,12);
            tx0=0.50;ty0=0.76;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, 9<P^{trg}_{T}<30 GeV/c",tsize*0.8,1,12);

	        tx0=0.50, ty0=0.69;
	        if(_form==0) myTextF(tx0,ty0,"h^{#pm}+jet,0-10%",tsize*0.8,1,12);
            if(_form==1) myTextF(tx0,ty0,"h^{#pm}+jet,60-80%",tsize*0.8,1,12);


            leg = mylegF(0.50,0.45,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"R=0.2","lp");

            leg->AddEntry(hhtem[1],"R=0.3","lp");

            leg->AddEntry(hhtem[2],"R=0.4","lp");
            
            leg->DrawClone();
            //sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy

}

void drawSESubME_R_Cen(int _form){
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "SESubME_R_Cen_%i",_form);
    
    can[0] = newDivCan2( name, Lmrg,llmrg, ratx,  raty, nx, ny, 400, 400 );
    
    for(int iy=0; iy<ny; iy++){
        for(int ix=0; ix<nx; ix++){
            
            int ipad = ix   + iy* nx;
            
            sprintf(name,"%s_pad_%i_%i",can[0]->GetName(),ix, iy);
            pad[ix][iy] = (TPad*) gROOT->FindObject(name);
            pad[ix][iy]->cd();
            
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            gPad->SetLogy(1);
            gStyle->SetOptStat(0);      //remove the entries,mean,RMS in the upper right.
            gStyle->SetOptTitle(0);
            
            
            hhtem[0] =  (TH1D*)    hJet_ME_Sub_SE_trg7_C[_form]  ->Clone();
            hhtem[1] =  (TH1D*)    hJet_ME_Sub_SE_trg9_C[_form]  ->Clone();
            
            
            
            
            
            char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="(1/N_{trig})d^{2}N_{jets}/(dp_{T,jet}^{reco,ch}d#eta) (Gev/c)^{-1} ";
            
    hhtem[0]->SetMarkerStyle(20);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(kRed);
    hhtem[0]->SetLineColor(kRed);

    hhtem[1]->SetMarkerStyle(25);
    hhtem[1]->SetMarkerSize(0.7);
    hhtem[1]->SetMarkerColor(kBlue+1);
    hhtem[1]->SetLineColor(kBlue+1);

    
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetXaxis()->SetRangeUser(0,1100);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,1);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
           
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("Psame");
            
            
            tx0=0.50, ty0=0.84;
            myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.8,1,12);
            tx0=0.50;ty0=0.76;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, 0-10%",tsize*0.8,1,12);

	        tx0=0.50, ty0=0.69;
	        if(_form==0) myTextF(tx0,ty0,"h^{#pm}+jet,R=0.2",tsize*0.8,1,12);
            if(_form==1) myTextF(tx0,ty0,"h^{#pm}+jet,R=0.3",tsize*0.8,1,12);
            if(_form==2) myTextF(tx0,ty0,"h^{#pm}+jet,R=0.4",tsize*0.8,1,12);


            leg = mylegF(0.50,0.45,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"7<P^{trg}_{T}<30 GeV/c","lp");

            leg->AddEntry(hhtem[1],"9<P^{trg}_{T}<30 GeV/c","lp");
            
            leg->DrawClone();
            //sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy

}

void drawSESubME_R_Per(int _form){
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "SESubME_R_Per_%i",_form);
    
    can[0] = newDivCan2( name, Lmrg,llmrg, ratx,  raty, nx, ny, 400, 400 );
    
    for(int iy=0; iy<ny; iy++){
        for(int ix=0; ix<nx; ix++){
            
            int ipad = ix   + iy* nx;
            
            sprintf(name,"%s_pad_%i_%i",can[0]->GetName(),ix, iy);
            pad[ix][iy] = (TPad*) gROOT->FindObject(name);
            pad[ix][iy]->cd();
            
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            gPad->SetLogy(1);
            gStyle->SetOptStat(0);      //remove the entries,mean,RMS in the upper right.
            gStyle->SetOptTitle(0);
            
            
            hhtem[0] =  (TH1D*)    hJet_ME_Sub_SE_trg7_P[_form]  ->Clone();
            hhtem[1] =  (TH1D*)    hJet_ME_Sub_SE_trg9_P[_form]  ->Clone();
            
            
            
            
            
            char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="(1/N_{trig})d^{2}N_{jets}/(dp_{T,jet}^{reco,ch}d#eta) (Gev/c)^{-1} ";
            
    hhtem[0]->SetMarkerStyle(20);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(kRed);
    hhtem[0]->SetLineColor(kRed);

    hhtem[1]->SetMarkerStyle(25);
    hhtem[1]->SetMarkerSize(0.7);
    hhtem[1]->SetMarkerColor(kBlue+1);
    hhtem[1]->SetLineColor(kBlue+1);

    
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetXaxis()->SetRangeUser(0,1100);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,1);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
           
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("Psame");
            
            
            tx0=0.50, ty0=0.84;
            myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.8,1,12);
            tx0=0.50;ty0=0.76;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, 60-80%",tsize*0.8,1,12);

	        tx0=0.50, ty0=0.69;
	        if(_form==0) myTextF(tx0,ty0,"h^{#pm}+jet,R=0.2",tsize*0.8,1,12);
            if(_form==1) myTextF(tx0,ty0,"h^{#pm}+jet,R=0.3",tsize*0.8,1,12);
            if(_form==2) myTextF(tx0,ty0,"h^{#pm}+jet,R=0.4",tsize*0.8,1,12);


            leg = mylegF(0.50,0.45,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"7<P^{trg}_{T}<30 GeV/c","lp");

            leg->AddEntry(hhtem[1],"9<P^{trg}_{T}<30 GeV/c","lp");
            
            leg->DrawClone();
            //sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy

}


 void drawJetSub_com(){
    
    readin();
    //0 central 1 peripheral
/*    for(int i=0;i<2;i++){
    
    drawSESubME_trg7(i);
    drawSESubME_trg9(i);
    
    }
    */
    //0-0.2 1-0.3 2-0.4
    for(int j=0;j<3;j++){
        drawSESubME_R_Cen(j);
        drawSESubME_R_Per(j);
    }
   
}
   
