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

double Jet_R=0.2;
const int regiontype=3;
double pt_region_C[regiontype][2] = { {-4,0}, {-4,-1}, {0,2} };
double pt_region_P[regiontype][2] = { {0,3}, {0,2}, {0,0} };

TFile *fin;
TFile *fin_mix;

vector<int> pt_regionC;
vector<int> pt_regionP;

TFile *fout;

TLine* line;

//________0 central 1 peripheral_______________ 
TH1D *hJetPt_SE[2];
TH1D *hNtrigger_SE[2];


TH1D *hJetPt_ME[2];
TH1D *hNtrigger_ME[2];

TH1D *hJetPt_pythia[2];
TH1D *hNtrigger_pythia[2];

TH1D *hRatio[2];

TH1D *hRatio_pythia[2];


//_______region test________________________________________
TH1D *hJetPt_ME_scaled[2][regiontype];
TH1D *hRatio_scaled[2][regiontype];


TGraph *norm_region[2][regiontype];
TGraph *norm_region_scaled[2][regiontype];

TH1D *hJet_ME_Sub_SE[2][regiontype];
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
    
    //===============================  read in =====================================
    sprintf(name,"jet_SE_R2.root");
    cout<<name<<endl;
    fin = TFile::Open(name);

    float nraw[2];
    for(int icent=0;icent<2;icent++){
        sprintf(name,"recoil_jet_cent%d",icent);
        hJetPt_SE[icent]     = (TH1D*)fin->Get(name);

        sprintf(name,"number_of_trigger_cent%d",icent);
        hNtrigger_SE[icent]    = (TH1D*)fin->Get(name);

        nraw[icent] = hNtrigger_SE[icent] -> GetEntries();
        hJetPt_SE[icent] -> Scale(1./nraw[icent]);

    }
    //===============================  read in =====================================
    sprintf(name,"jet_ME_R2.root");
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    float nmix[2];
    for(int icent=0;icent<2;icent++){
        sprintf(name,"recoil_jet_cent%d",icent);
        hJetPt_ME[icent]     = (TH1D*)fin->Get(name);

        sprintf(name,"number_of_trigger_cent%d",icent);
        hNtrigger_ME[icent]    = (TH1D*)fin->Get(name);

        nmix[icent] = hNtrigger_ME[icent] -> GetEntries();
        hJetPt_ME[icent] -> Scale(1./nmix[icent]);

        //hJetPt_ME_scaled[icent] = (TH1D*) hJetPt_ME[icent]->Clone();      
    }

    
    cout<<"mix event "<<nmix[0]<<" same event "<<nraw[0]<<endl;
    cout<<"mix eventP "<<nmix[1]<<" same eventP "<<nraw[1]<<endl;
    //_____________________________________________________________________

    for(int cent=0;cent<2;cent++){
            Int_t binmax=hJetPt_SE[cent]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_SE[cent]->GetBinContent(ibin);
                double deltap = hJetPt_SE[cent]->GetBinWidth(ibin);
                double deltaeta = 2-(2*Jet_R);
                double binerr = hJetPt_SE[cent]->GetBinError(ibin);
                hJetPt_SE[cent]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_SE[cent]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }

            binmax=hJetPt_ME[cent]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_ME[cent]->GetBinContent(ibin);
                double deltap = hJetPt_ME[cent]->GetBinWidth(ibin);
                double deltaeta = 2-(2*Jet_R);
                double binerr = hJetPt_ME[cent]->GetBinError(ibin);
                hJetPt_ME[cent]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_ME[cent]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }


    }
    //_________________________without scale____________________________________
    
    for(int cent=0;cent<2;cent++){
        hRatio[cent] = (TH1D*) hJetPt_SE[cent]->Clone();
        Int_t binmax=hJetPt_SE[cent]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE[cent]->GetBinContent(ibin);
            double yref = hJetPt_ME[cent]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio[cent]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE[cent]->GetBinError(ibin);
            double binerref =hJetPt_ME[cent]->GetBinError(ibin);
            hRatio[cent]->SetBinContent(ibin,y/yref);
	    double err = sqrt( binerr*binerr/(yref*yref)+ binerref*binerref*y*y/pow(yref,4) );
            hRatio[cent]->SetBinError(ibin,err);

        }
    }

  


}



void drawratiofit(int _from){
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "ratio_%i",_from);
    
    can[0] = newDivCan2( name, Lmrg,llmrg, ratx,  raty, nx, ny, 400, 400 );
    
    for(int iy=0; iy<ny; iy++){
        for(int ix=0; ix<nx; ix++){
            
            int ipad = ix   + iy* nx;
            
            sprintf(name,"%s_pad_%i_%i",can[0]->GetName(),ix, iy);
            pad[ix][iy] = (TPad*) gROOT->FindObject(name);
            pad[ix][iy]->cd();
            
            gPad->SetTickx(1);
            gPad->SetTicky(1);
            //gPad->SetLogy(1);
            gStyle->SetOptStat(0);      //remove the entries,mean,RMS in the upper right.
            gStyle->SetOptTitle(0);
            
            hhtem[0] =  (TH1D*)    hRatio[_from]  ->Clone();
            
            int fbin = hhtem[0]->FindFirstBinAbove(0,1);
            int ran = (int) hhtem[0]->GetBinCenter(fbin);
            int ran_f = ran-1;
            if(_from==0) TF1 *f1= new TF1("f1","[0]",ran_f,-2);
            if(_from==1) TF1 *f1= new TF1("f1","[0]",ran_f,0);
            
            
            
            char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="SE/ME";
            
    hhtem[0]->SetMarkerStyle(25);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(1);
    hhtem[0]->SetLineColor(1);

            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            hhtem[0] -> GetXaxis()->SetRangeUser(-10,10);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,1);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
           
            
            
            hhtem[0]->DrawClone("P");
            hhtem[0]->Fit("f1","R");

            TLine *tl = new TLine(-10,1,10,1);
	        tl->Draw("same"); 
            
            tx0=0.50, ty0=0.80;
            myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.8,1,12);
            tx0=0.50;ty0=0.72;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, R=0.2",tsize*0.8,1,12);

	    tx0=0.50, ty0=0.65;
	    if(_from==0) myTextF(tx0,ty0,"h^{#pm}+jet, 0-10%",tsize*0.8,1,12);
        if(_from==1) myTextF(tx0,ty0,"h^{#pm}+jet, 60-80%",tsize*0.8,1,12);



/*
    
            leg = mylegF(0.50,0.45,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"0-10%","lp");

            leg->AddEntry(hhtem[1],"60-80%","lp");

            leg->AddEntry(hhtem[2],"Pythia6 STAR tune","lp");
            
            leg->DrawClone();
*/            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy

}

 void drawratio(){
    
    readin();
    for(int i=0;i<2;i++){
    drawratiofit(i);
    
    }
 
    
}
   
