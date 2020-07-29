#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h> // strtod(), string to double
#include "drawheader_v4.h"

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

TH1D *hRho[7];
TH1D *hNtrigger[7];

TH1D *hRho_ME[7];
TH1D *hNtrigger_ME[7];


TLegend *leg;

TGraphErrors *grtmpErr0[100];

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
int clr[] = {kBlue,kRed+1,kGreen+2,kYellow+2,kViolet+1};

const char* type[]={"same event", "mix event", "^{4}He+^{14}N 1 PeV", "^{4}He+^{14}N 10 PeV"};

void readin(){
    
    //===============================  read in =====================================
    for(int i=2;i<7;i++){
    sprintf(name,"SE_rho/jet_SE_R%i.root",i,i);
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    sprintf(name,"Hrho_all");
    hRho[i]     = (TH1D*)fin->Get(name);
    
    float nraw = hRho[i] -> GetEntries();
    
    
    hRho[i] -> Scale(1./nraw);
    
    }
    
    
    
}

void drawRhoCompare(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "rho");
    
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
            
            hhtem[0] =  (TH1D*)    hRho[2]  ->Clone();
            hhtem[1] =  (TH1D*)    hRho[3]  ->Clone();
            hhtem[2] =  (TH1D*)    hRho[4]  ->Clone();
            hhtem[3] =  (TH1D*)    hRho[5]  ->Clone();
            hhtem[4] =  (TH1D*)    hRho[6]  ->Clone();


            char *xtile ="#rho";
            char *ytile ="N_{evt}/N";
            
    hhtem[0]->SetMarkerStyle(25);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(1);
    hhtem[0]->SetLineColor(1);
    hhtem[0]->SetLineWidth(1.7);

    hhtem[1]->SetLineColor(clr[0]);
    hhtem[1]->SetLineWidth(1.7);

    hhtem[2]->SetLineColor(clr[1]);
    hhtem[2]->SetLineWidth(1.7);

    hhtem[3]->SetLineColor(clr[2]);
    hhtem[3]->SetLineWidth(1.7);

    hhtem[4]->SetLineColor(clr[3]);
    hhtem[4]->SetLineWidth(1.7);



            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            hhtem[0] -> GetXaxis()->SetRangeUser(0,50);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,1);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
           
            
            
            hhtem[0]->DrawClone("hist");
            hhtem[1]->DrawClone("hist same");
            hhtem[2]->DrawClone("hist same");
            hhtem[3]->DrawClone("hist same");
            hhtem[4]->DrawClone("hist same");

            

            
            tx0=0.50, ty0=0.80;
            myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.6,1,12);
            tx0=0.50;ty0=0.72;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"k_{T} algorithm",tsize*0.6,1,12);

	    tx0=0.50, ty0=0.65;
	    myTextF(tx0,ty0,"h^{#pm}+jet,same events",tsize*0.6,1,12);
        




    
            leg = mylegF(0.50,0.45,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"R=0.2","l");

            leg->AddEntry(hhtem[1],"R=0.3","l");

            leg->AddEntry(hhtem[2],"R=0.4","l");

            leg->AddEntry(hhtem[3],"R=0.5","l");

            leg->AddEntry(hhtem[4],"R=0.6","l");
            
            leg->DrawClone();
            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy
    
}

 void drawRhocom(){
    
    readin();
    
    drawRhoCompare();
    
}
   
