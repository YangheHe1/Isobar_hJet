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
TFile *fin_ref;


TFile *fout;

TLine* line;

TH1D *hEff;
TH1D *hEtaEff;
TH1D *hPhiEff;
//H1D *hNtrigger[7];

TH1D *hEff_ref;
TH1D *hEtaEff_ref;
TH1D *hPhiEff_ref;
//TH1D *hNtrigger_py[7];


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
    
    sprintf(name,"2u1g_PtEff/pyEmb_2000k_charged_R0.2_central.root");
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    sprintf(name,"hpt_effi_pTl0");
    hEff     = (TH1D*)fin->Get(name);

    sprintf(name,"heta_effi_pTl0");
    hEtaEff     = (TH1D*)fin->Get(name);

    sprintf(name,"hphi_effi_pTl0");
    hPhiEff     = (TH1D*)fin->Get(name);
    

    //_____________________________
    sprintf(name,"2u1g_3dEff/pyEmb_2000k_charged_R0.2_central.root");
    cout<<name<<endl;
    fin_ref = TFile::Open(name);
    
    sprintf(name,"hpt_effi_pTl0");
    hEff_ref     = (TH1D*)fin_ref->Get(name);

    sprintf(name,"heta_effi_pTl0");
    hEtaEff_ref     = (TH1D*)fin_ref->Get(name);

    sprintf(name,"hphi_effi_pTl0");
    hPhiEff_ref     = (TH1D*)fin_ref->Get(name);
    
    
    
}

void drawEffCompare(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "Eff");
    
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
            
            hhtem[0] =  (TH1D*)    hEff  ->Clone();
            hhtem[1] =  (TH1D*)    hEff_ref  ->Clone();
            


            char *xtile ="P_{T,jet}^{ch}";
            char *ytile ="jet matching efficiency";
            
    hhtem[0]->SetMarkerStyle(25);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(1);
    hhtem[0]->SetLineColor(1);
    hhtem[0]->SetLineWidth(1.7);

    hhtem[1]->SetLineColor(clr[0]);
    hhtem[1]->SetLineWidth(1.7);



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
            
            

            
            tx0=0.50, ty0=0.80;
            //myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.6,1,12);
            tx0=0.50;ty0=0.72;
            //sprintf(name,"centrality 0-20%");




    
            leg = mylegF(0.50,0.50,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"Pt efficiency","l");

            leg->AddEntry(hhtem[1],"3D efficiency","l");

            leg->DrawClone();
     
            
        }//ix
    }//iy
    

}


void drawEtaEffCompare(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "etaEff");
    
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
            
            hhtem[0] =  (TH1D*)    hEtaEff  ->Clone();
            hhtem[1] =  (TH1D*)    hEtaEff_ref  ->Clone();
            


            char *xtile ="#eta";
            char *ytile ="jet matching efficiency";
            
    hhtem[0]->SetMarkerStyle(25);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(1);
    hhtem[0]->SetLineColor(1);
    hhtem[0]->SetLineWidth(1.7);

    hhtem[1]->SetLineColor(clr[0]);
    hhtem[1]->SetLineWidth(1.7);



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
            
            

            
            tx0=0.50, ty0=0.80;
            //myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.6,1,12);
            tx0=0.50;ty0=0.72;
            //sprintf(name,"centrality 0-20%");




    
            leg = mylegF(0.50,0.50,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"Pt efficiency","l");

            leg->AddEntry(hhtem[1],"3D efficiency","l");

            leg->DrawClone();
     
            
        }//ix
    }//iy
    
}


void drawPhiEffCompare(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "phiEff");
    
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
            
            hhtem[0] =  (TH1D*)    hPhiEff  ->Clone();
            hhtem[1] =  (TH1D*)    hPhiEff_ref  ->Clone();
            


            char *xtile ="#phi";
            char *ytile ="jet matching efficiency";
            
    hhtem[0]->SetMarkerStyle(25);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(1);
    hhtem[0]->SetLineColor(1);
    hhtem[0]->SetLineWidth(1.7);

    hhtem[1]->SetLineColor(clr[0]);
    hhtem[1]->SetLineWidth(1.7);



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
            
            

            
            tx0=0.50, ty0=0.80;
            //myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.6,1,12);
            tx0=0.50;ty0=0.72;
            //sprintf(name,"centrality 0-20%");




    
            leg = mylegF(0.50,0.50,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"Pt efficiency","l");

            leg->AddEntry(hhtem[1],"3D efficiency","l");

            leg->DrawClone();
     
            
        }//ix
    }//iy
    
}


 void drawEffcom(){
    
    readin();
    
    drawEffCompare();
    drawEtaEffCompare();
    drawPhiEffCompare();
    
}
   
