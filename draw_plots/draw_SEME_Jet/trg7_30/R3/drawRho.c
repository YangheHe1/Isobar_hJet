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

TH1D *hRho_SE;
TH1D *hRho_SEP;
TH1D *hNtrigger_SE;

TH1D *hRho_ME;
TH1D *hRho_MEP;
TH1D *hNtrigger_ME;

TH1D *hRatio;

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
int clr[] = {kBlue,kRed+1, kBlue, 1,  kGreen+1,kGray+2, kCyan+1, kYellow+1,kViolet+1};

const char* type[]={"same event", "mix event", "^{4}He+^{14}N 1 PeV", "^{4}He+^{14}N 10 PeV"};

void readin(){
    
    //===============================  read in =====================================
    sprintf(name,"jet_SE_R3.root");
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    sprintf(name,"Rho_cent0");
    hRho_SE     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Rho_cent1");
    hRho_SEP     = (TH1D*)fin->Get(name);
    
    //===============================  read in =====================================
    sprintf(name,"jet_ME_R3.root");
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    sprintf(name,"Rho_cent0");
    hRho_ME     = (TH1D*)fin->Get(name);

    sprintf(name,"Rho_cent1");
    hRho_MEP     = (TH1D*)fin->Get(name);
    
    float nraw = hRho_SE -> GetEntries();
    float nmix = hRho_ME -> GetEntries();
    float nrawP = hRho_SEP -> GetEntries();
    float nmixP = hRho_MEP -> GetEntries();
    cout<<"mix event "<<nmix<<" same event "<<nraw<<endl;
    cout<<"mix eventP "<<nmixP<<" same eventP "<<nrawP<<endl;
    
    hRho_SE -> Scale(1./nraw);
    hRho_ME -> Scale(1./nmix);

    hRho_SEP -> Scale(1./nrawP);
    hRho_MEP -> Scale(1./nmixP);

    hRatio = (TH1D*) hRho_SE->Clone();
    Int_t binmax=hRho_SE->FindLastBinAbove(0,1);
    for(int ibin=1;ibin<=binmax;ibin++){
        double y = hRho_SE->GetBinContent(ibin);
        double yref = hRho_ME->GetBinContent(ibin);
        if(y==0||yref==0)
			{hRatio->SetBinContent(ibin,0);
				continue;}
        double binerr = hRho_SE->GetBinError(ibin);
        double binerref =hRho_ME->GetBinError(ibin);
        hRatio->SetBinContent(ibin,y/yref);
        double err = (binerr/yref)-(y*binerref/(yref*yref));
        hRatio->SetBinError(ibin,err);

    }
    
    
}

void drawRhoCompare(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "RhoCompar");
    
    can[0]= new TCanvas(name,"Graph",10,10,1100,900);    
	pad[0][0]=new TPad("pad1","pad1",0.06,0.4,1,1);
    pad[0][0]->Draw();
	pad[0][1]=new TPad("pad1","pad1",0.06,0,1,0.4);
    pad[0][1]->Draw();
    
    pad[0][0]->SetTopMargin(0.08);
	pad[0][0]->SetBottomMargin(0);
    pad[0][0]->cd();
            
            
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);      //remove the entries,mean,RMS in the upper right.
    gStyle->SetOptTitle(0);
            
    hhtem[0] =  (TH1D*)    hRho_SE  ->Clone();
    hhtem[1] =  (TH1D*)    hRho_ME  ->Clone();
    hhtem[4] =  (TH1D*)    hRho_SEP  ->Clone();
    hhtem[5] =  (TH1D*)    hRho_MEP  ->Clone();
            
    hhtem[1]->SetFillColor(46);
    hhtem[1]->SetFillStyle(3005); //<=============================
            
    hhtem[0] -> SetLineColor(clr[0]);
    hhtem[1] -> SetLineColor(clr[1]);
            
    hhtem[1] -> SetLineStyle(8);

    hhtem[5]->SetFillColor(46);
    hhtem[5]->SetFillStyle(3005); //<=============================
            
    hhtem[4] -> SetLineColor(clr[0]);
    hhtem[5] -> SetLineColor(clr[1]);
            
    hhtem[5] -> SetLineStyle(8);
            
            char *xtile ="#rho";
            char *ytile ="N_{evt}/N";
            
            //hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetRangeUser(0,35);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.2);
            hhtem[1]->SetTitleSize(0.05);
            hhtem[0]->GetYaxis()->CenterTitle(true);
            
            
            hhtem[0]->DrawClone("hist");
            hhtem[1]->DrawClone("same hist");
            hhtem[4]->DrawClone("same hist");
            hhtem[5]->DrawClone("same hist");

            //hhtem[1]->DrawClone("same");
            
            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
            
            tx0=0.3;ty0=0.72;
            //myTextF(tx0,ty0,"centrality 60-80%",tsize*0.8,1,12);
            
            tx0=0.65, ty0=0.8;
            //myTextF(tx0,ty0,"centrality 0-10%",tsize*0.8,1,12);
            
            for(int k=0; k<2; k++){
                float _yy = k*0.08 + 0.65;
                leg = mylegF(0.6,_yy,0.8,0.8,0.05);
                leg->AddEntry(hhtem[k],type[k],"l");
                leg->Draw("same");
            }
            
            
    pad[0][1]->SetTopMargin(0.00);
	pad[0][1]->SetBottomMargin(0.3);

    pad[0][1]->cd();

    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);    

    char *xtile ="#rho";
    char *ytile ="SE/ME";

    hhtem[3] =  (TH1D*)    hRatio ->Clone();

    hhtem[3] -> SetLineColor(clr[0]);
    hhtem[3] -> GetXaxis()->SetTitle(xtile);
    hhtem[3] -> GetYaxis()->SetTitle(ytile);
    hhtem[3] -> GetXaxis()->SetNdivisions(507);
    hhtem[3] -> GetYaxis()->SetNdivisions(507);
    hhtem[3] -> GetXaxis()->SetRangeUser(0,35);
            
    hhtem[3] -> GetYaxis()->SetTitleOffset(1.2); 
    hhtem[3]->GetXaxis()->SetTitleSize(0.08);
    hhtem[3]->GetYaxis()->CenterTitle(true);           
    hhtem[3]->DrawClone("hist");

    TLine *tl = new TLine(0,1,35,1);
	tl->Draw("same");    
            
   
    
}

 void drawRho(){
    
    readin();
    
    drawRhoCompare();
    
}
   
