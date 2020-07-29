#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TMath.h"
#include "TF1.h"
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h> // strtod(), string to double
#include "drawheader_v4.h"
#include <iostream>
#include <fstream>

using namespace std;


const int NBINS = 15;
const int NRT   = 4;

const float Xlow  = 1;
const float Xhigh = 500;

const double PI = acos(-1.0);

TFile *fin;

TFile *fout;

TLine* line;

TH2D *HTOF_Nch;

TH1D *HNCH[800];

TH1D *Hmean;
TH1D *Hhigh;
TH1D *Hlow;

TF1 *fit_mean;
TF1 *fit_low;
TF1 *fit_high;

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
int clr[] = {1,kRed+1, kBlue, 1,  kGreen+1,kGray+2, kCyan+1, kYellow+1,kViolet+1};

Double_t fit_function(Double_t* x,Double_t* par)
{
    // Track reconstruction efficiency parametrization
    Double_t n,y;
    Double_t A,B,C,D,E;A=par[0];B=par[1];C=par[2];D=par[3];E=par[4];
    n=x[0];
    y=A+B*n+C*n*n+D*pow(n,3)+E*pow(n,4);
    return y;
}

Double_t gausf(Double_t* x,Double_t* par){
    return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2]);
}

void readin(){
    sprintf(name,"Isobar_QA_fit.root");
    cout<<name<<endl;
    fin = TFile::Open(name);

    sprintf(name,"Htof_ntrk");
    HTOF_Nch     = (TH2D*)fin->Get(name);

    int npTbins=450;
	double pTmax=450;

    /*
    for(int i=0;i<npTbins;i++){
        sprintf(name,"htof_%i",i);
        HNCH[i] = new TH1D(name, name, npTbins, 0, pTmax);
    }
*/
    sprintf(name,"hmean");
    Hmean = new TH1D(name, name, npTbins, 0, pTmax);

    sprintf(name,"hlow");
    Hlow = new TH1D(name, name, npTbins, 0, pTmax);

    sprintf(name,"hhigh");
    Hhigh = new TH1D(name, name, npTbins, 0, pTmax);

    TF1 *func = new TF1("func",gausf,0,npTbins,3);

    int ntof = HTOF_Nch->FindLastBinAbove(0,1);
    int nnch = HTOF_Nch->FindLastBinAbove(0,2);

    int ftof = HTOF_Nch->FindFirstBinAbove(0,1);
    int fnch = HTOF_Nch->FindFirstBinAbove(0,2);

    double yl,yh;

    int binnum=23;
    double fitbin[23]={0,5,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,240,280,320,360,400,450};
        for(int k=0; k<binnum-1; k++){
            
            int i=fitbin[k]+1;
            int j=fitbin[k+1]+1;
            sprintf(name,"htof_%i",k);
            HNCH[k]= HTOF_Nch->ProjectionY(name,i,j);
            
            
            func->SetParameters(1,HNCH[k]->GetMean(),HNCH[k]->GetRMS());
            func->SetParNames("Constant","Mean_value","Sigma");
            if(HNCH[k]->Integral()==0) continue;
            HNCH[k]->Fit("func","0");
            TF1 *para = HNCH[k]->GetFunction("func");
            
            double mean = para->GetParameter(1);
            double smear= para->GetParameter(2);

            double xbins = (fitbin[k+1]-fitbin[k])/2+fitbin[k];
            int nbin=1+(Int_t)xbins;


            if(mean<=0) continue;
            Hmean ->SetBinContent(nbin,mean);
            
            if(k>0){ if((mean+4*smear)<(yh)) continue;}
            yh= mean+4*smear;
            Hhigh ->SetBinContent(nbin,yh);
            
            if(k>0){ if((mean-4*smear)<(yl)) continue;}
            yl= mean-4*smear;
            if(yl<=0) continue;
            Hlow ->SetBinContent(nbin,yl);
        }

    TF1 *fun_mean = new TF1("fun_mean",fit_function,0,npTbins,5);
    TF1 *fun_low = new TF1("fun_low",fit_function,0,npTbins,5);
    TF1 *fun_high = new TF1("fun_high",fit_function,0,npTbins,5);

    Hmean ->Fit("fun_mean","R");
    fit_mean = Hmean ->GetFunction("fun_mean");

    Hlow ->Fit("fun_low","R");
    fit_low = Hlow ->GetFunction("fun_low");

    Hhigh ->Fit("fun_high","R0");
    fit_high = Hhigh ->GetFunction("fun_high");

    ofstream outf("parametre.txt");

    for(int i=0;i<5;i++){
        double a=fit_mean->GetParameter(i);
        double l=fit_low->GetParameter(i);
        double h=fit_high->GetParameter(i);
        outf<<"parametre "<<i<<"    a: "<<a<<"  low: "<<l<<"  high: "<<h<<endl;

    }

}

void draw2d(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "00EtaPhi_raw");
    
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
            
            
            h2D[0] =  (TH2D*)    HTOF_Nch  ->Clone();
            
            char *xtile ="TOFMatch";
            char *ytile ="charged particles";
            //char *ztile ="C(#Delta#eta #Delta#phi)";
            
            h2D[0] -> GetXaxis()->SetTitle(xtile);
            h2D[0] -> GetYaxis()->SetTitle(ytile);

            h2D[0] -> GetXaxis()->SetRangeUser(0,450);
            h2D[0] -> GetYaxis()->SetRangeUser(0,700);
            /*
            h2D[0]->GetXaxis()->CenterTitle(true);
            h2D[0]->GetXaxis()->SetLabelFont(42);
            h2D[0]->GetXaxis()->SetLabelSize(0.045);
            h2D[0]->GetXaxis()->SetTitleSize(0.05);
            h2D[0]->GetXaxis()->SetTitleOffset(1);
            h2D[0]->GetXaxis()->SetTitleFont(42);
 
            h2D[0]->GetYaxis()->CenterTitle(true);
            h2D[0]->GetYaxis()->SetLabelFont(42);
            h2D[0]->GetYaxis()->SetLabelSize(0.045);
            h2D[0]->GetYaxis()->SetTitleSize(0.05);
            h2D[0]->GetYaxis()->SetTitleOffset(1);
            h2D[0]->GetYaxis()->SetTitleFont(42);

            h2D[0]->GetZaxis()->CenterTitle(true);
            h2D[0]->GetZaxis()->SetLabelFont(42);
            h2D[0]->GetZaxis()->SetLabelSize(0.045);
            h2D[0]->GetZaxis()->SetTitleSize(0.05);
            h2D[0]->GetZaxis()->SetTitleOffset(1.5);
            h2D[0]->GetZaxis()->SetTitleFont(42);
            */
    


            fit_mean->SetLineColor(clr[1]);
            fit_low->SetLineColor(clr[1]);
            fit_high->SetLineColor(clr[1]);
	/*
            fit_mean->SetRange(7,80);
            fit_low->SetRange(7,80);
            fit_high->SetRange(7,80);
	*/
            h2D[0]->DrawClone("colz");
            fit_mean->Draw("same");
            fit_low->Draw("same");
            fit_high->Draw("same");
            
            
        }//ix
    }//iy
    
}

void drawnFit(){
    
    
    readin();
    draw2d();

}
