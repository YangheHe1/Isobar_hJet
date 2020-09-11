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
TFile *fin_ref;

TFile *fout;

TLine* line;

TH2D *HJetRes;
TH2D *HJetRes_ref;

TH1D *HRes[800];
TH1D *HRes_ref[800];

TH1D *Hmean;
TH1D *Hmean_ref;

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


Double_t gausf(Double_t* x,Double_t* par){
    return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2]);
}

void readin(){

    //_read in______________
    sprintf(name,"2u1g_PtEff/pyEmb_2000k_charged_R0.2_central.root");
    cout<<name<<endl;
    fin = TFile::Open(name);

    sprintf(name,"hJER_pTl0");
    HJetRes     = (TH2D*)fin->Get(name);

    //read in________________
    sprintf(name,"2u1g_3dEff/pyEmb_2000k_charged_R0.2_central.root");
    cout<<name<<endl;
    fin_ref = TFile::Open(name);

    sprintf(name,"hJER_pTl0");
    HJetRes_ref     = (TH2D*)fin_ref->Get(name);

    int npTbins=100;
	double pTmax=100;



    /*
    for(int i=0;i<npTbins;i++){
        sprintf(name,"htof_%i",i);
        HNCH[i] = new TH1D(name, name, npTbins, 0, pTmax);
    }
*/
    sprintf(name,"hmean");
    Hmean = new TH1D(name, name, npTbins, 0, pTmax);

    sprintf(name,"hmean_ref");
    Hmean_ref = new TH1D(name, name, npTbins, 0, pTmax);

    TF1 *func = new TF1("func",gausf,0,npTbins,3);
    TF1 *func_ref = new TF1("func_ref",gausf,0,npTbins,3);

    int ntof = HJetRes->FindLastBinAbove(0,1);
    int nnch = HJetRes->FindLastBinAbove(0,2);

    int ftof = HJetRes->FindFirstBinAbove(0,1);
    int fnch = HJetRes->FindFirstBinAbove(0,2);

    double yl,yh;

    int binnum=28;
    double fitbin[28]={0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,24,28,32,36,40,45,50,60,70,80,90,99};
        for(int k=0; k<binnum-1; k++){
            
            int i=fitbin[k]+1;
            int j=fitbin[k+1]+1;
            sprintf(name,"hres_%i",k);
            HRes[k]= HJetRes->ProjectionX(name,i,j);

            sprintf(name,"hres_ref_%i",k);
            HRes_ref[k]= HJetRes_ref->ProjectionX(name,i,j);
            
            
            func->SetParameters(1,HRes[k]->GetMean(),HRes[k]->GetRMS());
            func->SetParNames("Constant","Mean_value","Sigma");
            if(HRes[k]->Integral()==0) continue;
            HRes[k]->Fit("func","0");
            TF1 *para = HRes[k]->GetFunction("func");
            
            double mean = para->GetParameter(1);

            double xbins = (fitbin[k+1]-fitbin[k])/2+fitbin[k];
            int nbin=1+(Int_t)xbins;

            Hmean ->SetBinContent(nbin,mean);


            func_ref->SetParameters(1,HRes_ref[k]->GetMean(),HRes_ref[k]->GetRMS());
            func_ref->SetParNames("Constant","Mean_value","Sigma");
            if(HRes_ref[k]->Integral()==0) continue;
            HRes_ref[k]->Fit("func_ref","0");
            TF1 *para_ref = HRes_ref[k]->GetFunction("func_ref");
            
            double mean_ref = para_ref->GetParameter(1);

            Hmean_ref ->SetBinContent(nbin,mean_ref);            
            
            
        }

}

void drawMeanCom(){
    
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
            
            
            hhtem[0] =  (TH1D*)    Hmean  ->Clone();
            hhtem[1] =  (TH1D*)    Hmean_ref  ->Clone();

            char *xtile ="P_{T,jet}^{part}";
            char *ytile ="<(P_{T,jet}^{dete}-P_{T,jet}^{part})/P_{T,jet}^{part}>";
            //char *ztile ="C(#Delta#eta #Delta#phi)";
            
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

            tx0=0.50, ty0=0.70;
            myTextF(tx0,ty0,"centrality 0-10%,R=0.2",tsize*0.6,1,12);
            

    
            leg = mylegF(0.50,0.50,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"Pt efficiency","l");

            leg->AddEntry(hhtem[1],"3D efficiency","l");

            leg->DrawClone();

            
            
        }//ix
    }//iy
    
}

void drawnResCom(){
    
    
    readin();
    drawMeanCom();

}
