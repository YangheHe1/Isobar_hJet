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

TH1D *hNch_raw;
TH1D *hEta_raw;
TH1D *hPhi_raw;
TH1D *hPt_raw;
TH2D *hEtaPhi_raw;

TH1D *hNch_mix;
TH1D *hEta_mix;
TH1D *hPhi_mix;
TH1D *hPt_mix;
TH2D *hEtaPhi_mix;

TH1D *hentry_raw;
TH1D *hentry_mix;
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

const char* type[]={"same event", "mix event", "^{4}He+^{14}N 1 PeV", "^{4}He+^{14}N 10 PeV"};

void readin(){
    
    hentry_mix= new TH1D("hentry_mix","mix class entry",20,0,20);
    hentry_raw= new TH1D("hentry_raw","class entry",20,0,20);

    float all_raw;
    float all_mix;

    for(int _from=0;_from<16;_from++){
    //===============================  read in =====================================
    sprintf(name,"ntrk_per/MC_Out_from%d_to%d.root",_from, _from+1);
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    sprintf(name,"Hntrk");
    hNch_raw     = (TH1D*)fin->Get(name);
    
    //===============================  read in =====================================
    sprintf(name,"ntrk_Mixper/MC_Out_from%d_to%d.root",_from, _from+1);
    cout<<name<<endl;
    fin_mix = TFile::Open(name);
    
    sprintf(name,"Hntrk");
    hNch_mix     = (TH1D*)fin_mix->Get(name);
    
    
    float nraw = hNch_raw -> GetEntries();
    float nmix = hNch_mix -> GetEntries();
    all_mix += nmix;
    all_raw += nraw;
    int ibin=_from+1;
    hentry_mix->SetBinContent(ibin,nmix);
    hentry_raw->SetBinContent(ibin,nraw);
    }
	cout<<"nraw"<<nraw<<endl;
    hentry_mix->Scale(1./all_mix);
    hentry_raw->Scale(1./all_raw);

     
}
    


void drawNchCompare(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "NchCompare");
    
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
            
            hhtem[0] =  (TH1D*)    hentry_raw  ->Clone();
            hhtem[1] =  (TH1D*)    hentry_mix  ->Clone();
            
            //hhtem[1]->SetFillColor(46);
            //hhtem[1]->SetFillStyle(3005); //<=============================
            
            hhtem[0] -> SetLineColor(clr[0]);
            hhtem[1] -> SetLineColor(clr[1]);
            
            //hhtem[1] -> SetLineStyle(8);
            
            char *xtile ="class index";
            char *ytile ="(Events in class)/(Events in the cenrality)";
            char *ztile ="C(#Delta#eta #Delta#phi)";
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            //hhtem[0] -> GetXaxis()->SetRangeUser(0,1100);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
            /*
             hhtem[0]->GetXaxis()->CenterTitle(true);
             hhtem[0]->GetXaxis()->SetLabelFont(42);
             hhtem[0]->GetXaxis()->SetLabelSize(0.045);
             hhtem[0]->GetXaxis()->SetTitleSize(0.05);
             hhtem[0]->GetXaxis()->SetTitleOffset(1.5);
             hhtem[0]->GetXaxis()->SetTitleFont(42);
             hhtem[0]->GetXaxis()->SetNdivisions(507);
             
             hhtem[0]->GetYaxis()->CenterTitle(true);
             hhtem[0]->GetYaxis()->SetLabelFont(42);
             hhtem[0]->GetYaxis()->SetLabelSize(0.045);
             hhtem[0]->GetYaxis()->SetTitleSize(0.05);
             hhtem[0]->GetYaxis()->SetTitleOffset(1.5);
             hhtem[0]->GetYaxis()->SetTitleFont(42);
             */
            
            
            hhtem[0]->DrawClone("hist");
            hhtem[1]->DrawClone("same hist");
            //hhtem[1]->DrawClone("same");
            
            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
            
            tx0=0.3;ty0=0.72;
            //sprintf(name,"235 <N_{tracks}< 400");
            //myTextF(tx0,ty0,name,tsize*0.8,1,12);
            
            tx0=0.65, ty0=0.8;
            //myTextF(tx0,ty0,"2Sub",tsize*0.8,1,12);
            
            for(int k=0; k<2; k++){
                float _yy = k*0.08 + 0.65;
                leg = mylegF(0.6,_yy,0.8,0.8,0.03);
                leg->AddEntry(hhtem[k],type[k],"l");
                leg->Draw("same");
            }
            
            
            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy
    
}

void draw2D(){
    
    
        readin();
        
        //drawEtaPhiraw(ifrom);
        //drawEtaPhimix(ifrom);
        
        //drawEtaCompare(ifrom);
        //drawPhiCompare(ifrom);
        //drawpTCompare(ifrom);
        drawNchCompare();
        
   
    
}







