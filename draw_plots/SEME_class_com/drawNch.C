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
TH1D *hNch_raw_TT730;
TH1D *hNch_raw_TT710;
TH1D *hNch_raw_TT1030;
TH1D *hNch_raw_TT47;
TH1D *hEta_raw;
TH1D *hPhi_raw;
TH1D *hPt_raw;
TH2D *hEtaPhi_raw;

TH1D *hNch_mix;
TH1D *hEta_mix;
TH1D *hPhi_mix;
TH1D *hPt_mix;
TH2D *hEtaPhi_mix;

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
int clr[] = {1,kRed+1, kBlue,kCyan+1,  kGreen+1,kGray+2, kYellow+1,kViolet+1};

const char* type[]={"same event", "mix event", "^{4}He+^{14}N 1 PeV", "^{4}He+^{14}N 10 PeV"};

void readin(){
    
    //===============================  read in =====================================
    sprintf(name,"cen.root");
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    sprintf(name,"Hntrk_class");
    hNch_raw     = (TH1D*)fin->Get(name);

    sprintf(name,"Hntrk_TT7_30");
    hNch_raw_TT730     = (TH1D*)fin->Get(name);

    sprintf(name,"Hntrk_TT7_10");
    hNch_raw_TT710     = (TH1D*)fin->Get(name);

    sprintf(name,"Hntrk_TT10_30");
    hNch_raw_TT1030     = (TH1D*)fin->Get(name);

    sprintf(name,"Hntrk_TT4_7");
    hNch_raw_TT47     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Heta");
    hEta_raw     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Hphi");
    hPhi_raw     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Hpt");
    hPt_raw     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Hphi_vs_eta");
    hEtaPhi_raw = (TH2D*)fin->Get(name);
    
    //===============================  read in =====================================
    sprintf(name,"Mixcen.root");
    cout<<name<<endl;
    fin_mix = TFile::Open(name);
    
    sprintf(name,"Hntrk_class");
    hNch_mix     = (TH1D*)fin_mix->Get(name);
    
    sprintf(name,"Heta");
    hEta_mix     = (TH1D*)fin_mix->Get(name);
    
    sprintf(name,"Hphi");
    hPhi_mix     = (TH1D*)fin_mix->Get(name);
    
    sprintf(name,"Hpt");
    hPt_mix     = (TH1D*)fin_mix->Get(name);
    
    sprintf(name,"Hphi_vs_eta");
    hEtaPhi_mix = (TH2D*)fin_mix->Get(name);
    
    float nraw = hNch_raw -> GetEntries();
    float nmix = hNch_mix -> GetEntries();
    float TT7_30 = hNch_raw_TT730 -> GetEntries();
    float TT7_15 = hNch_raw_TT710 -> GetEntries();
    float TT15_30 = hNch_raw_TT1030 -> GetEntries();
    float TT4_7 = hNch_raw_TT47 -> GetEntries();
    cout<<"mix event "<<nmix<<" same event "<<nraw<<endl;
    
    hEta_raw -> Scale(1./nraw);
    hEta_mix -> Scale(1./nmix);
    
    hPhi_raw -> Scale(1./nraw);
    hPhi_mix -> Scale(1./nmix);
    
    hPt_raw -> Scale(1./nraw);
    hPt_mix -> Scale(1./nmix);

    hNch_raw -> Scale(1./nraw);
    hNch_mix -> Scale(1./nmix);
    hNch_raw_TT47 -> Scale(1./TT4_7);
    hNch_raw_TT1030 -> Scale(1./TT15_30);
    hNch_raw_TT730 -> Scale(1./TT7_30);
    hNch_raw_TT710 -> Scale(1./TT7_15);

    hEtaPhi_raw -> Scale(1./nraw);
    hEtaPhi_mix-> Scale(1./nmix);
    
}

void drawNchCompare(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "NchCompare_%d");
    
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
            
            hhtem[0] =  (TH1D*)    hNch_raw  ->Clone();
            hhtem[1] =  (TH1D*)    hNch_mix  ->Clone();

            hhtem[5] =  (TH1D*)    hNch_raw_TT47  ->Clone();

            hhtem[2] =  (TH1D*)    hNch_raw_TT710  ->Clone();
            hhtem[3] =  (TH1D*)    hNch_raw_TT1030  ->Clone();
            hhtem[4] =  (TH1D*)    hNch_raw_TT730  ->Clone();
            //hhtem[1]->SetFillColor(46);
            //hhtem[1]->SetFillStyle(3005); //<=============================
            
            hhtem[0] -> SetLineColor(clr[0]);
            hhtem[0] -> SetMarkerColor(clr[0]);
            hhtem[0] -> SetMarkerStyle(20);
            hhtem[0] -> SetMarkerSize(0.6);
           hhtem[0] -> SetFillColor(clr[0]);

            hhtem[1] -> SetLineColor(clr[1]);
            hhtem[1] -> SetMarkerColor(clr[1]);
            hhtem[1] -> SetMarkerStyle(25);
            hhtem[1] -> SetMarkerSize(0.6);
            hhtem[1] -> SetFillColor(clr[1]);

            hhtem[2] -> SetLineColor(clr[2]);
            hhtem[2] -> SetMarkerColor(clr[2]);
            hhtem[2] -> SetMarkerStyle(20);
            hhtem[2] -> SetMarkerSize(0.6);
            hhtem[2] -> SetFillColor(clr[2]);


            hhtem[3] -> SetLineColor(clr[3]);
            hhtem[3] -> SetMarkerColor(clr[3]);
            hhtem[3] -> SetMarkerStyle(26);
            hhtem[3] -> SetMarkerSize(0.6);
            hhtem[3] -> SetFillColor(clr[3]);


            hhtem[4] -> SetLineColor(clr[4]);
            hhtem[4] -> SetMarkerColor(clr[4]);
            hhtem[4] -> SetMarkerStyle(24);
            hhtem[4] -> SetMarkerSize(0.6);
            hhtem[4] -> SetFillColor(clr[4]);


            hhtem[5] -> SetLineColor(clr[5]);
            hhtem[5] -> SetMarkerColor(clr[5]);
            hhtem[5] -> SetMarkerStyle(30);
            hhtem[5] -> SetMarkerSize(0.6);
            hhtem[5] -> SetFillColor(clr[5]);

            //hhtem[1] -> SetLineStyle(8);
            
            char *xtile ="N_{ch}";
            char *ytile ="N_{evt}/N";
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
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("same P");
            hhtem[5]->DrawClone("same P");
            hhtem[2]->DrawClone("same P");
            hhtem[3]->DrawClone("same P");
            hhtem[4]->DrawClone("same P");
            //hhtem[1]->DrawClone("same");
            
          /*  hhtem[0]->DrawClone("LF");
            hhtem[1]->DrawClone("same LF");
            hhtem[5]->DrawClone("same LF");
            hhtem[2]->DrawClone("same LF");
            hhtem[3]->DrawClone("same LF");
            hhtem[4]->DrawClone("same LF");
           */

            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
            
            tx0=0.3;ty0=0.72;
            //sprintf(name,"60-80%");
            myTextF(tx0,ty0,"centrality 0-10%",tsize*0.8,1,12);
            
            tx0=0.65, ty0=0.8;
            //myTextF(tx0,ty0,"2Sub",tsize*0.8,1,12);
            
            
                float _yy =0.6;
                leg = mylegF(0.6,_yy,0.8,0.8,0.03);
                leg->AddEntry(hhtem[0],"same events","lp");
                leg->AddEntry(hhtem[1],"mixed events","lp");
                leg->AddEntry(hhtem[5],"TT[4-7]","lp");
                leg->AddEntry(hhtem[2],"TT[7-10]","lp");
                leg->AddEntry(hhtem[3],"TT[10-30]","lp");
                leg->AddEntry(hhtem[4],"TT[7-30]","lp");
                leg->Draw("same");
            
            
            
            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy
    
}

void drawNch(){
    
    
        readin();
        
        //drawEtaPhiraw(ifrom);
        //drawEtaPhimix(ifrom);
        
        //drawEtaCompare(ifrom);
        //drawPhiCompare(ifrom);
        //drawpTCompare(ifrom);
        drawNchCompare();
        
    
    
    
}







