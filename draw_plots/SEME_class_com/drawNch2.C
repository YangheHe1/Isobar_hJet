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
    
    //===============================  read in =====================================
    sprintf(name,"ntrk_per.root");
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    sprintf(name,"Hntrk_class");
    hNch_raw     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Heta");
    hEta_raw     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Hphi");
    hPhi_raw     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Hpt");
    hPt_raw     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Hphi_vs_eta");
    hEtaPhi_raw = (TH2D*)fin->Get(name);
    
    //===============================  read in =====================================
    sprintf(name,"ntrk_Mixper.root");
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
    cout<<"mix event "<<nmix<<" same event "<<nraw<<endl;
    
    hEta_raw -> Scale(1./nraw);
    hEta_mix -> Scale(1./nmix);
    
    hPhi_raw -> Scale(1./nraw);
    hPhi_mix -> Scale(1./nmix);
    
    hPt_raw -> Scale(1./nraw);
    hPt_mix -> Scale(1./nmix);

    hNch_raw -> Scale(1./nraw);
    hNch_mix -> Scale(1./nmix);
    
    hEtaPhi_raw -> Scale(1./nraw);
    hEtaPhi_mix-> Scale(1./nmix);
    
}
    
void drawEtaPhiraw(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "00EtaPhi_raw%d");
    
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
            
            
            h2D[0] =  (TH2D*)    hEtaPhi_raw  ->Clone();
            
            char *ytile ="#eta";
            char *xtile ="#phi";
            //char *ztile ="C(#Delta#eta #Delta#phi)";
            
            h2D[0] -> GetXaxis()->SetTitle(xtile);
            h2D[0] -> GetYaxis()->SetTitle(ytile);
            
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
            
            h2D[0]->DrawClone("colz");
            
            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
            
            tx0=0.3;ty0=0.72;
            sprintf(name,"235 <N_{tracks}< 400");
            //myTextF(tx0,ty0,name,tsize*0.8,1,12);
            
            tx0=0.65, ty0=0.8;
            //myTextF(tx0,ty0,"2Sub",tsize*0.8,1,12);
            
            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy
    
}

void drawEtaPhimix(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "00EtaPhi_mix%d");
    
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
            
            
            h2D[0] =  (TH2D*)    hEtaPhi_raw  ->Clone();
            
            char *ytile ="#eta";
            char *xtile ="#phi";
            //char *ztile ="C(#Delta#eta #Delta#phi)";
            
            h2D[0] -> GetXaxis()->SetTitle(xtile);
            h2D[0] -> GetYaxis()->SetTitle(ytile);
            
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
            
            h2D[0]->DrawClone("colz");
            
            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
            
            tx0=0.3;ty0=0.72;
            sprintf(name,"235 <N_{tracks}< 400");
            //myTextF(tx0,ty0,name,tsize*0.8,1,12);
            
            tx0=0.65, ty0=0.8;
            //myTextF(tx0,ty0,"2Sub",tsize*0.8,1,12);
            
            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy
    
}

void drawEtaCompare(){
 
    int nx = 1;
    int ny = 1;
 
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
 
    sprintf(name, "EtaCompare_%d");
 
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
 
            hhtem[0] =  (TH1D*)    hEta_raw  ->Clone();
            hhtem[1] =  (TH1D*)    hEta_mix  ->Clone();
            
            hhtem[0] -> SetLineColor(clr[0]);
            hhtem[1] -> SetLineColor(clr[1]);
            
            hhtem[1]->SetFillColor(46);
            hhtem[1]->SetFillStyle(3005); //<=============================
            
            hhtem[1] -> SetLineStyle(8);
            
            char *xtile ="#eta";
            char *ytile ="counts(arb. units)";
            char *ztile ="C(#Delta#eta #Delta#phi)";
 
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            //hhtem[0] -> GetXaxis()->SetRangeUser(0,1100);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,2.5);
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
            sprintf(name,"235 <N_{tracks}< 400");
            //myTextF(tx0,ty0,name,tsize*0.8,1,12);
 
            tx0=0.65, ty0=0.8;
            //myTextF(tx0,ty0,"2Sub",tsize*0.8,1,12);
 
            for(int k=0; k<2; k++){
                float _yy = k*0.08 + 0.75;
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

void drawPhiCompare(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "PhiCompare_%d");
    
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
            
            hhtem[0] =  (TH1D*)    hPhi_raw  ->Clone();
            hhtem[1] =  (TH1D*)    hPhi_mix  ->Clone();
            
            hhtem[1]->SetFillColor(46);
            hhtem[1]->SetFillStyle(3005); //<=============================
            
            hhtem[0] -> SetLineColor(clr[0]);
            hhtem[1] -> SetLineColor(clr[1]);
            
            hhtem[1] -> SetLineStyle(8);
            
            char *xtile ="#phi";
            char *ytile ="counts(arb. units)";
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
            sprintf(name,"235 <N_{tracks}< 400");
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

void drawpTCompare(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "pTCompare_%d");
    
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
            
            hhtem[0] =  (TH1D*)    hPt_raw  ->Clone();
            hhtem[1] =  (TH1D*)    hPt_mix  ->Clone();
            
            
            hhtem[1]->SetFillColor(46);
            hhtem[1]->SetFillStyle(3005); //<=============================
            
            hhtem[1] -> SetLineStyle(8);
            
            hhtem[0] -> SetLineColor(clr[0]);
            hhtem[1] -> SetLineColor(clr[1]);
            
            char *xtile ="p_{T}";
            char *ytile ="counts(arb. units)";
            char *ztile ="C(#Delta#eta #Delta#phi)";
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            //hhtem[0] -> GetXaxis()->SetRangeUser(0,1100);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,1);
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
            sprintf(name,"235 <N_{tracks}< 400");
            //myTextF(tx0,ty0,name,tsize*0.8,1,12);
            
            tx0=0.65, ty0=0.8;
            //myTextF(tx0,ty0,"2Sub",tsize*0.8,1,12);
            
            for(int k=0; k<2; k++){
                float _yy = k*0.08 + 0.55;
                leg = mylegF(0.5,_yy,0.8,0.8,0.03);
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
            
            //hhtem[1]->SetFillColor(46);
            //hhtem[1]->SetFillStyle(3005); //<=============================
            
            hhtem[0] -> SetLineColor(clr[0]);
            hhtem[1] -> SetLineColor(clr[1]);
            
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
            
            
            hhtem[0]->DrawClone("hist");
            hhtem[1]->DrawClone("same hist");
            //hhtem[1]->DrawClone("same");
            
            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
            
            tx0=0.3;ty0=0.72;
            //sprintf(name,"60-80%");
            myTextF(tx0,ty0,"centrality 60-80%",tsize*0.8,1,12);
            
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

void drawNch2(){
    
    
        readin();
        
        //drawEtaPhiraw(ifrom);
        //drawEtaPhimix(ifrom);
        
        //drawEtaCompare(ifrom);
        //drawPhiCompare(ifrom);
        //drawpTCompare(ifrom);
        drawNchCompare();
        
    
    
    
}







