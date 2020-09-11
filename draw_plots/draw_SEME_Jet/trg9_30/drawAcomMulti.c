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
TFile *fin_ME;
TFile *fin_py;


TFile *fout;

TLine* line;

TH1D *hArea[7];
TH1D *hArea_Pt5[7];
TH1D *hNtrigger[7];

TH1D *hArea_ME[7];
TH1D *hArea_ME_Pt5[7];
TH1D *hNtrigger_ME[7];

TH1D *hArea_py[7];
TH1D *hArea_py_Pt5[7];
TH1D *hNtrigger_py[7];


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
    for(int i=2;i<6;i++){
    sprintf(name,"SE_A/jet_SE_R%i.root",i);
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    sprintf(name,"Harea7_30");
    hArea[i]     = (TH1D*)fin->Get(name);

    sprintf(name,"Harea7_30_Pt5");
    hArea_Pt5[i]     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Htrg7_30");
    hNtrigger[i]     = (TH1D*)fin->Get(name);

    float nraw = hNtrigger[i] -> GetEntries();
    
    
    hArea[i] -> Scale(1./nraw);
    hArea_Pt5[i] -> Scale(1./nraw);

    //=======================
    sprintf(name,"ME_A/jet_ME_R%i.root",i);
    cout<<name<<endl;
    fin_ME = TFile::Open(name);
    
    sprintf(name,"Harea");
    hArea_ME[i]     = (TH1D*)fin_ME->Get(name);

    sprintf(name,"Harea_Pt5");
    hArea_ME_Pt5[i]     = (TH1D*)fin_ME->Get(name);
    
    sprintf(name,"HNevt");
    hNtrigger_ME[i]     = (TH1D*)fin_ME->Get(name);

    float nmix = hNtrigger_ME[i] -> GetEntries();
    
    
    hArea_ME[i] -> Scale(1./nmix);
    hArea_ME_Pt5[i] -> Scale(1./nmix);

    //_____________________________
    sprintf(name,"Pythia6_A/jet_pythia6_R%i.root",i);
    cout<<name<<endl;
    fin_py = TFile::Open(name);
    
    sprintf(name,"Harea7_30");
    hArea_py[i]     = (TH1D*)fin_py->Get(name);

    sprintf(name,"Harea7_30_Pt5");
    hArea_py_Pt5[i]     = (TH1D*)fin_py->Get(name);
    
    sprintf(name,"Htrg7_30");
    hNtrigger_py[i]     = (TH1D*)fin_py->Get(name);

    float npy = hNtrigger_py[i] -> GetEntries();
    
    
    hArea_py[i] -> Scale(1./npy);
    hArea_py_Pt5[i] -> Scale(1./npy);


    //scale pythia area
    int maxBin;
    double SE_max;
    //maxBin = hArea[i]->GetMaximumBin();
    //SE_max = hArea[i]->GetBinContent(maxBin);

    maxBin = hArea_Pt5[i]->GetMaximumBin();
    SE_max = hArea[i]->GetBinContent(maxBin);
    double SE_pt5_max = hArea_Pt5[i]->GetBinContent(maxBin);
    double SE_pt5_scaler = SE_max/SE_pt5_max;
    hArea_Pt5[i]->Scale(SE_pt5_scaler);

    /*
    double ME_max = hArea_ME[i]->GetBinContent(maxBin);
    double ME_scaler = SE_max/ME_max;
    hArea_ME[i]->Scale(ME_scaler);
    */

    maxBin = hArea_ME_Pt5[i]->GetMaximumBin();
    SE_max = hArea[i]->GetBinContent(maxBin);
    double ME_pt5_max = hArea_ME_Pt5[i]->GetBinContent(maxBin);
    double ME_pt5_scaler = SE_max/ME_pt5_max;
    hArea_ME_Pt5[i]->Scale(ME_pt5_scaler);

    maxBin = hArea_py[i]->GetMaximumBin();
    SE_max = hArea[i]->GetBinContent(maxBin);
    double py_max = hArea_py[i]->GetBinContent(maxBin);
    double py_scaler = SE_max/py_max;
    hArea_py[i]->Scale(py_scaler);

    maxBin = hArea_py_Pt5[i]->GetMaximumBin();
    SE_max = hArea[i]->GetBinContent(maxBin);
    double py_pt5_max = hArea_py_Pt5[i]->GetBinContent(maxBin);
    double py_pt5_scaler = SE_max/py_pt5_max;
    hArea_py_Pt5[i]->Scale(py_pt5_scaler);

    
    
    }
    
    

    
}

void drawACompare(int _form){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "A_R%i",_form);
    
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
            
            hhtem[0] =  (TH1D*)    hArea[_form]  ->Clone();
            hhtem[1] =  (TH1D*)    hArea_ME[_form]  ->Clone();
            hhtem[2] =  (TH1D*)    hArea_py[_form]  ->Clone();

            hhtem[3] =  (TH1D*)    hArea_Pt5[_form]  ->Clone();
            hhtem[4] =  (TH1D*)    hArea_ME_Pt5[_form]  ->Clone();
            hhtem[5] =  (TH1D*)    hArea_py_Pt5[_form]  ->Clone();

            


            char *xtile ="A";
            char *ytile ="N_{jet}/N_{trg}";
            
    hhtem[0]->SetMarkerStyle(25);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(kRed);
    hhtem[0]->SetLineColor(kRed);
    hhtem[0]->SetLineWidth(1.7);

    hhtem[1]->SetLineColor(kGreen+2);
    hhtem[1]->SetLineWidth(1.7);
    hhtem[1]->SetFillColor(kGreen+2);
    hhtem[1]->SetFillStyle(3005);

    hhtem[2]->SetLineColor(kViolet+1);
    hhtem[2]->SetLineWidth(1.7);

    hhtem[3]->SetLineColor(kRed);
    hhtem[3]->SetLineWidth(1.7);
    hhtem[3]->SetLineStyle(2);

    hhtem[4]->SetLineColor(kGreen+2);
    hhtem[4]->SetLineWidth(1.7);
    hhtem[4]->SetLineStyle(2);

    hhtem[5]->SetLineColor(kViolet+1);
    hhtem[5]->SetLineWidth(1.7);
    hhtem[5]->SetLineStyle(2);    


            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            hhtem[0] -> GetXaxis()->SetRangeUser(0,50);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,1);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
           
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("hist same");
            hhtem[2]->DrawClone("hist same");
            hhtem[3]->DrawClone("hist same");
            hhtem[4]->DrawClone("hist same");
            hhtem[5]->DrawClone("hist same");
            

            
            tx0=0.50, ty0=0.80;
            myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.6,1,12);
            tx0=0.50;ty0=0.72;
            //sprintf(name,"centrality 0-20%");

            if(_form==2) myTextF(tx0,ty0,"anti_k_{T} algorithm,R=0.2",tsize*0.6,1,12);
            if(_form==3) myTextF(tx0,ty0,"anti_k_{T} algorithm,R=0.3",tsize*0.6,1,12);
            if(_form==4) myTextF(tx0,ty0,"anti_k_{T} algorithm,R=0.4",tsize*0.6,1,12);
            if(_form==5) myTextF(tx0,ty0,"anti_k_{T} algorithm,R=0.5",tsize*0.6,1,12);

	    tx0=0.50, ty0=0.65;
	    myTextF(tx0,ty0,"h^{#pm}+jet,7<P_{trg,T}<30 GeV/c",tsize*0.6,1,12);
        




    
            leg = mylegF(0.28,0.48,0.55,0.78,0.03);
            leg->AddEntry(hhtem[0],"same events","lp");
            leg->AddEntry(hhtem[1],"mixed events","lf");
            leg->AddEntry(hhtem[2],"pythia6 star tune","l");
            leg->AddEntry(hhtem[0],"same events,P^{ch}_{T,jet}>5GeV/c","l");
            leg->AddEntry(hhtem[1],"mixed events,P^{ch}_{T,jet}>5GeV/c","l");
            leg->AddEntry(hhtem[2],"pythia6 star tune,P^{ch}_{T,jet}>5GeV/c","l");


            leg->DrawClone();
            //sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy
    
}

 void drawAcomMulti(){
    
    readin();
    for(int i=2;i<6;i++){
    drawACompare(i);
    }
}
   
