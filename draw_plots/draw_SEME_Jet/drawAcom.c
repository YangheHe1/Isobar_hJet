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
TFile *fin_py;


TFile *fout;

TLine* line;

TH1D *hArea[7];
TH1D *hNtrigger[7];

TH1D *hArea_py[7];
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
    
    sprintf(name,"Harea9_30");
    hArea[i]     = (TH1D*)fin->Get(name);
    
    sprintf(name,"Htrg9_30");
    hNtrigger[i]     = (TH1D*)fin->Get(name);

    float nraw = hNtrigger[i] -> GetEntries();
    
    
    hArea[i] -> Scale(1./nraw);

    //_____________________________
    sprintf(name,"Pythia6_A/jet_pythia6_R%i.root",i);
    cout<<name<<endl;
    fin_py = TFile::Open(name);
    
    sprintf(name,"Harea9_30");
    hArea_py[i]     = (TH1D*)fin_py->Get(name);
    
    sprintf(name,"Htrg9_30");
    hNtrigger_py[i]     = (TH1D*)fin_py->Get(name);

    float npy = hNtrigger_py[i] -> GetEntries();
    
    
    hArea_py[i] -> Scale(1./npy);
    
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
            hhtem[1] =  (TH1D*)    hArea_py[_form]  ->Clone();
            


            char *xtile ="A";
            char *ytile ="N_{jet}/N_{trg}";
            
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
            myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.6,1,12);
            tx0=0.50;ty0=0.72;
            //sprintf(name,"centrality 0-20%");

            if(_form==2) myTextF(tx0,ty0,"anti_k_{T} algorithm,R=0.2",tsize*0.6,1,12);
            if(_form==3) myTextF(tx0,ty0,"anti_k_{T} algorithm,R=0.3",tsize*0.6,1,12);
            if(_form==4) myTextF(tx0,ty0,"anti_k_{T} algorithm,R=0.4",tsize*0.6,1,12);
            if(_form==5) myTextF(tx0,ty0,"anti_k_{T} algorithm,R=0.5",tsize*0.6,1,12);

	    tx0=0.50, ty0=0.65;
	    myTextF(tx0,ty0,"h^{#pm}+jet,7<P_{trg,T}<30 GeV/c",tsize*0.6,1,12);
        




    
            leg = mylegF(0.50,0.50,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"same events","l");

            leg->AddEntry(hhtem[1],"pythia6 star tune","l");

            leg->DrawClone();
            //sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy
    
}

 void drawAcom(){
    
    readin();
    for(int i=2;i<6;i++){
    drawACompare(i);
    }
}
   
