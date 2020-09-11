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

using namespace std;

const int NBINS = 15;
const int NRT   = 4;

const float Xlow  = 1;
const float Xhigh = 250;

TFile *fin;
TFile *fin1;
TLine* line;

TH1D *htmp;
TH1D *hSumPt[9];
TH1D *hSumN[9];
TH1D *hMeanPt[9];
TH1D *hMeanN[9];

TH1D *hecc[20];

TGraphErrors *grq2[9];
TGraphErrors *gr_ratio[9];
//===============================  cumulant =====================================

TGraphErrors *grtmpErr0[100];

char name[200];

TCanvas *can[100];
TPad* pad[20][20];

TH1D *hhtem[100];
TH2D *h2D[100];

double llmrg[] = {0.2, 0.2, 0.2, 0.2};
double Lmrg[]  = {0.2,0.2};

double ratx[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
double raty[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

int mrk[] = {20,24,21,20, 33, 27};
int clr[] = {kMagenta+1,kRed+1, kBlue, 1,  kGreen+1,kGray+2, kCyan+1, kYellow+1,kViolet+1};


void readin(){
    
    //===============================  read in =====================================
    sprintf(name, "Step2/QA_afterallcuts_Sep10.root");
    
    cout<<name<<endl;
    fin = TFile::Open(name);
    
//______________________________pt hisht___________________________________---
    for(int icent=0; icent!=4; icent++){
    sprintf(name,"Pthist_Cut%d",icent);
    htmp = (TH1D*)fin->Get(name);
    if( htmp==0 ) cout<<name<<endl;
    hecc[icent] = htmp;
    
    

	
    }
  


//===============================  zpct0 =====================================
    cout<<"finish read in "<<endl;
//=============================== finish read in =============================
    
}

void drawGrFit(){
    
    int nx = 1;
    int ny = 1;
    
    int icent = 8;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.04;
    
    //sprintf(name, "GrFit_cent%d", icent);
    
    can[0] = newDivCan2( "0-20%_jetPt", Lmrg,llmrg, ratx,  raty, nx, ny, 400, 400 );
    
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
            
            char *ytile ="counts";
            char *xtile ="p_{T}^{track} (GeV/c)";
            
            hhtem[1] =  (TH1D*)    hecc[0] ->Clone();
            hhtem[2] =  (TH1D*)    hecc[1] ->Clone();
            hhtem[3] =  (TH1D*)    hecc[2] ->Clone();
            hhtem[4] =  (TH1D*)    hecc[3] ->Clone();

            
	hhtem[1]->SetLineColor(2);
  hhtem[1]->SetLineWidth(2);
//	hhtem[1]->SetMarkerSize(0.6);

           //hhtem[2]->SetMarkerStyle(20);
           //hhtem[2]->SetMarkerColor(4);
        hhtem[2]->SetLineColor(4);
        hhtem[2]->SetLineWidth(2);
//	hhtem[2]->SetMarkerSize(0.6);

            //hhtem[3]->SetMarkerStyle(21);
            //hhtem[3]->SetMarkerColor(6);
         hhtem[3]->SetLineColor(6);
         hhtem[3]->SetLineWidth(2);
//	hhtem[3]->SetMarkerSize(0.6);

		//hhtem[4]->SetMarkerStyle(29);
          //  hhtem[4]->SetMarkerColor(12);
		        hhtem[4]->SetLineColor(12);
            hhtem[4]->SetLineWidth(2);





           hhtem[1] -> GetXaxis()->SetTitle(xtile);
            hhtem[1] -> GetYaxis()->SetTitle(ytile);

            hhtem[1]->GetYaxis()->SetTitleOffset(1.6);
	//hhtem[1]->GetXaxis()->SetTitleOffset(1.2);   
         hhtem[1]->GetXaxis()->SetNdivisions(505);
            hhtem[1]->GetYaxis()->SetRangeUser(1,1e12); 
            hhtem[1]->DrawClone("hist");
           /* ffit_BG->DrawClone("same");
            ffit_Gaus->DrawClone("same");
            ffit_Power->DrawClone("same");
            */
		hhtem[2]->DrawClone("same hist");
		hhtem[3]->DrawClone("same hist");
		hhtem[4]->DrawClone("same hist");
/*		hhtem[5]->DrawClone("same");            
		hhtem[6]->DrawClone("same");
*/	//	hhtem[7]->DrawClone("same");

            tx0=0.25, ty0=0.58;
           // myTextF(tx0,ty0,"primary track",tsize*0.8,1,12);
    	   //myTextF(tx0,ty0,"Dca<1&nHitsRatio>0.52&nHitsFit>15",tsize*0.8,1,12);        
/*            tx0=0.5;ty0=0.62;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T} charged Jets,R=0.2",tsize*0.8,1,12);
*/
            
            leg = mylegF(0.60,0.70,0.88,0.88,0.035);
            leg->AddEntry(hhtem[1],"primary tracks","l");
            //leg->DrawClone();
            
            //leg = mylegF(0.30,0.68,0.88,0.88,0.03);
            leg->AddEntry(hhtem[2],"|#eta|<1&nHitsFit>15&Dca<3","l");
            //leg->DrawClone();
           
            //leg = mylegF(0.30,0.48,0.88,0.88,0.03);
            leg->AddEntry(hhtem[3],"|#eta|<1&nHitsFit>15&Dca<1","l");
          //  leg->DrawClone();
            
           // leg = mylegF(0.30,0.38,0.88,0.88,0.03);
            leg->AddEntry(hhtem[4],"|#eta|<1&nHitsFit>15&Dca<1&nHitsRatio>0.52","l");
            leg->DrawClone();
    
 
            
        }//ix
    }//iy
    
    //sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
    
}



void drawpt(){
    
    readin();
    
    drawGrFit();
    
}







