#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include<stdio.h>
#include<algorithm>
#include<vector>
#include<iostream>
#include <string>
#include <stdlib.h> // strtod(), string to double
#include "drawheader_v4.h"
#include "TGraph.h"

using std::vector;
using std::string;

using namespace std;

const int NBINS = 15;
const int NRT   = 4;

const float Xlow  = 1;
const float Xhigh = 500;

const double PI = acos(-1.0);

double Jet_R=0.2;
const int regiontype=3;
double pt_region_C[regiontype][2] = { {-4,0}, {-4,-1}, {0,2} };
double pt_region_P[regiontype][2] = { {0,3}, {0,2}, {0,0} };

TFile *fin;
TFile *fin_mix;

vector<int> pt_regionC;
vector<int> pt_regionP;

TFile *fout;

TLine* line;

//________0 central 1 peripheral_______________ 
TH1D *hJetPt_SE[2];
TH1D *hNtrigger_SE[2];


TH1D *hJetPt_ME[2];
TH1D *hNtrigger_ME[2];

TH1D *hJetPt_pythia[2];
TH1D *hNtrigger_pythia[2];

TH1D *hRatio[2];

TH1D *hRatio_pythia[2];


//_______region test________________________________________
TH1D *hJetPt_ME_scaled[2][regiontype];
TH1D *hRatio_scaled[2][regiontype];


TGraph *norm_region[2][regiontype];
TGraph *norm_region_scaled[2][regiontype];

TH1D *hJet_ME_Sub_SE[2][regiontype];
//__________________________________________________________
TLegend *leg;

TGraphErrors *grtmpErr0[100];
TGraph *gr[100];
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

const char* type[]={"same event", "mix event","mix event(scaled)","norm.region", "^{4}He+^{14}N 1 PeV", "^{4}He+^{14}N 10 PeV"};

void readin(){
    
    //===============================  read in =====================================
    sprintf(name,"jet_SE_R2.root");
    cout<<name<<endl;
    fin = TFile::Open(name);

    float nraw[2];
    for(int icent=0;icent<2;icent++){
        sprintf(name,"recoil_jet_cent%d",icent);
        hJetPt_SE[icent]     = (TH1D*)fin->Get(name);

        sprintf(name,"number_of_trigger_cent%d",icent);
        hNtrigger_SE[icent]    = (TH1D*)fin->Get(name);

        nraw[icent] = hNtrigger_SE[icent] -> GetEntries();
        hJetPt_SE[icent] -> Scale(1./nraw[icent]);

    }
    //===============================  read in =====================================
    sprintf(name,"jet_ME_R2.root");
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    float nmix[2];
    for(int icent=0;icent<2;icent++){
        sprintf(name,"recoil_jet_cent%d",icent);
        hJetPt_ME[icent]     = (TH1D*)fin->Get(name);

        sprintf(name,"number_of_trigger_cent%d",icent);
        hNtrigger_ME[icent]    = (TH1D*)fin->Get(name);

        nmix[icent] = hNtrigger_ME[icent] -> GetEntries();
        hJetPt_ME[icent] -> Scale(1./nmix[icent]);

        //hJetPt_ME_scaled[icent] = (TH1D*) hJetPt_ME[icent]->Clone();      
    }

    //===============================  read in =====================================
    sprintf(name,"jet_pythia6_Jun05.root");
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    float npy[2];
    for(int icent=0;icent<1;icent++){
        sprintf(name,"recoil_jet_cent%d",icent);
        hJetPt_pythia[icent]     = (TH1D*)fin->Get(name);

        sprintf(name,"number_of_trigger_cent%d",icent);
        hNtrigger_pythia[icent]    = (TH1D*)fin->Get(name);

        npy[icent] = hNtrigger_pythia[icent] -> GetEntries();
        hJetPt_pythia[icent] -> Scale(1./npy[icent]);
        //hJetPt_ME_scaled[icent] = (TH1D*) hJetPt_ME[icent]->Clone();      
    }

    hJetPt_pythia[1] = (TH1D*)hJetPt_pythia[0]->Clone();
    
    cout<<"mix event "<<nmix[0]<<" same event "<<nraw[0]<<endl;
    cout<<"mix eventP "<<nmix[1]<<" same eventP "<<nraw[1]<<endl;
    //_____________________________________________________________________

    for(int cent=0;cent<2;cent++){
            Int_t binmax=hJetPt_SE[cent]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_SE[cent]->GetBinContent(ibin);
                double deltap = hJetPt_SE[cent]->GetBinWidth(ibin);
                double deltaeta = 2-(2*Jet_R);
                double binerr = hJetPt_SE[cent]->GetBinError(ibin);
                hJetPt_SE[cent]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_SE[cent]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }

            binmax=hJetPt_ME[cent]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_ME[cent]->GetBinContent(ibin);
                double deltap = hJetPt_ME[cent]->GetBinWidth(ibin);
                double deltaeta = 2-(2*Jet_R);
                double binerr = hJetPt_ME[cent]->GetBinError(ibin);
                hJetPt_ME[cent]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_ME[cent]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }

            binmax=hJetPt_pythia[cent]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_pythia[cent]->GetBinContent(ibin);
                double deltap = hJetPt_pythia[cent]->GetBinWidth(ibin);
                double deltaeta = 2-(2*Jet_R);
                double binerr = hJetPt_pythia[cent]->GetBinError(ibin);
                hJetPt_pythia[cent]->SetBinContent(ibin,y/(deltap*deltaeta));
                hJetPt_pythia[cent]->SetBinError(ibin,binerr/(deltap*deltaeta));

            }

    }
    //_________________________without scale____________________________________
    
    for(int cent=0;cent<2;cent++){
        hRatio[cent] = (TH1D*) hJetPt_SE[cent]->Clone();
        Int_t binmax=hJetPt_SE[cent]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE[cent]->GetBinContent(ibin);
            double yref = hJetPt_ME[cent]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio[cent]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE[cent]->GetBinError(ibin);
            double binerref =hJetPt_ME[cent]->GetBinError(ibin);
            hRatio[cent]->SetBinContent(ibin,y/yref);
            double err = (binerr/yref)-(y*binerref/(yref*yref));
            hRatio[cent]->SetBinError(ibin,err);

        }

        /*
        hRatio_pythia[cent] = (TH1D*) hJetPt_SE[cent]->Clone();
        Int_t binmax=hJetPt_SE[cent]->FindLastBinAbove(0,1);
        for(int ibin=1;ibin<=binmax;ibin++){
            double y = hJetPt_SE[cent]->GetBinContent(ibin);
            double yref = hJetPt_pythia[cent]->GetBinContent(ibin);
            if(y==0||yref==0)
			    {   hRatio_pythia[cent]->SetBinContent(ibin,0);
				    continue;}
            double binerr = hJetPt_SE[cent]->GetBinError(ibin);
            double binerref =hJetPt_pythia[cent]->GetBinError(ibin);
            hRatio_pythia[cent]->SetBinContent(ibin,y/yref);
            double err = (binerr/yref)-(y*binerref/(yref*yref));
            hRatio_pythia[cent]->SetBinError(ibin,err);

        }
        */
    
    }

    //______________calculate left region_______________________

    //________________calculat scaling factor__________________________
    double scaling_factor[2][regiontype];

    
    for(int cent=0;cent<2;cent++){

        for(int reg=0;reg<regiontype;reg++){

            norm_region[cent][reg] =new TGraph();
            sprintf(name,"norm_region_cent%d_reg%d",cent,reg);
            norm_region[cent][reg]->SetName(name);
            Int_t i_point_norm = 0;

            int region_begin;
            int region_stop;
            double integral_ME=0;
            double integral_SE=0;
            if(cent==0){
                region_begin= pt_region_C[reg][0];
                region_stop= pt_region_C[reg][1];
            }
            if(cent==1){
                region_begin= pt_region_P[reg][0];
                region_stop= pt_region_P[reg][1];
            }

            Int_t binmax=hJetPt_SE[cent]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<binmax;ibin++){
                double y = hJetPt_SE[cent]->GetBinContent(ibin);
                double y_ME = hJetPt_ME[cent]->GetBinContent(ibin);
                double x = hJetPt_SE[cent]->GetBinCenter(ibin);
                if(x<region_begin) continue;
                if(x>region_stop)  break;
                integral_SE += y;
                integral_ME += y_ME;

                double x_val_ME = hJetPt_ME[cent]->GetBinCenter(ibin);
                if(i_point_norm == 0) norm_region[cent][reg] ->SetPoint(i_point_norm,x_val_ME,0.0);
                norm_region[cent][reg] ->SetPoint(i_point_norm+1,x_val_ME,y_ME);
                norm_region[cent][reg] ->SetPoint(i_point_norm+2,x_val_ME,0.0);

                i_point_norm++;            
        }

            if(integral_ME>0) scaling_factor[cent][reg]=integral_SE/integral_ME;
        //___________________________________________________________________________________________

        }//regiontype        
    
    }//cent

	cout<<"scaling factor: CENTRAL  "<<scaling_factor[0][0]<<" PERIPHERAL "<<scaling_factor[1][0]<<endl;
    cout<<"scaling factor 1: CENTRAL  "<<scaling_factor[0][1]<<" PERIPHERAL "<<scaling_factor[1][1]<<endl;
    //_____________ME sclaing and ratio_____________________
    for(int cent=0;cent<2;cent++){

        for(int reg=0;reg<regiontype;reg++){

	        hJetPt_ME_scaled[cent][reg] = (TH1D*) hJetPt_ME[cent]->Clone(); 
            hJetPt_ME_scaled[cent][reg]->Scale(scaling_factor[cent][reg]);

            hRatio_scaled[cent][reg] = (TH1D*) hJetPt_SE[cent]->Clone();
            Int_t binmax=hJetPt_SE[cent]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_SE[cent]->GetBinContent(ibin);
                double yref = hJetPt_ME_scaled[cent][reg]->GetBinContent(ibin);
                if(y==0||yref==0){
			        hRatio_scaled[cent][reg]->SetBinContent(ibin,0);
				    continue;}
                double binerr = hJetPt_SE[cent]->GetBinError(ibin);
                double binerref =hJetPt_ME_scaled[cent][reg]->GetBinError(ibin);
                hRatio_scaled[cent][reg]->SetBinContent(ibin,y/yref);
                double err = (binerr/yref)-(y*binerref/(yref*yref));
                hRatio_scaled[cent][reg]->SetBinError(ibin,err);

            }

        }


        /*
        //norm region
        norm_region_scaled[cent] =new TGraph();
        sprintf(name,"norm_region_scaled_cent%d",cent);
        norm_region_scaled[cent]->SetName(name);
        int region_begin;
        int region_stop;
        if(cent==0){
            //region_begin= pt_regionC.at(0);
            region_stop= pt_regionC.size();
        }
        if(cent==1){
            //region_begin= pt_regionP.at(0);
            region_stop= pt_regionP.size();
        }
	//region_begin=2*region_begin+1;
	    region_stop=2*region_stop+2;
        for(int ibin=0;ibin<region_stop;ibin++){
            Double_t x_val_ME, y_val_ME;
            norm_region[cent]->GetPoint(ibin,x_val_ME,y_val_ME);
            norm_region_scaled[cent]->SetPoint(ibin,x_val_ME,y_val_ME*scaling_factor[cent]);
	    //cout<<"ibin: "<<ibin<<"x_val_ME: "<<x_val_ME<<"y_val_ME "<<y_val_ME*scaling_factor[cent]<<endl;
        }
        */
        
    }

    //____ME-SE jet___________________________________________
    for(int cent=0;cent<2;cent++){

        for(int reg=0;reg<regiontype;reg++){
            hJet_ME_Sub_SE[cent][reg] = (TH1D*) hJetPt_SE[cent]->Clone();
            Int_t binmax=hJetPt_SE[cent]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_SE[cent]->GetBinContent(ibin);
                double yref = hJetPt_ME_scaled[cent][reg]->GetBinContent(ibin);
                if(y-yref<0)
			        {hJet_ME_Sub_SE[cent][reg]->SetBinContent(ibin,0);
				    continue;}
                double binerr = hJetPt_SE[cent]->GetBinError(ibin);
                double binerref =hJetPt_ME_scaled[cent][reg]->GetBinError(ibin);
                hJet_ME_Sub_SE[cent][reg]->SetBinContent(ibin,y-yref);
                double err = binerr+binerref;
                hJet_ME_Sub_SE[cent][reg]->SetBinError(ibin,err);

            }

        }
    
    }



}


void drawPtCompare(int _form){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "PtCompar_%d",_form);
    
    can[0]= new TCanvas(name,"Graph",10,10,1100,900);    
	pad[0][0]=new TPad("pad1","pad1",0.06,0.4,0.94,0.94);
    pad[0][0]->Draw();
	pad[0][1]=new TPad("pad1","pad1",0.06,0.06,0.94,0.4);
    pad[0][1]->Draw();
    
    pad[0][0]->SetTopMargin(0.08);
	pad[0][0]->SetBottomMargin(0);
    pad[0][0]->cd();
            
            
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);      //remove the entries,mean,RMS in the upper right.
    gStyle->SetOptTitle(0);
            
    hhtem[0] =  (TH1D*)    hJetPt_SE[_form]  ->Clone();
    hhtem[1] =  (TH1D*)    hJetPt_ME[_form]  ->Clone();
    //hhtem[2] =  (TH1D*)    hJetPt_pythia[_form]  ->Clone();
    gr[0]    =  (TGraph*)  norm_region[_form][0] ->Clone();
    gr[1]    =  (TGraph*)  norm_region[_form][1] ->Clone();       
    
    hhtem[0]->SetMarkerStyle(29);
    hhtem[0]->SetMarkerSize(1);
    hhtem[0]->SetMarkerColor(2);
    hhtem[0]->SetLineColor(2);

    //hhtem[1]->SetMarkerStyle(20);
    //hhtem[1]->SetMarkerSize(0.7);
    //hhtem[1]->SetMarkerColor(4);
    hhtem[1]->SetLineColor(kGray+2);
    hhtem[1]->SetLineWidth(1.5);
    hhtem[1]->SetLineStyle(2);
    hhtem[1]->SetFillColor(kGray+2);
    hhtem[1]->SetFillStyle(3004);

    /*
    hhtem[2]->SetMarkerStyle(25);
    hhtem[2]->SetMarkerSize(0.7);
    hhtem[2]->SetMarkerColor(kGreen+4);
    hhtem[2]->SetLineColor(kGreen+4);
    */

    gr[0]->SetLineStyle(2);
    gr[0]->SetFillColor(4);
    gr[0]->SetFillStyle(3005);

    gr[1]->SetLineStyle(2);
    gr[1]->SetFillColor(kGreen+2);
    gr[1]->SetFillStyle(3002);

            
            //char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="(1/N_{trig})d^{2}N_{jets}/(dp_{T,jet}^{reco,ch}d#eta) (Gev/c)^{-1} ";
            
            //hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetRangeUser(-5,30);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.2);
            hhtem[1]->SetTitleSize(0.05);
            hhtem[0]->GetYaxis()->CenterTitle(true);
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("same hist");
            //hhtem[2]->DrawClone("same flx");
 //           gr[0]->DrawClone("same F");
 //           gr[1]->DrawClone("same F");
            

            //hhtem[1]->DrawClone("same");
            
            tx0=0.50, ty0=0.84;
            if(_form==0) myTextF(tx0,ty0,"Isobar 200GeV, 0-10%",tsize*1.,1,12);
            if(_form==1) myTextF(tx0,ty0,"Isobar 200GeV, 60-80%",tsize*1.,1,12);
            
            tx0=0.50;ty0=0.76;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, R=0.2",tsize*1.,1,12);

	        tx0=0.50, ty0=0.69;
	        myTextF(tx0,ty0,"h^{#pm}+jet,7<P_{T}^{trg}<30GeV/c",tsize*1.,1,12);
            
            float _yy = 0.60;
            leg = mylegF(0.6,_yy,0.8,0.88,0.05);
            
            leg->AddEntry(hhtem[0],type[0],"pl");
	        leg->AddEntry(hhtem[1],"mix event","l");
	    //leg->AddEntry(hhtem[2],"Pythia6 STAR tune","l");
//	        leg->AddEntry(gr[0],"Fill region range1","f");
//            leg->AddEntry(gr[1],"Fill region range2","f");
            leg->Draw("same");
            
            
            
    pad[0][1]->SetTopMargin(0.00);
	pad[0][1]->SetBottomMargin(0.3);

    pad[0][1]->cd();

    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);    

    char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
    char *ytile1 ="SE/ME";

    hhtem[3] =  (TH1D*)    hRatio[_form] ->Clone();
    //hhtem[4] =  (TH1D*)    hRatio_pythia[_form] ->Clone();
    hhtem[4] =  (TH1D*)    hRatio_scaled[_form][0] ->Clone();
    hhtem[5] =  (TH1D*)    hRatio_scaled[_form][1] ->Clone();

    hhtem[3]->SetMarkerStyle(29);
    hhtem[3]->SetMarkerSize(0.7);
    hhtem[3]->SetMarkerColor(2);
    hhtem[3]->SetLineColor(2);

    hhtem[4]->SetMarkerStyle(29);
    hhtem[4]->SetMarkerSize(0.7);
    hhtem[4]->SetMarkerColor(4);
    hhtem[4]->SetLineColor(4);


    hhtem[5]->SetMarkerStyle(29);
    hhtem[5]->SetMarkerSize(0.7);
    hhtem[5]->SetMarkerColor(kGreen+2);
    hhtem[5]->SetLineColor(kGreen+2);

    hhtem[3] -> GetXaxis()->SetTitle(xtile);
    hhtem[3] -> GetYaxis()->SetTitle(ytile1);
    hhtem[3] -> GetXaxis()->SetNdivisions(507);
    hhtem[3] -> GetYaxis()->SetNdivisions(507);
    hhtem[3] -> GetXaxis()->SetRangeUser(-5,30);
            
    hhtem[3] -> GetYaxis()->SetTitleOffset(1.2); 
    hhtem[3]->GetXaxis()->SetTitleSize(0.08);
    hhtem[3]->GetYaxis()->CenterTitle(true);           
    hhtem[3]->DrawClone("P");
 //   hhtem[4]->DrawClone("Psame");
 //   hhtem[5]->DrawClone("Psame");


    TLine *tl = new TLine(-5,1,30,1);
	tl->Draw("same");    
            
   
    
}

/*
void drawPtScaledCompare(int _form){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "PtScaaledCompar_%d",_form);
    
    can[0]= new TCanvas(name,"Graph",10,10,1100,900);    
	pad[0][0]=new TPad("pad1","pad1",0.06,0.4,0.94,0.94);
    pad[0][0]->Draw();
	pad[0][1]=new TPad("pad1","pad1",0.06,0,0.94,0.4);
    pad[0][1]->Draw();
    
    pad[0][0]->SetTopMargin(0.08);
	pad[0][0]->SetBottomMargin(0);
    pad[0][0]->cd();
            
            
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLogy(1);
    gStyle->SetOptStat(0);      //remove the entries,mean,RMS in the upper right.
    gStyle->SetOptTitle(0);
            
    hhtem[0] =  (TH1D*)    hJetPt_SE[_form]  ->Clone();
    hhtem[1] =  (TH1D*)    hJetPt_ME_scaled[_form]  ->Clone();
    gr[0]    =  (TGraph*)  norm_region_scaled[_form] ->Clone();      
    
    hhtem[0]->SetMarkerStyle(29);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetLineColor(2);
    hhtem[0]->SetMarkerColor(2);

    hhtem[1]->SetMarkerStyle(20);
    hhtem[1]->SetMarkerSize(0.7);
    hhtem[1]->SetMarkerColor(4);
    hhtem[1]->SetLineColor(4);
    //hhtem[1]->SetLineWidth(2);
    //hhtem[1]->SetLineStyle(1);
    hhtem[1]->SetFillColor(kGray+2);
    hhtem[1]->SetFillStyle(3004);

    gr[0]->SetLineStyle(2);
    gr[0]->SetFillColor(kAzure+2);
    gr[0]->SetFillStyle(3004);

            
            //char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="(1/N_{trig})d^{2}N_{jets}/(dp_{T,jet}^{reco,ch}d#eta) (Gev/c)^{-1} ";
            
            //hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetRangeUser(-10,35);
            //hhtem[0] -> GetYaxis()->SetRangeUser(-10,35);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.2);
            hhtem[1]->SetTitleSize(0.05);
            hhtem[0]->GetYaxis()->CenterTitle(true);
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("same Pflx");
            gr[0]->DrawClone("same Flx");
            

            //hhtem[1]->DrawClone("same");
            
            tx0=0.50, ty0=0.84;
            if(_form==0) myTextF(tx0,ty0,"Isobar 200GeV, 0-10%",tsize*1.,1,12);
            if(_form==1) myTextF(tx0,ty0,"Isobar 200GeV, 60-80%",tsize*1.,1,12);
            
            tx0=0.50;ty0=0.76;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, R=0.2",tsize*1.,1,12);

	        tx0=0.50, ty0=0.69;
	        myTextF(tx0,ty0,"h^{#pm}+jet",tsize*1.,1,12);
            
		float _yy = 0.60;
            leg = mylegF(0.6,_yy,0.8,0.88,0.05);

            leg->AddEntry(hhtem[0],type[0],"lp");
            leg->AddEntry(hhtem[1],type[2],"lp");
            //leg->AddEntry(hhtem[2],type[1],"l");
            leg->AddEntry(gr[0],type[3],"l");
            leg->Draw("same");
            
            
    pad[0][1]->SetTopMargin(0.00);
	pad[0][1]->SetBottomMargin(0.3);

    pad[0][1]->cd();

    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLogy(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);    

    char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
    char *ytile1 ="SE/ME";

    hhtem[3] =  (TH1D*)    hRatio_scaled[_form] ->Clone();

    hhtem[3]->SetMarkerStyle(29);
    hhtem[3]->SetMarkerSize(0.7);
    hhtem[3]->SetMarkerColor(2);
    hhtem[3]->SetLineColor(2);

    hhtem[3] -> GetXaxis()->SetTitle(xtile);
    hhtem[3] -> GetYaxis()->SetTitle(ytile1);
    hhtem[3] -> GetXaxis()->SetNdivisions(507);
    hhtem[3] -> GetYaxis()->SetNdivisions(507);
    hhtem[3] -> GetXaxis()->SetRangeUser(-10,35);
            
    hhtem[3] -> GetYaxis()->SetTitleOffset(1.2); 
    hhtem[3]->GetXaxis()->SetTitleSize(0.08);
    hhtem[3]->GetYaxis()->CenterTitle(true);           
    hhtem[3]->DrawClone("P");

    TLine *tl = new TLine(-10,1,35,1);
	tl->Draw("same");    
            
   
    
}
*/
void drawSESubME(){
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "SESubME");
    
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
            
            hhtem[0] =  (TH1D*)    hJet_ME_Sub_SE[0][1]  ->Clone();
            hhtem[1] =  (TH1D*)    hJet_ME_Sub_SE[1][1]  ->Clone();
            hhtem[2] =  (TH1D*)    hJetPt_pythia[0]  ->Clone();
            
            
            
            
            char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="(1/N_{trig})d^{2}N_{jets}/(dp_{T,jet}^{reco,ch}d#eta) (Gev/c)^{-1} ";
            
    hhtem[0]->SetMarkerStyle(25);
    hhtem[0]->SetMarkerSize(0.7);
    hhtem[0]->SetMarkerColor(1);
    hhtem[0]->SetLineColor(1);

    hhtem[1]->SetMarkerStyle(24);
    hhtem[1]->SetMarkerSize(0.7);
    hhtem[1]->SetMarkerColor(kRed);
    hhtem[1]->SetLineColor(kRed);

    hhtem[2]->SetMarkerStyle(20);
    hhtem[2]->SetMarkerSize(0.7);
    hhtem[2]->SetMarkerColor(kGreen+2);
    hhtem[2]->SetLineColor(kGreen+2);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetXaxis()->SetRangeUser(0,1100);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,1);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
           
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("Psame");
            hhtem[2]->DrawClone("Psame");
            
            tx0=0.50, ty0=0.84;
            myTextF(tx0,ty0,"Isobar,200GeV",tsize*0.8,1,12);
            tx0=0.50;ty0=0.76;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, R=0.2",tsize*0.8,1,12);

	    tx0=0.50, ty0=0.69;
	    myTextF(tx0,ty0,"h^{#pm}+jet",tsize*0.8,1,12);


            leg = mylegF(0.50,0.45,0.65,0.60,0.03);
            leg->AddEntry(hhtem[0],"0-10%","lp");

            leg->AddEntry(hhtem[1],"60-80%","lp");

            leg->AddEntry(hhtem[2],"Pythia6 STAR tune","lp");
            
            leg->DrawClone();
            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy

}

 void drawJetSub(){
    
    readin();
    for(int i=0;i<2;i++){
    drawPtCompare(i);
    //drawPtScaledCompare(i);
    
    }
 //   drawSESubME();
    
}
   
