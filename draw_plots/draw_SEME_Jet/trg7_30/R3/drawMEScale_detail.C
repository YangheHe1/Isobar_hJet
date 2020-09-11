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




TFile *fin;
TFile *fin_mix;
TFile *fout;
TLine* line;

// settings____________________________________
double Jet_R=0.3;

//__normalization range_________________________
double regionP[2]={0,0};
double regionC[2]={0,-2};

//________0 central 1 peripheral_______________ 
TH1D *hJetPt_SE[2];
TH1D *hNtrigger_SE[2];


TH1D *hJetPt_ME[2];
TH1D *hNtrigger_ME[2];

TH1D *hJetPt_pythia[2];
TH1D *hNtrigger_pythia[2];

TH1D *hRatio[2];

TH1D *hRatio_pythia[2];

TH1D *hratio_for_fit[2];


//_______region test________________________________________
TH1D *hJetPt_ME_scaled[2];
TH1D *hRatio_scaled[2];


TGraph *gJetPt_ME_scaled[2];
TGraph *norm_region[2];
TGraph *norm_region_scaled[2];

TH1D *hJet_ME_Sub_SE[2];
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
    sprintf(name,"jet_SE_R3.root");
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
    sprintf(name,"jet_ME_R3.root");
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
    sprintf(name,"jet_Pythia6_R3.root");
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

        
    
    }

    //______________calculate left region_______________________


    //________________calculat scaling factor__________________________
    double scaling_factor[2];
    
   
    for(int cent=0;cent<2;cent++){

        hratio_for_fit[cent] = (TH1D*)hRatio[cent]->Clone();
        int fbin = hratio_for_fit[cent]->FindFirstBinAbove(0,1);
        int ran = (int) hratio_for_fit[cent]->GetBinCenter(fbin);
        int ran_f = ran-1;
        int ran_l;
        if(cent==0) {
            ran_l = regionC[1];
            regionC[0] = ran_f;
            //regionC[1] = ran_l;

        }


        if(cent==1){
            ran_l = regionP[1];
            regionP[0] = ran_f;
            //regionP[1] = ran_l;

        }
        
        //if(cent==1) continue;

        TF1 *f1= new TF1("f1","[0]",ran_f,ran_l);
        hratio_for_fit[cent]->Fit("f1","R");

        
        double fit_par = f1->GetParameter(0);

        scaling_factor[cent]= fit_par;
        //___________________________________________________________________________________________

           
    
    }//cent

    scaling_factor[1] = 1;

	cout<<"scaling factor: CENTRAL  "<<scaling_factor[0]<<" PERIPHERAL "<<scaling_factor[1]<<endl;
   
    //_____________ME sclaing and ratio_____________________
    for(int cent=0;cent<2;cent++){

        

	    hJetPt_ME_scaled[cent] = (TH1D*) hJetPt_ME[cent]->Clone(); 
        hJetPt_ME_scaled[cent]->Scale(scaling_factor[cent]);

        hRatio_scaled[cent] = (TH1D*) hJetPt_SE[cent]->Clone();
        Int_t binmax=hJetPt_SE[cent]->FindLastBinAbove(0,1);
            for(int ibin=1;ibin<=binmax;ibin++){
                double y = hJetPt_SE[cent]->GetBinContent(ibin);
                double yref = hJetPt_ME_scaled[cent]->GetBinContent(ibin);
                if(y==0||yref==0){
			        hRatio_scaled[cent]->SetBinContent(ibin,0);
				    continue;}
                double binerr = hJetPt_SE[cent]->GetBinError(ibin);
                double binerref =hJetPt_ME_scaled[cent]->GetBinError(ibin);
                hRatio_scaled[cent]->SetBinContent(ibin,y/yref);
                double err = (binerr/yref)-(y*binerref/(yref*yref));
                hRatio_scaled[cent]->SetBinError(ibin,err);    

        }

        
        //norm region
        //if(cent==1) continue;

        norm_region_scaled[cent] =new TGraph();
        sprintf(name,"norm_region_scaled_cent%d",cent);
        norm_region_scaled[cent]->SetName(name);
        Int_t i_point_norm = 0;
/*
        gJetPt_ME_scaled[cent] =new TGraph();
        sprintf(name,"gJetPt_ME_scaled_cent%d",cent);
        gJetPt_ME_scaled[cent]->SetName(name);
        Int_t g_point_norm = 0;
*/
        int region_begin=0;
        int region_stop=0;
        
        if(cent==0){
            region_begin = regionC[0];
            region_stop= regionC[1];
        }
        
        if(cent==1){
            region_begin = regionP[0];
            region_stop= regionP[1];
        }
        

        Int_t r_binmax=hJetPt_ME_scaled[cent]->FindLastBinAbove(0,1);
        /*    for(int ibin=1;ibin<r_binmax;ibin++){
                double y_ME = hJetPt_ME_scaled[cent]->GetBinContent(ibin);
                double x_val_ME = hJetPt_ME_scaled[cent]->GetBinCenter(ibin);

                if(g_point_norm == 0) gJetPt_ME_scaled[cent] ->SetPoint(g_point_norm,x_val_ME,0.0);
                gJetPt_ME_scaled[cent] ->SetPoint(g_point_norm+1,x_val_ME,y_ME);
                gJetPt_ME_scaled[cent] ->SetPoint(g_point_norm+2,x_val_ME,0.0);

                g_point_norm++;
            }

        */
            for(int ibin=1;ibin<r_binmax;ibin++){
                double y_ME = hJetPt_ME_scaled[cent]->GetBinContent(ibin);
                double x_val_ME = hJetPt_ME_scaled[cent]->GetBinCenter(ibin);

                if(x_val_ME<region_begin) continue;
                if(x_val_ME>region_stop)  break;
                
                if(i_point_norm == 0) norm_region_scaled[cent] ->SetPoint(i_point_norm,x_val_ME,0.0);
                norm_region_scaled[cent] ->SetPoint(i_point_norm+1,x_val_ME,y_ME);
                norm_region_scaled[cent] ->SetPoint(i_point_norm+2,x_val_ME,0.0);

                i_point_norm++;            
        }

        //_______norm region_______________________________________
        norm_region[cent] =new TGraph();
        sprintf(name,"norm_region_cent%d",cent);
        norm_region[cent]->SetName(name);
        Int_t b_point_norm = 0;

         for(int j=1;j<r_binmax;j++){
                double y_ME = hJetPt_ME[cent]->GetBinContent(j);
                double x_val_ME = hJetPt_ME[cent]->GetBinCenter(j);

                if(x_val_ME<region_begin) continue;
                if(x_val_ME>region_stop)  break;
                
                if(b_point_norm == 0) norm_region[cent] ->SetPoint(b_point_norm,x_val_ME,0.0);
                norm_region[cent] ->SetPoint(b_point_norm+1,x_val_ME,y_ME);
                norm_region[cent] ->SetPoint(b_point_norm+2,x_val_ME,0.0);

                b_point_norm++;            
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
    hhtem[1] =  (TH1D*)    hJetPt_ME_scaled[_form]  ->Clone();
    //gr[1]    =  (TGraph*)  gJetPt_ME_scaled[_form] ->Clone();
    //hhtem[2] =  (TH1D*)    hJetPt_pythia[_form]  ->Clone();
    gr[0]    =  (TGraph*)  norm_region_scaled[_form] ->Clone();      
    
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
    gr[1]->SetLineColor(kGray+2);
    gr[1]->SetLineWidth(1.5);
    gr[1]->SetLineStyle(2);
    gr[1]->SetFillColor(kGray+2);
    gr[1]->SetFillStyle(3004);
*/    
    
    gr[0]->SetLineStyle(2);
    gr[0]->SetFillColor(kGreen+2);
    gr[0]->SetFillStyle(3002);;
    
    

            
            //char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="(1/N_{trig})d^{2}N_{jets}/(dp_{T,jet}^{reco,ch}d#eta) (Gev/c)^{-1} ";
            
            //hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            if(_form==0) hhtem[0] -> GetXaxis()->SetRangeUser(-7,10);
            if(_form==1) hhtem[0] -> GetXaxis()->SetRangeUser(-2,10);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            
            hhtem[0] -> GetYaxis()->SetLabelSize(0.05);
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.0);
            hhtem[0] -> GetYaxis()->SetTitleSize(0.05);
            hhtem[0]->GetYaxis()->CenterTitle(true);
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("hist same");
            //gr[1]->DrawClone("same F2");
            //hhtem[2]->DrawClone("same flx");
            gr[0]->DrawClone("same F");
 
            

            //hhtem[1]->DrawClone("same");
            
            tx0=0.50, ty0=0.84;
            if(_form==0) myTextF(tx0,ty0,"Isobar 200GeV, 0-10%",tsize*1.,1,12);
            if(_form==1) myTextF(tx0,ty0,"Isobar 200GeV, 60-80%",tsize*1.,1,12);
            
            tx0=0.50;ty0=0.76;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, R=0.3",tsize*1.,1,12);

	        tx0=0.50, ty0=0.69;
	        myTextF(tx0,ty0,"h^{#pm}+jet,7<P_{T}^{trig}<30GeV/c",tsize*1.,1,12);
            
            float _yy = 0.38;
            leg = mylegF(0.6,_yy,0.8,0.58,0.05);
            
            leg->AddEntry(hhtem[0],type[0],"pl");
	        leg->AddEntry(hhtem[1],"mix events(scaled)","l");
	    //leg->AddEntry(hhtem[2],"Pythia6 STAR tune","l");
	        leg->AddEntry(gr[0],"norm. region","f");
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

    hhtem[3] =  (TH1D*)    hRatio_scaled[_form] ->Clone();
    

    hhtem[3]->SetMarkerStyle(29);
    hhtem[3]->SetMarkerSize(0.7);
    hhtem[3]->SetMarkerColor(2);
    hhtem[3]->SetLineColor(2);

    

    hhtem[3] -> GetXaxis()->SetTitle(xtile);
    hhtem[3] -> GetYaxis()->SetTitle(ytile1);
    hhtem[3] -> GetXaxis()->SetNdivisions(507);
    hhtem[3] -> GetYaxis()->SetNdivisions(507);
    if(_form==0) hhtem[3] -> GetXaxis()->SetRangeUser(-7,10);
    if(_form==1) hhtem[3] -> GetXaxis()->SetRangeUser(-2,10);

    hhtem[3] -> GetYaxis()->SetTitleOffset(0.7); 
    hhtem[3] -> GetYaxis()->SetLabelSize(0.07);
    hhtem[3] -> GetXaxis()->SetLabelSize(0.07);
    hhtem[3]->GetYaxis()->SetTitleSize(0.07);             
    hhtem[3]->GetXaxis()->SetTitleSize(0.08);
    hhtem[3]->GetYaxis()->CenterTitle(true);           
    hhtem[3]->DrawClone("P");
 


    TLine *t1;
    if(_form==0) t1= new TLine(-7,1,11,1);
    if(_form==1) t1= new TLine(-2,1,11,1);
    t1->SetLineStyle(2);
	t1->Draw("same");  



    //______________________zoom in panel
    if(_form==0) pad[0][2]=new TPad("pad1","pad1",0.12,0.54,0.42,0.99);
    if(_form==1) pad[0][2]=new TPad("pad1","pad1",0.12,0.54,0.42,0.99);
    
    pad[0][2]->Draw();

    pad[0][2]->SetBottomMargin(0.35);
    pad[0][2]->SetLeftMargin(0.15);
    pad[0][2]->cd();

    gPad->SetFillStyle(4000); // make it transparent
    
    gPad->SetTickx(1);
    gPad->SetTicky(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);    

    char *xtile ="p_{T,jet}^{reco,ch} [Gev/c]";
    char *ytile1 ="SE/ME";

    hhtem[4] =  (TH1D*)    hRatio_scaled[_form] ->Clone();  

    hhtem[4]->SetMarkerStyle(29);
    hhtem[4]->SetMarkerSize(0.7);
    hhtem[4]->SetMarkerColor(2);
    hhtem[4]->SetLineColor(2);

    

    hhtem[4] -> GetXaxis()->SetTitle(xtile);
    hhtem[4] -> GetYaxis()->SetTitle(ytile1);
    hhtem[4] -> GetXaxis()->SetNdivisions(507);
    hhtem[4] -> GetYaxis()->SetNdivisions(507);
    if(_form==0)hhtem[4] -> GetXaxis()->SetRangeUser(regionC[0],regionC[1]-1);
    if(_form==1)hhtem[4] -> GetXaxis()->SetRangeUser(regionP[0],regionP[1]-1);

    hhtem[4] -> GetYaxis()->SetTitleOffset(0.6); 
    hhtem[4] -> GetXaxis()->SetTitleOffset(1.1);
    hhtem[4]->GetXaxis()->SetTitleSize(0.12);
    hhtem[4]->GetYaxis()->SetTitleSize(0.12);
    hhtem[4] -> GetYaxis()->SetLabelSize(0.11);
    hhtem[4] -> GetXaxis()->SetLabelSize(0.12);
    hhtem[4]->GetYaxis()->CenterTitle(true);  
    hhtem[4]->GetXaxis()->CenterTitle(true);          
    hhtem[4]->DrawClone("P");

    TLine *t2;
    if(_form==0) t2 = new TLine(regionC[0],1,regionC[1],1);
    if(_form==1) {
        t2 = new TLine(regionP[0],1,regionP[1],1);
        if(regionP[0]==0) t2 = new TLine(regionP[0],1,1,1);}
    t2->SetLineStyle(2);
	t2->Draw("same");
            
   
    
}

void drawPt_beforeScaled_Compare(int _form){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "Pt_beforeScaled_Compar_%d",_form);
    
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
    //gr[1]    =  (TGraph*)  gJetPt_ME_scaled[_form] ->Clone();
    //hhtem[2] =  (TH1D*)    hJetPt_pythia[_form]  ->Clone();
    gr[0]    =  (TGraph*)  norm_region[_form] ->Clone();      
    
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
    gr[1]->SetLineColor(kGray+2);
    gr[1]->SetLineWidth(1.5);
    gr[1]->SetLineStyle(2);
    gr[1]->SetFillColor(kGray+2);
    gr[1]->SetFillStyle(3004);
*/    
    
    gr[0]->SetLineStyle(2);
    gr[0]->SetFillColor(kGreen+2);
    gr[0]->SetFillStyle(3002);;
    
    

            
            //char *xtile ="p_{T,jet}^{reco,ch}( =p_{T,jet}^{raw,ch}-#rhoA ) [Gev/c]";
            char *ytile ="(1/N_{trig})d^{2}N_{jets}/(dp_{T,jet}^{reco,ch}d#eta) (Gev/c)^{-1} ";
            
            //hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            if(_form==0) hhtem[0] -> GetXaxis()->SetRangeUser(-7,10);
            if(_form==1) hhtem[0] -> GetXaxis()->SetRangeUser(-2,10);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            
            hhtem[0] -> GetYaxis()->SetLabelSize(0.05);
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.0);
            hhtem[0] -> GetYaxis()->SetTitleSize(0.05);
            hhtem[0]->GetYaxis()->CenterTitle(true);
            
            
            hhtem[0]->DrawClone("P");
            hhtem[1]->DrawClone("hist same");
            //gr[1]->DrawClone("same F2");
            //hhtem[2]->DrawClone("same flx");
            gr[0]->DrawClone("same F");
 
            

            //hhtem[1]->DrawClone("same");
            
            tx0=0.50, ty0=0.84;
            if(_form==0) myTextF(tx0,ty0,"Isobar 200GeV, 0-10%",tsize*1.,1,12);
            if(_form==1) myTextF(tx0,ty0,"Isobar 200GeV, 60-80%",tsize*1.,1,12);
            
            tx0=0.50;ty0=0.76;
            //sprintf(name,"centrality 0-20%");
            myTextF(tx0,ty0,"Anti-k_{T}, R=0.3",tsize*1.,1,12);

	        tx0=0.50, ty0=0.69;
	        myTextF(tx0,ty0,"h^{#pm}+jet,7<P_{T}^{trig}<30GeV/c",tsize*1.,1,12);
            
            float _yy = 0.38;
            leg = mylegF(0.6,_yy,0.8,0.58,0.05);
            
            leg->AddEntry(hhtem[0],type[0],"pl");
	        leg->AddEntry(hhtem[1],"mix events","l");
	    //leg->AddEntry(hhtem[2],"Pythia6 STAR tune","l");
	        leg->AddEntry(gr[0],"norm. region","f");
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
    

    hhtem[3]->SetMarkerStyle(29);
    hhtem[3]->SetMarkerSize(0.7);
    hhtem[3]->SetMarkerColor(2);
    hhtem[3]->SetLineColor(2);

    

    hhtem[3] -> GetXaxis()->SetTitle(xtile);
    hhtem[3] -> GetYaxis()->SetTitle(ytile1);
    hhtem[3] -> GetXaxis()->SetNdivisions(507);
    hhtem[3] -> GetYaxis()->SetNdivisions(507);
    if(_form==0) hhtem[3] -> GetXaxis()->SetRangeUser(-7,10);
    if(_form==1) hhtem[3] -> GetXaxis()->SetRangeUser(-2,10);

    hhtem[3] -> GetYaxis()->SetTitleOffset(0.7); 
    hhtem[3] -> GetYaxis()->SetLabelSize(0.07);
    hhtem[3] -> GetXaxis()->SetLabelSize(0.07);
    hhtem[3]->GetYaxis()->SetTitleSize(0.07);             
    hhtem[3]->GetXaxis()->SetTitleSize(0.08);
    hhtem[3]->GetYaxis()->CenterTitle(true);           
    hhtem[3]->DrawClone("P");
 


    TLine *t1;
    if(_form==0) t1= new TLine(-7,1,11,1);
    if(_form==1) t1= new TLine(-2,1,11,1);
    t1->SetLineStyle(2);
	t1->Draw("same");  



    //______________________zoom in panel
    if(_form==0) pad[0][2]=new TPad("pad1","pad1",0.12,0.54,0.42,0.99);
    if(_form==1) pad[0][2]=new TPad("pad1","pad1",0.12,0.54,0.42,0.99);
    
    pad[0][2]->Draw();

    pad[0][2]->SetBottomMargin(0.35);
    pad[0][2]->SetLeftMargin(0.15);
    pad[0][2]->cd();

    gPad->SetFillStyle(4000); // make it transparent
    
    gPad->SetTickx(1);
    gPad->SetTicky(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);    

    char *xtile ="p_{T,jet}^{reco,ch} [Gev/c]";
    char *ytile1 ="SE/ME";

    hhtem[4] =  (TH1D*)    hRatio[_form] ->Clone();  

    hhtem[4]->SetMarkerStyle(29);
    hhtem[4]->SetMarkerSize(0.7);
    hhtem[4]->SetMarkerColor(2);
    hhtem[4]->SetLineColor(2);

    

    hhtem[4] -> GetXaxis()->SetTitle(xtile);
    hhtem[4] -> GetYaxis()->SetTitle(ytile1);
    hhtem[4] -> GetXaxis()->SetNdivisions(507);
    hhtem[4] -> GetYaxis()->SetNdivisions(507);
    if(_form==0)hhtem[4] -> GetXaxis()->SetRangeUser(regionC[0],regionC[1]-1);
    if(_form==1)hhtem[4] -> GetXaxis()->SetRangeUser(regionP[0],regionP[1]-1);

    hhtem[4] -> GetYaxis()->SetTitleOffset(0.6); 
    hhtem[4] -> GetXaxis()->SetTitleOffset(1.1);
    hhtem[4]->GetXaxis()->SetTitleSize(0.12);
    hhtem[4]->GetYaxis()->SetTitleSize(0.12);
    hhtem[4] -> GetYaxis()->SetLabelSize(0.11);
    hhtem[4] -> GetXaxis()->SetLabelSize(0.12);
    hhtem[4]->GetYaxis()->CenterTitle(true);  
    hhtem[4]->GetXaxis()->CenterTitle(true);          
    hhtem[4]->DrawClone("P");

    TLine *t2;
    if(_form==0) t2 = new TLine(regionC[0],1,regionC[1],1);
    if(_form==1) {
        t2 = new TLine(regionP[0],1,regionP[1],1);
        if(regionP[0]==0) t2 = new TLine(regionP[0],1,1,1);}
    t2->SetLineStyle(2);
	t2->Draw("same");
            
   
    
}

 void drawMEScale_detail(){
    
    readin();
    for(int i=0;i<2;i++){
        drawPtCompare(i);
        drawPt_beforeScaled_Compare(i);
    
    }
    
    
}
   
