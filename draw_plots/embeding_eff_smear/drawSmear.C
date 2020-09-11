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
const float Xhigh = 500;

const double PI = acos(-1.0);

TFile *fin;
TFile *fin_mix;


TFile *fout;

TLine* line;

TH1D *hEff_P[2];
TH1D *hEff_Pi[2];
TH1D *hEff_K[2];

TH2D *hTrkresponse_P[2];
TH2D *hTrkresponse_Pi[2];
TH2D *hTrkresponse_K[2];

//_________________________
TH2D *hResponse_P;
TH2D *hResponse_Pi;
TH2D *hResponse_K;
TH2D *htmp_response;

TH1D *htmp_sigma;

TH1D *hpT_dete[200];
//______________________________

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


Double_t Eff_track_rec_function(Double_t* x,Double_t* par)
{
    // Track reconstruction efficiency parametrization
    Double_t pt,y;
    Double_t A,B,C;A=par[0];B=par[1];C=par[2];
    pt=x[0];
    y=A*(exp(-pow(B/pt,C)));
    return y;
}

Double_t mom_res(Double_t* x,Double_t* par)
{
    Double_t pt,y;
    Double_t A,B,C;A=par[0];B=par[1];C=par[2];
    pt=x[0];
    y=A+B*pt+C*pt*pt;
    return y;
}

Double_t gausf(Double_t* x,Double_t* par){
    return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2/par[2]/par[2]);
}

void readin(int _from){
    
    //===============================  read in =====================================
    for(int icent=0;icent<2;icent++){
        sprintf(name,"embeding/MC_P_cent%i.root",icent);
        cout<<name<<endl;
        fin = TFile::Open(name);
    
        sprintf(name,"heffi");
        hEff_P[icent]     = (TH1D*)fin->Get(name);
    
        sprintf(name,"hpTRec_pTMc");
        hTrkresponse_P[icent]     = (TH2D*)fin->Get(name);
        //_____________________________________________________________
        sprintf(name,"embeding/MC_Pi_cent%i.root",icent);
        cout<<name<<endl;
        fin = TFile::Open(name);
    
        sprintf(name,"heffi");
        hEff_Pi[icent]     = (TH1D*)fin->Get(name);
    
        sprintf(name,"hpTRec_pTMc");
        hTrkresponse_Pi[icent]     = (TH2D*)fin->Get(name);

        //______________________________________________________
        sprintf(name,"embeding/MC_K_cent%i.root",icent);
        cout<<name<<endl;
        fin = TFile::Open(name);
    
        sprintf(name,"heffi");
        hEff_K[icent]     = (TH1D*)fin->Get(name);
    
        sprintf(name,"hpTRec_pTMc");
        hTrkresponse_K[icent]     = (TH2D*)fin->Get(name);
    
    
    }
    //===========================================================================================
    //_____________________________________________________________________________________________
    sprintf(name,"embeding/MC_P.root");
    cout<<name<<endl;
    fin = TFile::Open(name);

    sprintf(name,"hpTRec_pTMc");
    hResponse_P     = (TH2D*)fin->Get(name);
    //hResponse_P->Scale(1./2.);

    //_____________________________________________________________________________________________
    sprintf(name,"embeding/MC_Pi.root");
    cout<<name<<endl;
    fin = TFile::Open(name);

    sprintf(name,"hpTRec_pTMc");
    hResponse_Pi     = (TH2D*)fin->Get(name);
    //hResponse_Pi ->Scale(1./2.);

    //_____________________________________________________________________________________________
    sprintf(name,"embeding/MC_K.root");
    cout<<name<<endl;
    fin = TFile::Open(name);

    sprintf(name,"hpTRec_pTMc");
    hResponse_K     = (TH2D*)fin->Get(name);
    //hResponse_K->Scale(1./2.);

    //============================================
    

    int npTbins=200;
	double pTmax=20;
    for(int i=0;i<npTbins;i++){
        sprintf(name,"hpTMC_RecP_%d",i);
        hpT_dete[i] = new TH1D(name, name, npTbins, 0, pTmax);
    }  

    sprintf(name,"hsigma_mean_cent");
    htmp_sigma = new TH1D(name, name, npTbins, 0, pTmax);

    if(_from==0){
        htmp_response = (TH2D*) hResponse_P->Clone();
        
    }

    if(_from==1){
        htmp_response = (TH2D*) hResponse_Pi->Clone();
        
    }

    if(_from==2){
        htmp_response = (TH2D*) hResponse_K->Clone();
        
    } 
    

    TF1 *func = new TF1("func",gausf,0.2,20,3); 
    //_________________________________________________________________
    //cal mean and sigma
    
    int nptbins = htmp_response->FindLastBinAbove(0,2);
    int nptDbins = htmp_response->FindLastBinAbove(0,1);
    for(int j=1; j<nptbins+1; j++){
        for(int i=1; i<nptDbins+1; i++){
            double binc = htmp_response->GetBinContent(i,j);
            hpT_dete[j]->SetBinContent(i,binc);
        }
            
        func->SetParameters(1,hpT_dete[j]->GetMean(),hpT_dete[j]->GetRMS());
        func->SetParNames("Constant","Mean_value","Sigma");
        if(hpT_dete[j]->Integral()==0) continue;
        hpT_dete[j]->Fit("func");
        TF1 *para = hpT_dete[j]->GetFunction("func");
            
        double smear= para->GetParameter(2);
        htmp_sigma ->SetBinContent(j,smear);
    }
    

    
}


void drawPEffFit(int _from){
 
    int nx = 1;
    int ny = 1;
 
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
 
    sprintf(name, "efficiency_P_cent%d",_from);
 
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
 
            hhtem[0] =  (TH1D*)    hEff_P[_from]  ->Clone();
            
            
            hhtem[0] -> SetLineColor(clr[0]);
            
            
            
            char *xtile ="P_{T,true}";
            char *ytile ="#epsilon";
            char *ztile ="C(#Delta#eta #Delta#phi)";
 
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetRangeUser(0,5);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,2.5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
            TF1 *feff = new TF1("feff",Eff_track_rec_function,0.2,5,3);
            feff->SetParameters(0.7,0.1,2);
            feff->SetParNames("A","B","C");
            cout<<"fit for P eff: "<<endl;
           
 
            hhtem[0]->DrawClone("hist");
            hhtem[0]->Fit("feff","R");
            
            //hhtem[1]->DrawClone("same");
 
            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
 
            tx0=0.3;ty0=0.6;
            sprintf(name,"Proton");
            myTextF(tx0,ty0,name,tsize*0.8,1,12);
 
            tx0=0.65, ty0=0.6;
            if(_from==0) myTextF(tx0,ty0,"60-80%",tsize*0.8,1,12);
            if(_from==1) myTextF(tx0,ty0,"0-10%",tsize*0.8,1,12);

            
            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
 
            //delete can[0];
 
        }//ix
    }//iy
 

}

void drawPiEffFit(int _from){
 
    int nx = 1;
    int ny = 1;
 
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
 
    sprintf(name, "efficiency_Pi_cent%d",_from);
 
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
 
            hhtem[0] =  (TH1D*)    hEff_Pi[_from]  ->Clone();
            
            
            hhtem[0] -> SetLineColor(clr[0]);
            
            
            
            char *xtile ="P_{T,true}";
            char *ytile ="#epsilon";
            char *ztile ="C(#Delta#eta #Delta#phi)";
 
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetRangeUser(0,5);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,2.5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            TF1 *feff = new TF1("feff",Eff_track_rec_function,0.2,5,3);
            feff->SetParameters(0.7,0.1,2);
            feff->SetParNames("A","B","C");
            cout<<"fit for Pi eff: "<<endl;

        
            hhtem[0]->DrawClone("hist");
            hhtem[0]->Fit("feff","R");
            
            //hhtem[1]->DrawClone("same");
 
            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
 
            tx0=0.3;ty0=0.6;
            sprintf(name,"Pion");
            myTextF(tx0,ty0,name,tsize*0.8,1,12);
 
            tx0=0.65, ty0=0.6;
            if(_from==0) myTextF(tx0,ty0,"60-80%",tsize*0.8,1,12);
            if(_from==1) myTextF(tx0,ty0,"0-10%",tsize*0.8,1,12);

            
            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
 
            //delete can[0];
 
        }//ix
    }//iy
 
}

void drawKEffFit(int _from){
 
    int nx = 1;
    int ny = 1;
 
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
 
    sprintf(name, "efficiency_P_cent%d",_from);
 
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
 
            hhtem[0] =  (TH1D*)    hEff_K[_from]  ->Clone();
            
            
            hhtem[0] -> SetLineColor(clr[0]);
            
            
            
            char *xtile ="P_{T,true}";
            char *ytile ="#epsilon";
            char *ztile ="C(#Delta#eta #Delta#phi)";
 
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetRangeUser(0,5);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,2.5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            TF1 *feff = new TF1("feff",Eff_track_rec_function,0.2,5,3);
            feff->SetParameters(0.7,0.1,2);
            feff->SetParNames("A","B","C");
            cout<<"fit for K eff: "<<endl;
            

            hhtem[0]->DrawClone("hist");
            hhtem[0]->Fit("feff","R");
            
            //hhtem[1]->DrawClone("same");
 
            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
 
            tx0=0.3;ty0=0.6;
            sprintf(name,"Kaon");
            myTextF(tx0,ty0,name,tsize*0.8,1,12);
 
            
            tx0=0.65, ty0=0.6;
            if(_from==0) myTextF(tx0,ty0,"60-80%",tsize*0.8,1,12);
            if(_from==1) myTextF(tx0,ty0,"0-10%",tsize*0.8,1,12);

            
            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
 
            //delete can[0];
 
        }//ix
    }//iy
 
}


void drawSmearing(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "PTsmearing");
    
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
            
            hhtem[0] =  (TH1D*)    htmp_sigma  ->Clone();
            
            
            hhtem[0] -> SetLineColor(clr[0]);
            
            char *xtile ="P_{T,true}";
            char *ytile ="momentum resolution";
            char *ztile ="C(#Delta#eta #Delta#phi)";
            
            hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetRangeUser(0,5);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.5);
            
            
            TF1 *fsmear = new TF1("fsmear",mom_res,0.3,5,3);
            fsmear->SetParameters(0.04,-0.006,0.003);
            fsmear->SetParNames("A","B","C");
            
            
            hhtem[0]->DrawClone("hist");
            
            hhtem[0]->Fit("fsmear","R");
            //hhtem[1]->DrawClone("same");
            
            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
            
            tx0=0.3;ty0=0.6;
            sprintf(name,"Proton");
            myTextF(tx0,ty0,name,tsize*0.8,1,12);
            
            tx0=0.65, ty0=0.6;
            //if(_from==0) myTextF(tx0,ty0,"60-80%",tsize*0.8,1,12);
            //if(_from==1) myTextF(tx0,ty0,"0-10%",tsize*0.8,1,12);
            
            
            
            
            sprintf(name,"%s.pdf",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            //sprintf(name,"%s.gif",can[0]->GetName()); can[0]->SaveAs(name);
            
            //delete can[0];
            
        }//ix
    }//iy
    
}


void drawSmear(){
    
    //smear 0 p 1pi 2k
    readin(0);

    //1c 0p
    for(int ifrom = 0; ifrom != 1; ifrom++){
        
    
        //drawPEffFit(ifrom);
        //drawPiEffFit(ifrom);
        drawKEffFit(ifrom);
       
        
        
    }
    
     //drawSmearing();
    
}







