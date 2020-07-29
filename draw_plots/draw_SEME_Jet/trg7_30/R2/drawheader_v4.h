/*************************************************************************
  > File Name: drawheader_v3.h
 ************************************************************************/
#include <iostream>
#include <string>
#include <vector>
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TTree.h"
#include "TPad.h"
#include "TFile.h"
#include "TIterator.h"
#include "TObject.h"
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGaxis.h>
#include <fstream>
#include <TPaveText.h>
#include <TVirtualPad.h>
#include <TPaveLabel.h>
#include <TGraph.h>
#include <TFrame.h>
#include <TRandom3.h>
#include <TAttFill.h>
#include "TAttMarker.h"
#include "TMarker.h"
#include "TPave.h"
#include <algorithm>    // std::min_element, std::max_element

using namespace std;
TGraphErrors* myGrE_GrE(TGraphErrors* g1,TGraphErrors* g2) {

    const Int_t debug=0;

    if (!g1) printf("**myTGraphErrorsDivide: g1 does not exist !  \n");
    if (!g2) printf("**myTGraphErrorsDivide: g2 does not exist !  \n");


    Int_t n1=g1->GetN();
    Int_t n2=g2->GetN();

    if (n1!=n2) {
        printf("**myTGraphErrorsDivide: vector do not have same number of entries !  \n");
    }

    TGraphErrors* g3= new TGraphErrors();

    Double_t  x1=0., y1=0., x2=0., y2=0.;
    Double_t dx1=0.,dy1=0.,       dy2=0.;

    Int_t iv=0;
    for (Int_t i1=0; i1<n1; i1++) {
        for (Int_t i2=0; i2<n2; i2++) {

            g1->GetPoint(i1,x1,y1);
            g2->GetPoint(i2,x2,y2);
            //            cout<<x1<<" x "<<x2<<endl;
            if (x1!=x2) {
            }else{

                dx1  = g1->GetErrorX(i1);
                dy1  = g1->GetErrorY(i1);
                dy2  = g2->GetErrorY(i2);

                cout<<x1<<" "<<x2<<" "<<y1<<" "<<y2<<endl;
                g3->SetPoint(iv, x1,y1 - y2);

                Double_t e=0.;
                e=std::sqrt(dy1*dy1+dy2*dy2);
                g3->SetPointError(iv,dx1,e);

                iv++;
            }
        }
    }
    return g3;

}

TH1* pr2hist(TProfile* pr){
    int Nbx = pr->GetNbinsX();
    double* arrayx = new double[Nbx+1];

    for(int ib=0; ib<Nbx; ib++){
        arrayx[ib] = pr->GetBinCenter(ib+1) - 0.5* pr->GetBinWidth(ib+1);
    }
    arrayx[Nbx] = pr->GetBinCenter(Nbx)+0.5* pr->GetBinWidth(Nbx);

    char name[200];
    sprintf(name,"%s_h", pr->GetName());
    TH1* hist = new TH1D(name,"",Nbx,arrayx);
    for(int ib=0; ib<Nbx; ib++){
        hist ->SetBinContent( ib+1, pr->GetBinContent(ib+1));
        hist ->SetBinError( ib+1, pr->GetBinError(ib+1));
    }
    delete [] arrayx;
    return hist;
}


Int_t Tclr[] = {kMagenta+1,kRed+1,kBlue, 1,  kGreen+1,kGray+2, kCyan+1, kYellow+1,kViolet+1};

/*

   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetErrorX(0);

   TLatex text;
   text.SetNDC(1);text.SetTextFont(43);text.SetTextSize(25);
   text.DrawLatex(0.26,0.06,etaname[fb][ieta]);

   double xx = 0.17;
   double yy = 0.91;
   double itv = 0.05;
   int fsize = 23;
   SetTextAlign(11); //1,2,3
   text.DrawLatex(xx,yy,"ATLAS");
   text.SetTextFont(43);text.SetTextSize(fsize); text.DrawLatex(xx+0.16,yy,"internal"); 
   text.SetTextFont(43);text.SetTextSize(fsize); text.DrawLatex(xx,yy-itv,"287281"); 
   sprintf(name,"Centrality %d-%d%%", ic, ic*5+5);
   text.SetTextFont(43);text.SetTextSize(fsize); text.DrawLatex(xx,yy-2*itv,name); 

   TCanvas* can;

   TF1 *myfunc = h->GetFunction("myfunc");

   TPad* pad[10][10];
   sprintf(name,"%s_pad_%i_%i",can->GetName(),ix, iy);
   pad[ix][iy] = (TPad*) gROOT->FindObject(name);
   pad[ix][iy]->cd();
   */


void setstyle1F(TH1*h, float tsize=0.07, double Xtofset=1.1,  double Ytofset=1.1, float lsize=0.06, double lbfset=0.01, int ndx=505, int ndy=505){
    //text font are given in pixels pad is 400*400
    //Rectangular and
    int font=42;
    // h->GetYaxis()->CenterTitle();            h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetNdivisions(ndx);       h->GetXaxis()->SetNdivisions(ndy);
    h->GetYaxis()->SetTitleOffset(Ytofset);  h->GetXaxis()->SetTitleOffset(Xtofset);
    h->GetYaxis()->SetLabelOffset(lbfset);   h->GetXaxis()->SetLabelOffset(lbfset);
    h->GetYaxis()->SetTitleFont(font);       h->GetXaxis()->SetTitleFont(font);
    h->GetYaxis()->SetTitleSize(tsize);      h->GetXaxis()->SetTitleSize(tsize);
    h->GetYaxis()->SetLabelFont(font);       h->GetXaxis()->SetLabelFont(font);
    h->GetYaxis()->SetLabelSize(lsize);      h->GetXaxis()->SetLabelSize(lsize);
    h->GetXaxis()->SetTickLength(0.02);      h->GetYaxis()->SetTickLength(0.02);
}

void setstyle2F(TH2*h, float tsize=0.07, double Xtofset=1.1, double Ytofset=1.1, float lsize=0.06, double lbfset=0.01, int ndx=505, int ndy=505){
    //text font are given in pixels pad is 400*400
    //Rectangular and
    int font=42;
    h->GetYaxis()->CenterTitle();            h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetNdivisions(ndx);       h->GetXaxis()->SetNdivisions(ndy);
    h->GetYaxis()->SetTitleOffset(Ytofset);  h->GetXaxis()->SetTitleOffset(Xtofset);
    h->GetYaxis()->SetLabelOffset(lbfset);   h->GetXaxis()->SetLabelOffset(lbfset);
    h->GetYaxis()->SetTitleFont(font);       h->GetXaxis()->SetTitleFont(font);
    h->GetYaxis()->SetTitleSize(tsize);      h->GetXaxis()->SetTitleSize(tsize);
    h->GetYaxis()->SetLabelFont(font);       h->GetXaxis()->SetLabelFont(font);
    h->GetYaxis()->SetLabelSize(lsize);      h->GetXaxis()->SetLabelSize(lsize);
    h->GetXaxis()->SetTickLength(0.02);      h->GetYaxis()->SetTickLength(0.02);
}

void setstyle1(TH1*h, int tsize=34, double Xtofset=2.2,  double Ytofset=2.8, int lsize=30, double lbfset=0.01, int ndx=505, int ndy=505){
    //text font are given in pixels pad is 400*400
    //Rectangular and 
    int font=43;
    //h->GetYaxis()->CenterTitle();            h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetNdivisions(ndx);       h->GetXaxis()->SetNdivisions(ndy);
    h->GetYaxis()->SetTitleOffset(Ytofset);  h->GetXaxis()->SetTitleOffset(Xtofset);
    h->GetYaxis()->SetLabelOffset(lbfset);   h->GetXaxis()->SetLabelOffset(lbfset);
    h->GetYaxis()->SetTitleFont(font);       h->GetXaxis()->SetTitleFont(font);
    h->GetYaxis()->SetTitleSize(tsize);      h->GetXaxis()->SetTitleSize(tsize);
    h->GetYaxis()->SetLabelFont(font);       h->GetXaxis()->SetLabelFont(font);
    h->GetYaxis()->SetLabelSize(lsize);      h->GetXaxis()->SetLabelSize(lsize);
    h->GetXaxis()->SetTickLength(0.02);      h->GetYaxis()->SetTickLength(0.02);
}

void setstyle2(TH2*h, int tsize=34, double Xtofset=2.2,  double Ytofset=2.8, int lsize=30, double lbfset=0.01, int ndx=505, int ndy=505){
    //text font are given in pixels pad is 400*400
    //Rectangular and
    int font=43;
    //h->GetYaxis()->CenterTitle();            h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetNdivisions(ndx);       h->GetXaxis()->SetNdivisions(ndy);
    h->GetYaxis()->SetTitleOffset(Ytofset);  h->GetXaxis()->SetTitleOffset(Xtofset);
    h->GetYaxis()->SetLabelOffset(lbfset);   h->GetXaxis()->SetLabelOffset(lbfset);
    h->GetYaxis()->SetTitleFont(font);       h->GetXaxis()->SetTitleFont(font);
    h->GetYaxis()->SetTitleSize(tsize);      h->GetXaxis()->SetTitleSize(tsize);
    h->GetYaxis()->SetLabelFont(font);       h->GetXaxis()->SetLabelFont(font);
    h->GetYaxis()->SetLabelSize(lsize);      h->GetXaxis()->SetLabelSize(lsize);
    h->GetXaxis()->SetTickLength(0.02);      h->GetYaxis()->SetTickLength(0.02);
}

void setGrErr(TGraphErrors*h, int tsize=34, double Xtofset=2.2,  double Ytofset=2.8, int lsize=30, double lbfset=0.01, int ndx=505, int ndy=505){
    //text font are given in pixels pad is 400*400
    //Rectangular and
    int font=43;
    //h->GetYaxis()->CenterTitle();            h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetNdivisions(ndx);       h->GetXaxis()->SetNdivisions(ndy);
    h->GetYaxis()->SetTitleOffset(Ytofset);  h->GetXaxis()->SetTitleOffset(Xtofset);
    h->GetYaxis()->SetLabelOffset(lbfset);   h->GetXaxis()->SetLabelOffset(lbfset);
    h->GetYaxis()->SetTitleFont(font);       h->GetXaxis()->SetTitleFont(font);
    h->GetYaxis()->SetTitleSize(tsize);      h->GetXaxis()->SetTitleSize(tsize);
    h->GetYaxis()->SetLabelFont(font);       h->GetXaxis()->SetLabelFont(font);
    h->GetYaxis()->SetLabelSize(lsize);      h->GetXaxis()->SetLabelSize(lsize);
    h->GetXaxis()->SetTickLength(0.02);      h->GetYaxis()->SetTickLength(0.02);
}

void setGrErrF(TGraphErrors*h, float tsize=0.07, double Xtofset=1.1,  double Ytofset=1.1, float lsize=0.07, double lbfset=0.01, int ndx=505, int ndy=505){
    //text font are given in pixels pad is 400*400
    //Rectangular and
    int font=42;
    // h->GetYaxis()->CenterTitle();            h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetNdivisions(ndx);       h->GetXaxis()->SetNdivisions(ndy);
    h->GetYaxis()->SetTitleOffset(Ytofset);  h->GetXaxis()->SetTitleOffset(Xtofset);
    h->GetYaxis()->SetLabelOffset(lbfset);   h->GetXaxis()->SetLabelOffset(lbfset);
    h->GetYaxis()->SetTitleFont(font);       h->GetXaxis()->SetTitleFont(font);
    h->GetYaxis()->SetTitleSize(tsize);      h->GetXaxis()->SetTitleSize(tsize);
    h->GetYaxis()->SetLabelFont(font);       h->GetXaxis()->SetLabelFont(font);
    h->GetYaxis()->SetLabelSize(lsize);      h->GetXaxis()->SetLabelSize(lsize);
    h->GetXaxis()->SetTickLength(0.02);      h->GetYaxis()->SetTickLength(0.02);
}

void setGr(TGraph*h, int tsize=34, double Xtofset=2.2,  double Ytofset=2.8, int lsize=30, double lbfset=0.01, int ndx=505, int ndy=505){
    //text font are given in pixels pad is 400*400
    //Rectangular and
    int font=43;
    //h->GetYaxis()->CenterTitle();            h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetNdivisions(ndx);       h->GetXaxis()->SetNdivisions(ndy);
    h->GetYaxis()->SetTitleOffset(Ytofset);  h->GetXaxis()->SetTitleOffset(Xtofset);
    h->GetYaxis()->SetLabelOffset(lbfset);   h->GetXaxis()->SetLabelOffset(lbfset);
    h->GetYaxis()->SetTitleFont(font);       h->GetXaxis()->SetTitleFont(font);
    h->GetYaxis()->SetTitleSize(tsize);      h->GetXaxis()->SetTitleSize(tsize);
    h->GetYaxis()->SetLabelFont(font);       h->GetXaxis()->SetLabelFont(font);
    h->GetYaxis()->SetLabelSize(lsize);      h->GetXaxis()->SetLabelSize(lsize);
    h->GetXaxis()->SetTickLength(0.02);      h->GetYaxis()->SetTickLength(0.02);
}

void setGrF(TGraph*h, float tsize=0.07, double Xtofset=1.1,  double Ytofset=1.1, float lsize=0.06, double lbfset=0.01, int ndx=505, int ndy=505){
    //text font are given in pixels pad is 400*400
    //Rectangular and
    int font=42;
    // h->GetYaxis()->CenterTitle();            h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetNdivisions(ndx);       h->GetXaxis()->SetNdivisions(ndy);
    h->GetYaxis()->SetTitleOffset(Ytofset);  h->GetXaxis()->SetTitleOffset(Xtofset);
    h->GetYaxis()->SetLabelOffset(lbfset);   h->GetXaxis()->SetLabelOffset(lbfset);
    h->GetYaxis()->SetTitleFont(font);       h->GetXaxis()->SetTitleFont(font);
    h->GetYaxis()->SetTitleSize(tsize);      h->GetXaxis()->SetTitleSize(tsize);
    h->GetYaxis()->SetLabelFont(font);       h->GetXaxis()->SetLabelFont(font);
    h->GetYaxis()->SetLabelSize(lsize);      h->GetXaxis()->SetLabelSize(lsize);
    h->GetXaxis()->SetTickLength(0.02);      h->GetYaxis()->SetTickLength(0.02);
}


TCanvas* newDivCan2(char* name, double* Lmrg, double* llmrg, double* ratx, double* raty,  int Nx=1, int Ny=1, int x_length=400, int  y_width = 400){
    double Lmargin = Lmrg[0]; double Bmargin =  Lmrg[1];
    double lmarg = llmrg[2];  double rmarg = llmrg[3]; double tmarg = llmrg[0]; double bmarg = llmrg[1];

    //ratio plots
    //  double   ratx[]={1,1,1,1,1,1,1,1,1,1,1,1,1};
    //  double   raty[]={1,0.35,0.35,1,0.35,0.35,1,1,1,1,1,1,1,1,1};
    //  double   raty[]={1,0.8,0.8,1,0.8,0.8,1,1,1,1,1,1,1,1,1};
    //   ; // 1st pad
    //   ; // Ny pad
    //   ; //begin from 2-->Nx pad
    //for all pad
    //from 1-->Ny-1 pad
    //from 1-->Ny-1 pad

    //if same pad, then  Lmargin = lmarg
    //if only left bar, diable

    double x0[100] = {0};
    double x1[100] = {0};
    double y0[100] = {0};
    double y1[100] = {0};

    double hspace = 0;
    double vspace = 0;

    for(int ix=0; ix<Nx;  ix++){
        if(ix==0){
            x0[ix] = 0;
            x1[ix] = Lmargin+ ratx[ix] + rmarg;
        }else{
            x0[ix] = x1[ix-1]+hspace;
            x1[ix] = x0[ix] + lmarg + ratx[ix] + rmarg;
        }
    }

    for(int iy=0; iy<Ny; iy++){
        if( iy==0){
            y0[Ny-iy-1] = 0;
            y1[Ny-iy-1] = y0[Ny-iy-1] + Bmargin + raty[Ny-iy-1] +tmarg ;
        }else{

            y0[Ny-iy-1] = y1[Ny-iy] + vspace;
            y1[Ny-iy-1] = y0[Ny-iy-1] + bmarg+raty[Ny-iy-1] + tmarg;
        }
    }

    double sumX = x1[Nx-1];
    double sumY = y1[0];

    int lenght  =sumX *x_length;
    int withd  = sumY*y_width;


    double hx[100] = {0};
    double hy[100] = {0};

    for(int ix = 0; ix<Nx; ix++){ hx[ix] = x1[ix] - x0[ix];  x0[ix] /= sumX; x1[ix] /= sumX;
    }
    for(int iy = 0; iy<Ny; iy++){ hy[iy] = y1[iy] - y0[iy];
        y0[iy] /= sumY; y1[iy] /= sumY;
    }

    TCanvas*   can = new TCanvas(name,"", lenght+4, withd+28);
    cout<<"The width "<<lenght<<" "<<withd<<endl;

    for(int ix=0; ix<Nx; ix++){
        for(int iy=0; iy<Ny; iy++){

            can->cd( ix+iy*Nx+1);
            double right = rmarg;
            double left = lmarg;
            double top = tmarg;
            double bottom  = bmarg;

            if(ix==0) left = Lmargin;
            if(iy==Ny-1) bottom = Bmargin;
            left /= hx[ix]; right /= hx[ix];
            top /= hy[iy]; bottom /= hy[iy];

            char namex[200];
            sprintf(namex,"pad_%i_%i",ix,iy);
            TPad *pad = (TPad*) gROOT->FindObject(namex);
            if(pad) delete pad;
            sprintf(namex,"%s_pad_%i_%i",can->GetName(),ix, iy);
            pad = new TPad(namex,"",x0[ix],y0[iy],x1[ix],y1[iy]);
            //pad->SetPad( x0[ix],y0[iy],x1[ix],y1[iy]);

            pad->SetLeftMargin(left);
            pad->SetRightMargin(right);
            pad->SetTopMargin(top);
            pad->SetBottomMargin(bottom);

            //pad->SetFillStyle(3001);
            //pad->SetFillColor(ix+iy+1);
            pad->Draw();
        }
    }

    can->SetWindowSize(lenght + (lenght - can->GetWw()), withd + (withd - can->GetWh()));
    cout<<"The width "<<can->GetWw()<<" "<<can->GetWh()<<endl;
    return can;

}

TCanvas*  SetCan2D2( char*name,double* llmrg, int nx=1, int ny=1, int x_length=400, int  y_width = 400 ){

    double lmarg = llmrg[2];  double rmarg = llmrg[3]; double tmarg = llmrg[0]; double bmarg = llmrg[1];

    int withd  = x_length/(1-lmarg-rmarg);
    int lenght = y_width/(1-tmarg-bmarg);

    TCanvas*    can = new TCanvas(name,"", nx*withd, ny*lenght);

    if( nx==1 && ny==1){
        cout<<"We set margin"<<endl;
        gPad->SetLeftMargin(lmarg);
        gPad->SetRightMargin(rmarg);
        gPad->SetTopMargin(tmarg);
        gPad->SetBottomMargin(bmarg);
    }else{

        can->Divide( nx, ny,0,0);
        for(int i=0; i<nx*ny; i++){
            can->cd(i+1);

            gPad->SetLeftMargin(lmarg);
            gPad->SetRightMargin(rmarg);
            gPad->SetTopMargin(tmarg);
            gPad->SetBottomMargin(bmarg);

        }
    }


    return can;
}

TCanvas*  SetCan1Rec( char *name){

    TCanvas*    can = new TCanvas(name,"", 800, 600);
    std::cout<<"We make: "<<can->GetName()<<std::endl;
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.16);
    gPad->SetLeftMargin(0.16);

    return can;
}

TCanvas*  SetCan1Box( char *name){

    TCanvas*    can = new TCanvas(name,"", 600, 600);
    std::cout<<"We make: "<<can->GetName()<<std::endl;
    gPad->SetTopMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.16);
    gPad->SetLeftMargin(0.16);


    return can;
}

TLegend* mylegF(float x0=0.6, float y0=0.75, float x1=0.95, float y1 = 0.95, float fsize=0.06 ){
    TLegend*  leg = new TLegend(x0, y0, x1, y1);
    leg->SetTextFont(42);
    leg->SetTextSize(fsize);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    return leg;
}

TLegend* myleg(float x0=0.6, float y0=0.75, float x1=0.95, float y1 = 0.95, int fsize=23 ){
    TLegend*  leg = new TLegend(x0, y0, x1, y1);
    leg->SetTextFont(43);
    leg->SetTextSize(fsize);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    return leg;
}

void GetGraAsyRg(TGraphAsymmErrors* gr[100], int n, float& min, float& max, float& distance){

    max=0;
    distance = 1;
    min=1000000000;

    for(int i=0; i<n; i++){
        int Npoit = gr[i]->GetN();
        double* yval = gr[i]->GetY();
        double* yvallowE = gr[i]->GetEYlow();
        double* yvalhighE = gr[i]->GetEYhigh();

        double yvalUpp[100] = {0};
        double yvalLow[100] = {0};

        for (int k=0; k<Npoit; k++) {
            yvalUpp[k] = yval[k] + yvalhighE[k];
            yvalLow[k] = yval[k] - yvallowE[k];
        }

        double max1 =  *std::max_element(yvalUpp,yvalUpp+Npoit);
        double min1 =  *std::min_element(yvalLow,yvalLow+Npoit);
        if(max1>=max) max = max1;
        if(min>=min1) min = min1;
    }

    distance  = max-min;
    if(min<=0) cout<<"Catious when SetLogy "<<min<<"  "<<endl;
}

void GetGraRg(TGraphErrors* gr[100], int n, float& min, float& max, float& distance){

    max=0;
    distance = 1;
    min=1000000000;

    for(int i=0; i<n; i++){
        int Npoit = gr[i]->GetN();
        double* yval = gr[i]->GetY();
        double* yvalE = gr[i]->GetEY();

        double yvalUpp[100] = {0};
        double yvalLow[100] = {0};

        for (int k=0; k<Npoit; k++) {
            yvalUpp[k] = yval[k] + yvalE[k];
            yvalLow[k] = yval[k] - yvalE[k];
        }

        double max1 =  *std::max_element(yvalUpp,yvalUpp+Npoit);
        double min1 =  *std::min_element(yvalLow,yvalLow+Npoit);
        if(max1>=max) max = max1;
        if(min>=min1) min = min1;
    }

    distance  = max-min;
    if(min<=0) cout<<"Catious when SetLogy "<<min<<"  "<<endl;
}

void GetGrRg(TGraph* gr[100], int n, float& min, float& max, float& distance){

    max=0;
    distance = 1;
    min=1000000000;

    for(int i=0; i<n; i++){
        int Npoit = gr[i]->GetN();
        double* yval = gr[i]->GetY();

        double yvalUpp[100] = {0};
        double yvalLow[100] = {0};

        for (int k=0; k<Npoit; k++) {
            yvalUpp[k] = yval[k] ;
            yvalLow[k] = yval[k] ;
        }

        double max1 =  *std::max_element(yvalUpp,yvalUpp+Npoit);
        double min1 =  *std::min_element(yvalLow,yvalLow+Npoit);
        if(max1>=max) max = max1;
        if(min>=min1) min = min1;
    }

    distance  = max-min;
    if(min<=0) cout<<"Catious when SetLogy "<<min<<"  "<<endl;
}

void GetHistRg(TH1* httmp1[100], int n, float& min, float& max, float& distance){

    max=0;
    distance = 1;
    min=1000000000;


    for(int i=0; i<n; i++){
        double max1 = httmp1[i]->GetBinContent( httmp1[i] ->GetMaximumBin());
        double min1 = httmp1[i]->GetBinContent( httmp1[i] ->GetMinimumBin());
        if(i==0){
            max = max1;
            min = min1;
        }
        else{
            if(max1>=max) max = max1;
            if(min>=min1) min = min1;
        }

    }
    distance  = max-min;

    if(min<=0) cout<<"Catious when SetLogy "<<min<<"  "<<endl;

}

void mydrawHist( TH1* htmp1[100], int Nfile, int* mrk, int* clr, char* xtile, char* ytile, double rgd = 0.4, double rgu=0.4){
    gPad->cd();
    gPad-> SetTicks(1,1);

    float max=1; float min=0; float distance =2;
    GetHistRg( htmp1, Nfile,  min,  max, distance);
    min -= rgd*distance;
    max += rgu*distance;

    for (int k=0; k<Nfile; k++) {
        htmp1[k] ->GetXaxis()->SetTitle(xtile);
        htmp1[k] ->GetYaxis()->SetTitle(ytile);

        htmp1[k]->SetMarkerColor(clr[k]);
        htmp1[k]->SetLineColor(clr[k]);
        htmp1[k]->SetMarkerStyle(mrk[k]);

        htmp1[k]->SetMarkerSize(1.3);
        htmp1[k]->GetYaxis()->SetRangeUser(min, max);

        //setstyle11( htmp1[k], 34,4.0,2.8,30,0.01);
        if(k==0)  htmp1[k]->DrawCopy("p");
        else htmp1[k]->DrawCopy("psame");
    }
}

void mydrawPr( TProfile* hpr[100], int Nfile, int* mrk, int* clr, char* xtile, char* ytile,  double rgd = 0.4, double rgu=0.4){
    gPad->cd();
    gPad-> SetTicks(1,1);

    TH1* htmp1x[100];

    for(int k=0; k<Nfile; k++) htmp1x[k] = pr2hist( hpr[k]);

    float max=1; float min=-1; float distance =0;
    GetHistRg( htmp1x, Nfile,  min,  max, distance);
    min -= rgd*distance;
    max += rgu*distance;

    for (int k=0; k<Nfile; k++) {
        htmp1x[k] ->GetXaxis()->SetTitle(xtile);
        htmp1x[k] ->GetYaxis()->SetTitle(ytile);

        htmp1x[k]->SetMarkerColor(clr[k]);
        htmp1x[k]->SetLineColor(clr[k]);
        htmp1x[k]->SetMarkerStyle(mrk[k]);

        htmp1x[k]->SetMarkerSize(1.3);
        htmp1x[k]->GetYaxis()->SetRangeUser(min, max);

        if(k==0)  htmp1x[k]->DrawCopy("p");
        else htmp1x[k]->DrawCopy("psame");

    }

    for(int k=0; k<Nfile; k++) delete htmp1x[k];
}

void mydrawGr( TH1* h1, TGraph* hgrr[100], int Nfile, int* mrk, int* clr, char* xtile, char* ytile, bool drawh1=true, double rgd = 0.2, double rgu = 0.2, float mrksize=1.3, char* ctr="pl"){

    gPad->cd();
    gPad-> SetTicks(1,1);

    float max=1; float min=-1; float distance =0;
    GetGrRg( hgrr, Nfile,  min,  max, distance);
    min -= rgd*distance;
    max += rgu*distance;

    if( drawh1) {
        h1->GetYaxis()->SetRangeUser(min, max);
        h1->SetXTitle(xtile);
        h1->SetYTitle(ytile);
        h1->DrawCopy();
    }

    for (int k=0; k<Nfile; k++) {
        hgrr[k] ->GetXaxis()->SetTitle(xtile);
        hgrr[k] ->GetYaxis()->SetTitle(ytile);

        hgrr[k]->SetMarkerColor(clr[k]);
        hgrr[k]->SetLineColor(clr[k]);
        hgrr[k]->SetMarkerStyle(mrk[k]);

        hgrr[k]->SetMarkerSize(mrksize);

        if(k==0)  hgrr[k]->DrawClone(ctr);
        else hgrr[k]->DrawClone(ctr);

    }

}
void mydrawGrErr( TH1* h1, TGraphErrors* hgrr[100], int Nfile, int* mrk, int* clr, char* xtile, char* ytile, bool drawh1=true,double rgd = 0.2, double rgu = 0.2, float mrksize=1.3, char* ctr="pl"){

    gPad->cd();
    gPad-> SetTicks(1,1);

    float max=1; float min=-1; float distance =0;
    GetGraRg( hgrr, Nfile,  min,  max, distance);
    min -= rgd*distance;
    max += rgu*distance;

    if( drawh1) {
        h1->GetYaxis()->SetRangeUser(min, max);
        h1->SetXTitle(xtile);
        h1->SetYTitle(ytile);
        h1->DrawCopy();
    }

    for (int k=0; k<Nfile; k++) {
        hgrr[k] ->GetXaxis()->SetTitle(xtile);
        hgrr[k] ->GetYaxis()->SetTitle(ytile);

        hgrr[k]->SetMarkerColor(clr[k]);
        hgrr[k]->SetLineColor(clr[k]);
        hgrr[k]->SetMarkerStyle(mrk[k]);

        hgrr[k]->SetMarkerSize(mrksize);

        char name[200];

        sprintf(name, "a%s", ctr);

        if(k==0) {if (drawh1) hgrr[k]->DrawClone(ctr);
            else  hgrr[k]->DrawClone(name); }
        else hgrr[k]->DrawClone(ctr);

    }

}
void mydrawGrAsyErr( TH1* h1, TGraphAsymmErrors* hgrr[100], int Nfile, int* mrk, int* clr, char* xtile, char* ytile,bool drawh1=true,double rgd = 0.2, double rgu = 0.2,float mrksize=1.3){
    //use h1 as template and draw the tgaph together
    gPad->cd();
    gPad-> SetTicks(1,1);

    float max=1; float min=-1; float distance =0;
    GetGraAsyRg( hgrr, Nfile,  min,  max, distance);
    min -= rgd*distance;
    max += rgu*distance;

    h1->GetYaxis()->SetRangeUser(min, max);

    if( drawh1) {
        h1->GetYaxis()->SetRangeUser(min, max);
        h1->SetXTitle(xtile);
        h1->SetYTitle(ytile);
        h1->DrawCopy();
    }

    for (int k=0; k<Nfile; k++) {
        hgrr[k] ->GetXaxis()->SetTitle(xtile);
        hgrr[k] ->GetYaxis()->SetTitle(ytile);

        hgrr[k]->SetMarkerColor(clr[k]);
        hgrr[k]->SetLineColor(clr[k]);
        hgrr[k]->SetMarkerStyle(mrk[k]);

        hgrr[k]->SetMarkerSize(mrksize);

        if(k==0)  hgrr[k]->DrawClone("p");
        else hgrr[k]->DrawClone("p");
    }

}


TGraphErrors* myTGraphErrorsDivide(TGraphErrors* g1,TGraphErrors* g2) {

    const Int_t debug=0;

    if (!g1) printf("**myTGraphErrorsDivide: g1 does not exist !  \n");
    if (!g2) printf("**myTGraphErrorsDivide: g2 does not exist !  \n");


    Int_t n1=g1->GetN();
    Int_t n2=g2->GetN();

    if (n1!=n2) {
        printf("**myTGraphErrorsDivide: vector do not have same number of entries !  \n");
    }

    TGraphErrors* g3= new TGraphErrors();

    Double_t  x1=0., y1=0., x2=0., y2=0.;
    Double_t dx1=0.,dy1=0.,       dy2=0.;

    Int_t iv=0;
    for (Int_t i1=0; i1<n1; i1++) {
        for (Int_t i2=0; i2<n2; i2++) {
            //if (debug) printf("**myTGraphErrorsDivide: %d  %d !  \n",i1,i2);

            g1->GetPoint(i1,x1,y1);
            g2->GetPoint(i2,x2,y2);
            if (x1!=x2) {
                //printf("**myTGraphErrorsDivide: %d x1!=x2  %f %f  !  \n",iv,x1,x2);
            }else{
                //if (debug) printf("**myTGraphErrorsDivide: %d x1=x2  %f %f  !  \n",iv,x1,x2);
                dx1  = g1->GetErrorX(i1);
                if (y1!=0) dy1  = g1->GetErrorY(i1)/y1;
                if (y2!=0) dy2  = g2->GetErrorY(i2)/y2;

                if (debug)
                    printf("**myTGraphErrorsDivide: %d x1=%f x2=%f y1=%f y2=%f  \n",iv,x1,x2,y1,y2);

                if (y2!=0.) g3->SetPoint(iv, x1,y1/y2);
                else        g3->SetPoint(iv, x1,y2);

                Double_t e=0.;
                if (y1!=0 && y2!=0) e=std::sqrt(dy1*dy1+dy2*dy2)*(y1/y2);
                g3->SetPointError(iv,dx1,e);


                if (debug) {
                    //Double_t g3y, g3x,g3e;
                    //g3->GetPoint(iv, g3y,g3x);
                    //g3e=g3->GetErrorY(iv);
                    //printf("%d g3y= %f g3e=%f  \n",iv,g3y,g3e);
                }
                iv++;
            }
            //    printf("**myTGraphErrorsDivide: ...next  \n");
        }
    }
    return g3;

}

TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2) {

    const Int_t debug=0;

    TGraphAsymmErrors* g3= new TGraphAsymmErrors();
    Int_t n1=g1->GetN();
    Int_t n2=g2->GetN();

    if (n1!=n2) {
        printf(" vectors do not have same number of entries !  \n");
        return g3;
    }

    Double_t   x1=0.,   y1=0., x2=0., y2=0.;
    Double_t dx1h=0., dx1l=0.;
    Double_t dy1h=0., dy1l=0.;
    Double_t dy2h=0., dy2l=0.;

    Double_t* X1 = g1->GetX();
    Double_t* Y1 = g1->GetY();
    Double_t* EXhigh1 = g1->GetEXhigh();
    Double_t* EXlow1 =  g1->GetEXlow();
    Double_t* EYhigh1 = g1->GetEYhigh();
    Double_t* EYlow1 =  g1->GetEYlow();

    Double_t* X2 = g2->GetX();
    Double_t* Y2 = g2->GetY();
    Double_t* EXhigh2 = g2->GetEXhigh();
    Double_t* EXlow2 =  g2->GetEXlow();
    Double_t* EYhigh2 = g2->GetEYhigh();
    Double_t* EYlow2 =  g2->GetEYlow();

    for (Int_t i=0; i<g1->GetN(); i++) {
        g1->GetPoint(i,x1,y1);
        g2->GetPoint(i,x2,y2);
        dx1h  = EXhigh1[i];
        dx1l  = EXlow1[i];
        if (y1!=0.) dy1h  = EYhigh1[i]/y1;
        else        dy1h  = 0.;
        if (y2!=0.) dy2h  = EYhigh2[i]/y2;
        else        dy2h  = 0.;
        if (y1!=0.) dy1l  = EYlow1 [i]/y1;
        else        dy1l  = 0.;
        if (y2!=0.) dy2l  = EYlow2 [i]/y2;
        else        dy2l  = 0.;

        //if (debug)
        //printf("%d x1=%f x2=%f y1=%f y2=%f  \n",i,x1,x2,y1,y2);
        if (debug)
            printf("%d dy1=%f %f dy2=%f %f sqrt= %f %f \n",i,dy1l,dy1h,dy2l,dy2h,
                    std::sqrt(dy1l*dy1l+dy2l*dy2l), std::sqrt(dy1h*dy1h+dy2h*dy2h));

        if (y2!=0.) g3->SetPoint(i, x1,y1/y2);
        else       g3->SetPoint(i, x1,y2);
        Double_t el=0.; Double_t eh=0.;

        if (y1!=0. && y2!=0.) el=std::sqrt(dy1l*dy1l+dy2l*dy2l)*(y1/y2);
        if (y1!=0. && y2!=0.) eh=std::sqrt(dy1h*dy1h+dy2h*dy2h)*(y1/y2);

        if (debug) printf("dx1h=%f  dx1l=%f  el=%f  eh=%f \n",dx1h,dx1l,el,eh);
        g3->SetPointError(i,dx1h,dx1l,el,eh);

    }
    return g3;

}

TGraphAsymmErrors* myMakeBand(TGraphErrors* g0, TGraphErrors* g1,TGraphErrors* g2) {
    // default is g0
    //const Int_t debug=0;

    TGraphAsymmErrors* g3= new TGraphAsymmErrors();

    Double_t  x1=0., y1=0., x2=0., y2=0., y0=0, x3=0.;
    //Double_t dx1=0.;
    Double_t dum;
    for (Int_t i=0; i<g1->GetN(); i++) {
        g0->GetPoint(i, x1,y0);
        g1->GetPoint(i, x1,y1);
        g2->GetPoint(i, x1,y2);

        // if (y1==0) y1=1;
        //if (y2==0) y2=1;

        if (i==g1->GetN()-1) x2=x1;
        else                 g2->GetPoint(i+1,x2,dum);

        if (i==0)            x3=x1;
        else                 g2->GetPoint(i-1,x3,dum);

        Double_t tmp=y2;
        if (y1<y2) {y2=y1; y1=tmp;}
        //Double_t y3=1.;
        Double_t y3=y0;
        g3->SetPoint(i,x1,y3);

        Double_t binwl=(x1-x3)/2.;
        Double_t binwh=(x2-x1)/2.;
        if (binwl==0.)  binwl= binwh;
        if (binwh==0.)  binwh= binwl;
        g3->SetPointError(i,binwl,binwh,(y3-y2),(y1-y3));

    }
    return g3;

}

void myAddtoBand(TGraphErrors* g1, TGraphAsymmErrors* g2) {

    Double_t  x1=0., y1=0.,  y2=0., y0=0;
    //Double_t dx1=0.;
    //Double_t dum;

    if (g1->GetN()!=g2->GetN())
        std::cout << " graphs have not the same # of elements " << std::endl;

    Double_t* EYhigh = g2-> GetEYhigh();
    Double_t* EYlow  = g2-> GetEYlow();

    for (Int_t i=0; i<g1->GetN(); i++) {
        g1->GetPoint(i, x1,y1);
        g2->GetPoint(i, x1,y2);

        if ( y1==0 || y2==0 ) {
            std::cerr << "check these points very carefully : myAddtoBand() : point " << i << std::endl;
        }
        //    if (y1==0) y1=1;
        //    if (y2==0) y2=1;

        //    if (i==g1->GetN()-1) x2=x1;
        //    else                 g2->GetPoint(i+1,x2,dum);
        //    if (i==0)            x3=x1;
        //    else                 g2->GetPoint(i-1,x3,dum);

        Double_t eyh=0., eyl=0.;
        //if (y1<y2) {y2=y1; y1=tmp;}
        //Double_t y3=1.;

        //printf("%d: y1=%f y2=%f Eyhigh= %f Eylow= %f \n",i,y1,y2,EYhigh[i],EYlow[i]);

        y0=y1-y2;
        if (y0!=0) {
            if (y0>0){
                eyh=EYhigh[i];
                eyh=std::sqrt(eyh*eyh+y0*y0);
                //printf("high: %d: y0=%f eyh=%f  \n",i,y0,eyh);
                g2->SetPointEYhigh(i,eyh);
            } else {
                eyl=EYlow[i];
                eyl=std::sqrt(eyl*eyl+y0*y0);
                // printf("low: %d: y0=%f eyl=%f  \n",i,y0,eyl);
                g2->SetPointEYlow (i,eyl);
            }
        }
    }
    return;

}

TGraphErrors* TH1TOTGraph(TH1 *h1){


    if (!h1) std::cout << "TH1TOTGraph: histogram not found !" << std::endl;

    TGraphErrors* g1= new TGraphErrors();

    Double_t x, y, ex, ey;
    for (Int_t i=1 ; i<=h1->GetNbinsX(); i++) {
        y=h1->GetBinContent(i);
        ey=h1->GetBinError(i);
        x=h1->GetBinCenter(i);
        ex=h1->GetBinWidth(i);

        //   cout << " x,y = " << x << " " << y << " ex,ey = " << ex << " " << ey << endl;

        g1->SetPoint(i-1,x,y);
        g1->SetPointError(i-1,ex,ey);

    }

    //g1->Print();

    return g1;
}


void myTextF(Double_t x,Double_t y, const char *text, float size=0.06, Color_t color=1, int align=12) {

    //Double_t tsize=0.05;
    TLatex l; l.SetTextAlign(align);// l.SetTextSize(tsize);
    l.SetNDC();
    l.SetTextColor(color);
    l.SetTextFont(42);
    l.SetTextSize(size);
    l.DrawLatex(x,y,text);
}

void myText(Double_t x,Double_t y, const char *text, int size=23, Color_t color=1, int align=12) {
    //Double_t tsize=0.05;
    TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);
    l.SetTextAlign(align);
    l.SetNDC();
    l.SetTextColor(color);
    l.SetTextFont(43);
    l.SetTextSize(size);
    l.DrawLatex(x,y,text);
}


void myBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor,const char *text)
{

    Double_t tsize=0.06;

    TLatex l; l.SetTextAlign(12); //l.SetTextSize(tsize);
    l.SetNDC();
    l.DrawLatex(x,y,text);
    l.SetTextFont(42);

    Double_t y1=y-0.25*tsize;
    Double_t y2=y+0.25*tsize;
    Double_t x2=x-0.3*tsize;
    Double_t x1=x2-boxsize;

    printf("x1= %f x2= %f y1= %f y2= %f \n",x1,x2,y1,y2);

    TPave *mbox= new TPave(x1,y1,x2,y2,0,"NDC");

    mbox->SetFillColor(mcolor);
    mbox->SetFillStyle(1001);
    mbox->Draw();

    TLine mline;
    mline.SetLineWidth(4);
    mline.SetLineColor(1);
    mline.SetLineStyle(1);
    Double_t y_new=(y1+y2)/2.;
    mline.DrawLineNDC(x1,y_new,x2,y_new);

}

void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle, const char *text,Float_t msize, Int_t  size=26, double itv = 0.4)
{

    float pixels = gPad->GetWh()*gPad->GetAbsHNDC();
    float pixelsw = gPad->GetWw()*gPad->GetAbsWNDC();
    float tsize = size*1.0 / pixels;

    cout<<tsize<<" "<<pixels<<endl;

    TMarker *marker = new TMarker(x-(itv*tsize),y,8);
    marker->SetMarkerColor(color);  marker->SetNDC();
    marker->SetMarkerStyle(mstyle);
    marker->SetMarkerSize(msize);
    marker->Draw();

    TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize);
    l.SetTextFont(43);
    l.SetNDC();
    l.DrawLatex(x,y,text);
}


void myMarkerTextF(Double_t x,Double_t y,Int_t color,Int_t mstyle, const char *text,Float_t msize, double tsize=0.06, double itv = 0.4)
{
    TMarker *marker = new TMarker(x-(itv*tsize),y,8);
    marker->SetMarkerColor(color);  marker->SetNDC();
    marker->SetMarkerStyle(mstyle);
    marker->SetMarkerSize(msize);
    marker->Draw();

    TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize);
    l.SetNDC();
    l.SetTextFont(42);

    l.DrawLatex(x,y,text);
}

void findRangeGr(TGraphErrors* gr,  double rg0, double rg1,  double &rgLow, double &rgUp)
{
    
    
    
    int N = gr->GetN();
    double* x = gr->GetX();
    double* y = gr->GetY();
    double* ye = gr->GetEY();
    
    double mmin=0;
    double mmax=0;
    
    int xx=0;
    
    for(int k=0; k<N; k++){
        
        if( x[k] <rg0 || x[k] >rg1 ) continue;
        
        if( xx==0){
            mmax = y[k] + ye[k];
            mmin = y[k] - ye[k];
            xx++;
            
        }else{
            if( mmax< y[k] + ye[k] ) mmax = y[k] + ye[k];
            if( mmin> y[k] - ye[k] ) mmin = y[k] - ye[k];
        }
        //cout<< x[k] <<" "<<y[k] - ye[k] <<" "<< y[k] + ye[k] <<"  min "<<mmin<<" max "<<mmax<<endl;
        
    }
   // cout<<"Final " <<mmin<<" "<<mmax<<endl;
    rgUp = mmax;
    rgLow = mmin;
    
    
}
