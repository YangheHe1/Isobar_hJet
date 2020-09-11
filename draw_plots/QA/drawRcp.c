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
#include <fstream>

#include <iostream>


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

TH1D *hPt_C;
TH1D *hPt_P;
TH1D *hNC;
TH1D *hNP;

TH1D *hRatio;

TGraphErrors *gr_rcp;

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
int clr[] = {kBlue,kRed+1, kBlue, 1,  kGreen+1,kGray+2, kCyan+1, kYellow+1,kViolet+1};

const char* type[]={"same event", "mix event", "^{4}He+^{14}N 1 PeV", "^{4}He+^{14}N 10 PeV"};

void readin(){
    
    //===============================  read in =====================================
    sprintf(name,"Step2/Step2_QA_trig_badrun_pile.root");
    cout<<name<<endl;
    fin = TFile::Open(name);
    
    sprintf(name,"Pthist_C");
    hPt_C     = (TH1D*)fin->Get(name);
    
    sprintf(name,"eventC");
    hNC     = (TH1D*)fin->Get(name);
    
    //===============================  read in =====================================
    
    sprintf(name,"Pthist_P");
    hPt_P     = (TH1D*)fin->Get(name);
    
    sprintf(name,"eventP");
    hNP     = (TH1D*)fin->Get(name);

    float nc = hNC -> GetEntries();
    float np = hNP -> GetEntries();
    
    cout<<"central "<<nc<<" Peripheral "<<np<<endl;
 
    hPt_C -> Scale(1./nc);
    hPt_P -> Scale(1./np);

	//hPt_C ->Rebin(2);
        //hPt_P ->Rebin(2);
    hPt_C -> Scale(1./359.05);
    hPt_P -> Scale(1./11.85);

/*
   TH1D *hp = (TH1D*) hPt_P->Clone();
    hPt_C ->Divide(hp);
    hRatio = (TH1D*) hPt_C ->Clone();
*/
    hRatio = (TH1D*) hPt_C ->Clone();
    Int_t binmax=hPt_C ->FindLastBinAbove(0,1);
    for(int ibin=1;ibin<=binmax;ibin++){
        double y = hPt_C ->GetBinContent(ibin);
        double yref = hPt_P ->GetBinContent(ibin);
        if(y==0||yref==0)
			{hRatio->SetBinContent(ibin,0);
			hRatio->SetBinError(ibin,0);
				continue;}
        double binerr = hPt_C ->GetBinError(ibin);
        double binerref =hPt_P ->GetBinError(ibin);
        hRatio->SetBinContent(ibin,y/yref);
        double err = sqrt( binerr*binerr/(yref*yref)+ binerref*binerref*y*y/pow(yref,4) );
        hRatio->SetBinError(ibin,err);
	if(ibin==10) cout<<"value : "<<y/yref<<"  error : "<<err<<endl;
    }
  //hRatio->Rebin(4);
  cout<<"Error(bin10): "<<hRatio->GetBinError(10)<<endl;
  cout<<"Content(bin10): "<<hRatio->GetBinContent(10)<<endl;
   hRatio->Draw(); 

    double xx[15],yy[15],x_err[15],x_err2[15],y_err[15],y_err2[15];
    int nn=15;
    ifstream f_au("AuAu_Rcp.txt");
    for(int i=0;i<15;i++){
        f_au>>xx[i]>>yy[i]>>x_err[i]>>x_err2[i]>>y_err[i]>>y_err2[i];
    }

    f_au.close();

    gr_rcp =  new TGraphErrors(nn,xx,yy,x_err,y_err);


}

void drawRatio(){
    
    int nx = 1;
    int ny = 1;
    
    double tx0=0.1, ty0=0.1;
    float tsize = 0.05;
    
    sprintf(name, "RhoCompar");
    
    can[0]= new TCanvas(name,"Graph",10,10,1100,900);    
	pad[0][0]=new TPad("pad1","pad1",0.06,0.06,1,1);
    pad[0][0]->Draw();
	
    
    pad[0][0]->SetTopMargin(0.08);
	pad[0][0]->SetBottomMargin(0.08);
    pad[0][0]->cd();
            
            
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    //gPad->SetLogy(1);
    gStyle->SetOptStat(0);      //remove the entries,mean,RMS in the upper right.
    gStyle->SetOptTitle(0);
            
    hhtem[0] =  (TH1D*)    hRatio  ->Clone();
    grtmpErr0[0 ]= (TGraphErrors*)  gr_rcp ->Clone();
            
    
    hhtem[0] -> SetLineColor(clr[0]);
    hhtem[0] -> SetMarkerStyle(20);
    hhtem[0] -> SetMarkerColor(clr[0]);

    grtmpErr0[0]-> SetLineColor(clr[1]);
    grtmpErr0[0] -> SetMarkerStyle(24);
    grtmpErr0[0] -> SetMarkerColor(clr[1]);

            char *xtile ="P_{T} (GeV/c)";
            //char *ytile ="R_{cp}=Yeild^{0-10%}/Yeild^{60-80%}";
            char *ytile ="R_{cp}=(Yeild/<N_{binary}>)^{0-10%}/(Yeild/<N_{binary}>)^{60-80%}";
            
            //hhtem[0] -> GetXaxis()->SetTitle(xtile);
            hhtem[0] -> GetYaxis()->SetTitle(ytile);
            //hhtem[0] -> GetZaxis()->SetTitle(ztile);
            
            hhtem[0] -> GetXaxis()->SetRangeUser(0,35);
            //hhtem[0] -> GetYaxis()->SetRangeUser(0,5);
            hhtem[0] -> GetXaxis()->SetNdivisions(507);
            hhtem[0] -> GetYaxis()->SetNdivisions(507);
            
            hhtem[0] -> GetYaxis()->SetTitleOffset(1.2);
            //hhtem[1]->SetTitleSize(0.05);
            hhtem[0]->GetYaxis()->CenterTitle(true);
            
            
            hhtem[0]->DrawClone("P");
            grtmpErr0[0]->DrawClone("Psame");
            //hhtem[1]->DrawClone("same");
            
            tx0=0.3, ty0=0.77;
            //myTextF(tx0,ty0,"p+Pb 5.02 TeV",tsize*0.8,1,12);
            
            tx0=0.3;ty0=0.72;
            //myTextF(tx0,ty0,"centrality 60-80%",tsize*0.8,1,12);
            
            tx0=0.65, ty0=0.8;
            //myTextF(tx0,ty0,"centrality 0-10%",tsize*0.8,1,12);
            
            
                float _yy = 0*0.08 + 0.65;
                leg = mylegF(0.6,_yy,0.8,0.8,0.05);
                leg->AddEntry(hhtem[0],"Zr+Zr","lp");
                leg->AddEntry(grtmpErr0[0],"Au+Au,#pi^{+}+#pi^{-}","lp");
                leg->Draw("same");
            
            
        
    
    TLine *tl = new TLine(0,1,35,1);
	tl->Draw("same");    
            
   
    
}

 void drawRcp(){
    
    readin();
    
    drawRatio();
    
}
   
