#include "TF1.h"
#include "TFile.h"
#include "TH3D.h"
//==============================================================
//efficiency functions
//==============================================================

//----------------------------------------------
//general functions
//----------------------------------------------
Double_t Eff_track_rec_function(Double_t* x,Double_t* par)
{
    // Track reconstruction efficiency parametrization
    Double_t pt,y;
    Double_t A,B,C;A=par[0];B=par[1];C=par[2];
    pt=x[0];
    y=A*(exp(-pow(B/pt,C)));
    return y;
}

//----------------------------------------------

Double_t RatioTsalis(Double_t* x_val, Double_t* par)
{
   double A1=par[0];
	double B1=par[1];
	double n1=par[2];
   double T1=par[3];
	double A2=par[4];
	double B2=par[5];
	double n2=par[6];
   double T2=par[7];
   double pT=x_val[0];

	double denom=A1*TMath::Power(B1,pT)*TMath::Power((1.+pT/(n1*T1)),-n1);
	double num=A2*TMath::Power(B2,pT)*TMath::Power((1.+pT/(n2*T2)),-n2);
   double y=num/denom;
   return y;
}
//----------------------------------------------

Double_t FitFunc(Double_t* x_val, Double_t* par)
{
   double A=par[0];
   double B=par[1];
   double C=par[2];
   double pT=x_val[0];

   double y=A*TMath::Exp(-TMath::Power(B/pT,C));
   return y;
}

//----------------------------------------------
//efficiencies
//----------------------------------------------
//
//----------------------------------------------
double efficiency11(double pt, bool central=1, TString ratio="AuAu", short cutset=1, double increment=0)
{

	const int npar=6; //number of parameters
	const int ncuts=3; //number of track cut sets

	//parameter values for different systems and sets of track cuts
	//central:0 hadron ratio:AuAu cutset:1-2 
	double par_peri_ratAuAu[ncuts][npar]={
		{0.964565,0.191389,3.0984,1.0283,0.000256939,0.285302}, //cut set 1
		{0.945533,0.191465,3.0695,1.01543,0.00116032,0.349932}, //cut set 2
		{0.99967,0.162127,4.63077,0.994843,0.0322634,4.61138} //cut set3 = cut set1 with global tracks	
	};

	//central:1 ratio:AuAu cutset:1-2 
	double par_cent_ratAuAu[ncuts][npar]={
		{0.809761,0.206039,1.90741,0.967221,0.00161199,0.193891},
		{0.764473,0.207159,1.97916,0.914462,0.00332684,0.225286},
		{0.902374,0.18195,5.50839,0.882829,0.00203566,2.37983} //cut set3 = cut set1 with global tracks
	};

	//central:0 ratio:pp cutset:1-2 
	double par_peri_ratpp[ncuts][npar]={
		{0.947572,0.193031,3.13351,1.03095,0.000357085,0.290538},
		{0.928739,0.193316,3.11881,1.0081,0.00305566,0.417085},
		{0.999824,0.162929,4.54422,0.994827,0.0315042,4.61133} //cut set3 = cut set1 with global tracks  
	};

	//central:1 ratio:pp cutset:1-2 
	double par_cent_ratpp[ncuts][npar]={
		{0.811528,0.208149,1.80443,0.913892,0.00415818,0.264741},
		{0.764165,0.208713,1.88208,0.866476,0.00809367,0.308728},
		{0.902645,0.184137,5.55228,0.883398,0.00199387,2.37534} //cut set3 = cut set1 with global tracks 
	};

	//fill parameter array
	double par[npar];
	for(int i=0; i<npar;i++)
	{
		if(ratio=="AuAu")
		{
			if(central) par[i]=par_cent_ratAuAu[cutset-1][i];
			else par[i]=par_peri_ratAuAu[cutset-1][i];
		}
		else if(ratio=="pp")
		{

			if(central) par[i]=par_cent_ratpp[cutset-1][i];
			else par[i]=par_peri_ratpp[cutset-1][i];
		}
	}//loop over parameters	

	float pTx=0.6;
	if(cutset==3 && central==1) pTx=2.0;
	TF1* f_Efficiency1 = new TF1("f_Efficiency1",Eff_track_rec_function,0,pTx,3);
	TF1* f_Efficiency2 = new TF1("f_Efficiency2",Eff_track_rec_function,pTx,8.0,3);
	
	f_Efficiency1->SetParameters(par[0],par[1],par[2]);
	f_Efficiency2->SetParameters(par[3],par[4],par[5]);
	
	double epsilon=0;
	if(pt<pTx)
		epsilon=f_Efficiency1->Eval(pt);
	else
		epsilon=f_Efficiency2->Eval(pt);
	delete f_Efficiency1;
	delete f_Efficiency2;
	epsilon=epsilon+increment;
	if (epsilon>1) epsilon=1.0;
	return epsilon;
}

//----------------------------------------------
//OLD - run11 TPC tracking efficiency from Stephen
Double_t efficiencyStephen(Double_t pt, TF1* effLow, TF1* effHigh, Double_t increment,bool kcentral)
{
   Double_t eff;
   if(pt<=1.2)eff = effLow->Eval(pt);
   else eff = effHigh->Eval(pt);
   eff=eff*(24.0/22.0); //efficiency was calculated as N_reco(in 22 TPC sectors)/N_emb(24 sectors), but we need it as  N_reco(22)/N_emb(22)
	if(eff<=0)eff=0.0001;
	Double_t scl=1;
	if(!kcentral)
	{
		Double_t A,B,C,D,E;
		A=5.41860e-01;
		B=4.98016;
		C= -9.63255e-02;
		D=-1.74537e-03;
		E=1.22717;
		if(pt<1)
    		scl=A*(exp(-B*pt+C))+D*pt+E;
		else
		 	scl=D*pt+E;
	}
	eff=eff*scl;
	eff=eff+increment;
return eff;
}

//---------------------------------------------
//OLD - run11 TPC tracking efficiency from Alex
Double_t efficiencyAlex(Double_t pt,bool kcentral=1,TString ratio="AuAu",Double_t increment=0.0)
{
	//if(pt>9.9)pt=9.9;
	if(pt>4.9)pt=4.9;
	TF1* f_Efficiency = new TF1("f_EfficiencyCent",Eff_track_rec_function,0,5.0,3);
	if(kcentral)
	{
		if(ratio=="pp") //pp like ratio of hadrons (pi/K/p)
			f_Efficiency->SetParameters(7.45643e-01,1.43725e-01,2.02904e+00);
		else //AuAu like ratio of hadrons (pi/K/p)
			f_Efficiency->SetParameters(7.76285e-01,1.66520e-01,1.66966e+00);
	}
	else
	{
		if(ratio=="pp")
			f_Efficiency->SetParameters(9.06946e-01,1.45242e-01,2.87409e+00);
		else
			f_Efficiency->SetParameters(8.51638e-01,1.65388e-01,1.95871e+00);
	}
	Double_t eff=f_Efficiency->Eval(pt);
	eff=eff+increment;
	delete f_Efficiency;
return eff;
}

double efficiencyTOFvBEMC(double pt,double increment)
{
	if(pt>25.0) pt=25.0;

	TF1* fitr1 = new TF1("fitr1",FitFunc,0.1,5.1,3);
	fitr1->SetNpx(10000);
	fitr1->SetParameters(9.18693e-01,1.03313e-03,1.85958e-01);

	TF1* fitr2 = new TF1("fitr2",RatioTsalis,5,26,8);
	fitr2->SetNpx(10000);
	fitr2->SetParameters(1.21012e+10,1.64290e+00,1.43554e+01,1.70386e-01,8.57414e+09,1.57782e+00,1.40643e+01,1.72733e-01);

	float epsilon=1;
	if(pt<=5)
		epsilon=fitr1->Eval(pt);
	else
		epsilon=fitr2->Eval(pt);

	epsilon=epsilon+increment;
	delete fitr1;
	delete fitr2;
	return epsilon;	

}
//======================================
//momentum resolution (from hadron embedding)
//======================================
Double_t mom_res(double pT, short type)
{
	if(pT<0) return 0;
	double sigma=0;
	if(type==0) sigma=0.01*pT*pT;
	else if(type==1) sigma=0.005*pT*pT;
	else if (type==2)sigma=0.036+0.00047*pT+0.0028*pT*pT;
	else if (type==3)sigma=0.003*pT*pT;
	if(sigma<0) sigma=0.005*pT*pT;
	return sigma;
}
/*
Double_t efficiency3d(double pT, double Eta, double Phi, bool kcentral=1,TH3D* h3dEff=NULL ){
	
	int xbin = h3dEff->GetXaxis()->FindBin(pT);
	int ybin = h3dEff->GetYaxis()->FindBin(Eta);
	int zbin = h3dEff->GetZaxis()->FindBin(Phi);
	int bin = h3dEff->GetBin(xbin,ybin,zbin);
	double eff = h3dEff->GetBinContent(bin);

	return eff;
}*/
