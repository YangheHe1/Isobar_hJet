#include "dNdpT.h"
#include "TMath.h"

dNdpT::dNdpT(int npar, int hjet, bool pTdepRAA)
{
  fnpar = npar;
  hjettype = hjet;
  fpTdepRAA=pTdepRAA;
}

dNdpT::~dNdpT()
{
}

double dNdpT::operator() (double *x, double *p)
{
  // dNdpT = fb*B(pT) + fs*S(pT)*C(pT)
  // B(pT) = b^2 * pT * exp(-b*pT)
  // S(pT) = 1/pT^n
  // C(pT) = 1 / {1 + exp[-(pT - a)/R]}
  
  // par[0] = fb
  // par[1] = b
  // par[2] = fs
  // par[3] = n
  // par[4] = a 
  // par[5] = R
  
  double pT = *x;
  double par[fnpar];

  for(int i = 0; i < fnpar; i++)
    par[i] = p[i];

  Double_t J;
  Double_t C = 1/(1 + TMath::Exp(-(pT - 2.8)/0.3)); //cutoff function
  Double_t A=par[0]; //amplitude
  Double_t RAAreq=par[fnpar-1]; //required RAA
  Double_t RAA=RAAreq; //final RAA (can be pT dependent)
//pT-dependent RAA
  if(fpTdepRAA)
	{
		Double_t RAA_lowpT=0.2;
		if(RAAreq<RAA_lowpT)RAA_lowpT=RAAreq;
		if(pT<5)RAA=RAA_lowpT;
		//else if(pT<10)RAA=1.0;
		else if(pT<15)RAA=(RAA_lowpT-RAAreq)*(15-pT)/10+RAAreq;
		else RAA=RAAreq;
	}

//POWER-LAW JET SPECTRUM  
	if(hjettype==0)
	{
	  J = 1./TMath::Power(pT, par[1]); 
	  return A*J; //hard jet spectrum
	}

//double Levy fit to Pythia (full) jet spectrum
	else if(hjettype==1) 
	{
	    Double_t y1,y2, B1, T1, n1, m1, mu1,B2, T2, n2, m2, mu2;
	    B1    = par[1];
   	 T1    = par[2];
	    n1    = par[3];
	    m1   = par[4];
		 mu1  = par[5];
		 B2    = par[6];
   	 T2    = par[7];
	    n2    = par[8];
   	 m2   = par[9];
		 mu2   = par[10];

	    Double_t mT1 = TMath::Sqrt((pT-mu1)*(pT-mu1)+m1*m1);
		 Double_t mT2 = TMath::Sqrt((pT-mu2)*(pT-mu2)+m2*m2);
		 
		y1 =B1/TMath::Power(1.0+(mT1-m1)/(n1*T1),n1);
		y2=B2/TMath::Power(1.0+(mT2-m2)/(n2*T2),n2);
		
		if(pT<15)
			J=y1;
		else if(pT<25)
		{
			Double_t c1=(pT-15.0)/10.0;
			J=((1-c1)*y1+c1*y2);
		}
	   else J=y2;
	
		if(J<0)J=0;
		return A*J*RAA; //hard jet spectrum
	}//hjettype==1 

	//Tsalis fit to pythia (full) jet spectrum
	else if(hjettype==2)
	{
	   double n=par[1];
	   double T=par[2];
   	
	   double J=pT*TMath::Power((1.+pT/(n*T)),-n);
		return A*J*RAA; //hard jet spectrum
	}

}//dNdpT
