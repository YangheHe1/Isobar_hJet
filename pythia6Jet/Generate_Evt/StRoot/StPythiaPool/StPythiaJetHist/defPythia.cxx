#include"defPythia.h"

#include "TAxis.h"
#include "TMath.h"
double partonicWeight(double pt, double p0, double p1, double p2, double p3)
{
  double wght = 0.;
  wght = 1./(1+(p0+p1*(pt-2)+p2*TMath::Power(pt-2, 2))*TMath::Exp(-p3*(pt-2)));
  return wght;
}
int indexDet(double pt, double eta)
{
  TAxis xdet(nbins, ptbins);
  int x1 = xdet.FindBin(pt);
  TAxis ydet(netabins, etabins);
  int y1 = ydet.FindBin(eta);
  return x1+y1*(xdet.GetNbins()+2);
}
int indexPar(double pt, double eta)
{
  TAxis xdet(nbins, parptbins);
  int x1 = xdet.FindBin(pt);
  TAxis ydet(netabins, etabins);
  int y1 = ydet.FindBin(eta);
  return x1+y1*(xdet.GetNbins()+2);
}

