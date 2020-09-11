#ifndef __dNdpT__hh
#define __dNdpT__hh

class dNdpT
{
 public:
  dNdpT(int npar, int hjet, bool pTdepRAA=0); 
  ~dNdpT();

  double operator() (double *x, double *p);

  int fnpar;
  int hjettype;
  bool ffullspec;
  bool fpTdepRAA;
};

#endif
