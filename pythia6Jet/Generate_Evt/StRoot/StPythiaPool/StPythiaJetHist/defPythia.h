#ifndef DEF_PYTHIA
#define DEF_PYTHIA

const int Npar = 2;
const int Ncone = 3;
const int Nproc = 3;
const int nbins = 15;
const double ptbins[nbins+1] = {5.0, 6.0, 7.0, 8.2, 9.6, 11.2, 13.1, 15.3, 17.9, 20.9, 24.5, 28.7, 33.6, 39.3, 46.0, 53.8};
const double parptbins[nbins+1] = {0, 6, 7.1, 8.3, 9.8, 11.6,
                             13.7, 16.1, 19.0, 22.4, 26.4,
                             31.1, 36.6, 43.2, 50.9, 60.0};
const int netabins = 2;
const double etabins[netabins+1] = {0, 0.5, 0.9};

int indexDet(double pt, double eta);
int indexPar(double pt, double eta);
double partonicWeight(double pt, double p0 = 1.987, double p1 = -0.677, double p2 = 0.265, double p3 = 0.759);
#endif
