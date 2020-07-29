#include "histResponse.h"

histResponse::histResponse(const char *name)
{
  hparticle = new TH2F(Form("%sParticleJetPtEta", name), ";p_{T} [GeV]; #eta", nbins, parptbins,netabins, etabins);
  hparton = new TH2F(Form("%sPartonJetPtEta", name), ";p_{T} [GeV]; #eta", nbins, parptbins,netabins, etabins);
  int N2d = (nbins+2)*(netabins+2);
  hresponse = new TH2F(Form("%sJetRes", name), ";particle jet (p_{T}, #eta) bin index;parton jet (p_{T}, #eta) bin index", N2d, 0, N2d, N2d, 0, N2d);
}
void histResponse::Add(const histResponse *hist, double w)
{
  hparticle->Add(hist->hparticle, w);
  hparton->Add(hist->hparton, w);
  hresponse->Add(hist->hresponse, w);
}
void histResponse::FillResponse(double pt1, double eta1, double pt2, double eta2, double weight){
  int aa = indexPar(pt1, eta1);
  int bb = indexPar(pt2, eta2);
  hresponse->Fill(aa, bb, weight);
}
