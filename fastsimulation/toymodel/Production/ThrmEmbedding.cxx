#include "ThrmEmbedding.h"

ClassImp(ThrmEmbedding)

//____________________________________________________________________________
ThrmEmbedding::ThrmEmbedding()
{
  fPtEmb = 0;
  fetaEmb = 0;
  fphiEmb = 0;
  fAreaEmb = 0;
  fPtLeading = 0;
  fNpartEmb = 0;
  fPtEmbReco = 0;
  fetaEmbReco = 0;
  fphiEmbReco = 0;
  fPtLeadingReco = 0;
  fNpartEmbReco = 0;
  fAreaReco = 0;
  fdeltapT = 0;
  ffoundJet = 0;
}

//____________________________________________________________________________
ThrmEmbedding::ThrmEmbedding(Int_t foundJet, Double_t pT, Double_t eta, Double_t phi, Double_t Area, Double_t deltaPt, Double_t pTreco, Double_t pTleading)
{
  ffoundJet = foundJet;
  fPtEmb = pT;
  fetaEmb = eta ;
  fphiEmb = phi;
  fAreaEmb = Area;
  fdeltapT = deltaPt;
  fPtEmbReco = pTreco;
  fPtLeading = pTleading;
}

//____________________________________________________________________________
ThrmEmbedding::~ThrmEmbedding()
{
}

//____________________________________________________________________________
void ThrmEmbedding::Set(Double_t pT, Double_t eta, Double_t phi, Double_t pTleading, Double_t area, Int_t N,  Double_t pTreco, Double_t etaReco, Double_t phiReco, Double_t pTleadingReco, Double_t areaReco, Int_t Nreco, Double_t deltaPt)
{
  fPtEmb = pT; 
  fetaEmb = eta;
  fphiEmb = phi;
  fPtLeading = pTleading;
  fAreaEmb = area;
  fNpartEmb = N;
  fPtEmbReco = pTreco;
  fetaEmbReco = etaReco;
  fphiEmbReco = phiReco;
  fPtLeadingReco = pTleadingReco;
  fAreaReco = areaReco;
  fNpartEmbReco = Nreco;
  fdeltapT = deltaPt;
}


