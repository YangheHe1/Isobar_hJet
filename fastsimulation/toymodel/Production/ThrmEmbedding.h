#ifndef __ThrmEmbedding__hh
#define __ThrmEmbedding__hh

#include "TObject.h"
#include "Rtypes.h"

class ThrmEmbedding : public TObject
{
 public:
  ThrmEmbedding();
  ThrmEmbedding(Int_t foundJet, Double_t pT, Double_t eta, Double_t phi, Double_t Area, Double_t deltaPt, Double_t pTreco, Double_t pTleading);
  void Set(Double_t pT, Double_t eta, Double_t phi, Double_t pTleading, Double_t area, Int_t N, Double_t pTreco, Double_t etaReco, Double_t phiReco, Double_t pTleadingReco, Double_t areaReco, Int_t Nreco, Double_t deltaPt);
  ~ThrmEmbedding();

  // EMBEDDING
  Int_t ffoundJet; 

  Double32_t fPtEmb;     //[0, 0, 16]
  Double32_t fetaEmb;    //[0, 0,  8]
  Double32_t fphiEmb;    //[0, 0,  8]
  Double32_t fAreaEmb;   //[0, 0,  8]
  Int_t 		 fNpartEmb;
  Double32_t fPtLeading; //[0, 0, 16]
  Double32_t fdeltapT;   //[0, 0, 16]
  Double32_t fPtEmbReco; //[0, 0, 16]
  Double32_t fetaEmbReco;    //[0, 0,  8]
  Double32_t fphiEmbReco;    //[0, 0,  8]
  Double32_t fPtLeadingReco; //[0, 0, 16]
  Double32_t fAreaReco;   //[0, 0,  8]
  Int_t 		 fNpartEmbReco;
 
  ClassDef(ThrmEmbedding, 3);
};

#endif
