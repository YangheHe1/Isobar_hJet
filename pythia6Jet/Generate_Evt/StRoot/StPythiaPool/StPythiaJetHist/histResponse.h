#ifndef HIST_RESPONSE
#define HIST_RESPONSE

#include "TObject.h"

#include "defPythia.h"

#include "TH2F.h"

class histResponse : public TObject
{
 public:
  histResponse(){}
  histResponse(const char *name);
  void Add(const histResponse* hist, double w = 1.0);
  void FillResponse(double pt1, double eta1, double pt2, double eta2, double weight = 1.0);
  TH2F *hparticle;
  TH2F *hresponse;
  TH2F *hparton;
  
 private:
  ClassDef(histResponse, 1);

};
#endif
