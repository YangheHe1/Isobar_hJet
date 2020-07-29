#ifndef STPYTHIAJETEVENT
#define STPYTHIAJETEVENT

#include "TObject.h"
//#include "StJetFinder/StProtoJet.h"
#include "TClonesArray.h"
#include "StPythiaJetCandidate.h"

class StPythiaJetEvent : public TObject{
 public:
  StPythiaJetEvent(){
    _id = -1;
    jets = new TClonesArray("StPythiaJetCandidate", 50);
  }
  ~StPythiaJetEvent(){
    if(jets) delete jets;
  }
  TClonesArray *getJets() const{ return jets;}
  StPythiaJetCandidate *constructJet() { 
    int ii = jets->GetEntriesFast();
    return (StPythiaJetCandidate*)jets->ConstructedAt(ii);
  }
//  void fill(const StProtoJet *pj);
  void setid(int ii) { _id = ii; }
  int eventId() const { return _id; }
   void clear();
 private:
  int _id;
  TClonesArray *jets;
  ClassDef(StPythiaJetEvent, 1);
};
//ClassImp(StPythiaJetEvent);
#endif
