#ifndef READ_PYTHIA
#define READ_PYTHIA

#include "TObject.h"
#include "StPythiaPool/StPythiaJetHist/defPythia.h"

class TFile;
class histPythia;
class histResponse;
//
class StPythiaEvent;
class StPythiaJetEvent;
class StPythiaJetCandidate;
class readPythia : public TObject
{
 public:
  readPythia(StPythiaEvent *, StPythiaJetEvent *parton, StPythiaJetEvent *particle);
  int Make(int iEvent);
  int Init();
  int Finish();
  void SetFile(const char *file) { mFilename = file; }
 protected:
  double UnderlyingEventCone(const StPythiaJetCandidate *jet, bool flag_match, double Rcone = 0.5);
  int ProcessIndex();

 private:
  StPythiaEvent *mPythia;
  StPythiaJetEvent *mPythiaJet[Npar];
  TFile *mFile;
  const char *mFilename;
  histPythia *mHist;
  histResponse *mHistResponse;
  ClassDef(readPythia, 0);
};
#endif
