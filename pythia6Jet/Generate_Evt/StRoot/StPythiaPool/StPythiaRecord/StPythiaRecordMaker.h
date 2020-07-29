#ifndef STPYTHIARECORDMAKER
#define STPYTHIARECORDMAKER

#include <iostream>

#include "TObject.h"
#include "TTree.h"
#include "TFile.h"
#include "TPythia6.h"
#include "TParticle.h"

#include "StSpinPool/StJetSkimEvent/StPythiaEvent.h"

using std::cout;
using std::endl;

class StPythiaRecordMaker : public TObject{
 public:
  StPythiaRecordMaker(){
    _skim = 0;
    record = new StPythiaEvent;
  }
  void SetFile(const char* fname){ file = fname; }
  void SetSkim(int skim) { _skim = skim;}
  int Init();
  int Make(int eventid, TPythia6 *pythia);
  int Finish();
 private:
  const char *file;
  TFile *fin;
  TTree *tree;
  StPythiaEvent *record;
  int _skim;
  ClassDef(StPythiaRecordMaker, 1);
};
#endif
