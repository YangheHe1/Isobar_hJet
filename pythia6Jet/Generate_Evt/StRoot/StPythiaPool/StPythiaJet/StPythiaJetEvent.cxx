#include "StPythiaJetEvent.h"
//#include <iostream>
//using namespace std;
/*
void StPythiaJetEvent::fill(const StProtoJet *pj)
{
  TClonesArray &array = *jets;
  new(array[(int)jets->GetEntriesFast()]) StPythiaJetCandidate(pj);
  //  cout<<jets->GetEntriesFast()<<endl;
}
*/
//
void StPythiaJetEvent::clear()
{
  int entries = jets->GetEntriesFast();
  for(int ii = 0; ii < entries; ii++){
	StPythiaJetCandidate *jetcnd = (StPythiaJetCandidate *)jets->At(ii);
	jetcnd->clear();
  }
  jets->Clear();
}
  
