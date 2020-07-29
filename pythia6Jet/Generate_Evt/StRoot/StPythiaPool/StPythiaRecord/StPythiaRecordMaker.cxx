#include "StPythiaRecordMaker.h"

int StPythiaRecordMaker::Init()
{
  if(!file) file = "pythia.root";
  fin = new TFile(file, "recreate");
  tree = new TTree("PythiaTree", "Pythia Record");
  tree->Branch("PythiaBranch", "StPythiaEvent", &record);

  return 1;
}
int StPythiaRecordMaker::Make(int eventid, TPythia6 *pythia)
{
  record->setRunId(1);
  record->setEventId(eventid);
  record->setProcessId(pythia->GetMSTI(1));
  record->setS(pythia->GetPARI(14));
  record->setT(pythia->GetPARI(15));
  record->setU(pythia->GetPARI(16));
  record->setPt(pythia->GetPARI(17));
  record->setCosTheta(pythia->GetPARI(41));
  record->setX1(pythia->GetPARI(33));
  record->setX2(pythia->GetPARI(34));
  record->setVertex(TVector3(0,0,0));
  record->setMstu72(pythia->GetMSTU(72));
  record->setMstu73(pythia->GetMSTU(73));
  record->setMstp111(pythia->GetMSTP(111));
  int Nlines = pythia->GetN();
  if(_skim == 1) Nlines = 8;
  for(int index = 1; index <= Nlines; index++){
    if(pythia->GetK(index, 2) == 88) continue;
    TParticle particle(pythia->GetK(index, 2),
		       pythia->GetK(index, 1),
		       pythia->GetK(index, 3),
		       -1,
		       pythia->GetK(index, 4),
		       pythia->GetK(index, 5),
		       pythia->GetP(index, 1),
		       pythia->GetP(index, 2),
		       pythia->GetP(index, 3),
		       pythia->GetP(index, 4),
		       pythia->GetV(index, 1),
		       pythia->GetV(index, 2),
		       pythia->GetV(index, 3),
		       pythia->GetV(index, 4));
  
    record->addParticle(particle);
  }
  tree->Fill();
  
  record->Clear();

  return 1;
}
int StPythiaRecordMaker::Finish()
{
  fin->Write();
  fin->Close();

  return 1;
}
