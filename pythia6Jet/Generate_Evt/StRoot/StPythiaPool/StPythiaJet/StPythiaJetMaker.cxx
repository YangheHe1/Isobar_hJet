#include "StPythiaJetMaker.h"
#include "StJetVector.h"
int StPythiaJetMaker::Init()
{
  fin = new TFile(file, "recreate");
  tree = new TTree("jet", "jetTree");
  for(size_t i = 0; i < jets.size(); i++){
    tree->Branch(jets[i]->branch(), "StPythiaJetEvent", jets[i]->jet());
  }
  return 1;
}

int StPythiaJetMaker::Make(int eventid, TPythia6 *pythia)
{
  for(size_t i = 0; i < jets.size(); i++){ 
    //set event id
    jets[i]->jet()->setid(eventid);
    //fill input particles to jet finder
    StProtoJet::FourVecList parjets;
    for (int ip = 0; ip < pythia->GetN(); ip++){
      bool pass = 0;
      StPythiaJetStatusCut *sc = jets[i]->statuscut();
      if(sc) pass = (*sc)(pythia, ip);
      if(pass) continue;
      int index = ip+1;
      StJetVector *vect = new StJetVector(pythia->GetK(index, 2),
					  pythia->GetK(index, 1),
					  pythia->GetK(index, 3),
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
      
      parjets.push_back(vect);
    }

    //jets found in jet finder
    StJetFinder::JetList protojets;
    jets[i]->jetfinder()->findJets(protojets, parjets);
    if(eventid % 1000 == 1){
       cout<<"eventid "<<eventid<<endl;
       cout<<"BEFORE proto size:"<<parjets.size()<<" jet size:"<<protojets.size()<<endl;
       cout<<"BRANCH "<<i<<" radius:"<<jets[i]->jetpars()->Rparam()<<endl;
       cout<<"AFTER proto size:"<<parjets.size()<<" jet size:"<<protojets.size()<<endl;
    }
    for(StJetFinder::JetList::const_iterator ijet = protojets.begin(); ijet != protojets.end(); ijet++){
      bool flag = 0;
      for(size_t icut = 0; icut < jets[i]->jetcut().size(); icut++){
	if((*(jets[i]->jetcut())[icut])(*ijet)) flag = 1;
      }
      if(flag) continue;
      //      jets[i]->jet()->fill(&(*ijet));
      Fill(jets[i]->jet(), &(*ijet));
      if(eventid % 1000 == 1) cout<<"eventId = "<<eventid<<"jet pt="<<(*ijet).pt()<<" found"<<endl;
    }
  }
  
  //  StPythiaJetCandidate *jc = (StPythiaJetCandidate *)(jet->getJets()->At(0));
  //  cout<<"jet pt="<<((StPythiaJetCandidate *)(jet->getJets()->At(0)))->pt()<<endl;
  //  cout<<"par px="<<((StPythiaJetParticle *)jc->getParticles()->At(0))->px()<<endl;
  tree->Fill();

  for(size_t i = 0; i < jets.size(); i++){
    jets[i]->jet()->clear();
  }

  return 1;
}
void StPythiaJetMaker::Fill(StPythiaJetEvent *jetevent, const StProtoJet *protojet)
{
  //TClonesArray *jets = jetevent->atJets();
  //  TClonesArray &array = *jets;
  //int index = jets->GetEntriesFast();
  //StPythiaJetCandidate *jetcnd = (StPythiaJetCandidate*)jets->ConstructedAt(index);
  //new(array[(int)jets->GetEntriesFast()]) StPythiaJetCandidate;
  StPythiaJetCandidate *jetcnd = jetevent->constructJet();
  
  jetcnd->setPx(protojet->px());
  jetcnd->setPy(protojet->py());
  jetcnd->setPz(protojet->pz());

  jetcnd->setPt(protojet->pt());
  jetcnd->setEta(protojet->eta());
  jetcnd->setPhi(protojet->phi());

  jetcnd->setArea(protojet->area());
  //cout<<"pt = "<<protojet->pt()<<" area = "<<protojet->area()<<endl;
  jetcnd->setAreaError(protojet->areaError());
  jetcnd->setE(protojet->e());
  jetcnd->setEt(protojet->eT());
  jetcnd->setMass(protojet->mass());
  jetcnd->setCharge(protojet->charge());

  //particles = jetcnd->atParticles();
  //    new TClonesArray("StPythiaJetParticle", 50);
  //  TClonesArray &array = *particles;
  //int index_p = particles->GetFastEntries();
  //StPythiaJetParticle *jetpar = (StPythiaJetParticle*)particles->ConstructAt(index_p);
  //cout<<"mult:"<<protojet->list().size()<<endl;
  for(size_t ii = 0; ii < protojet->list().size(); ii++){
    StPythiaJetParticle *jetpar = jetcnd->constructParticle();
    StJetVector* vec = (StJetVector*)(protojet->list()[ii]);
    //    cout<<"px="<<vec->px()<<endl;
    //    StPythiaJetParticle *p = new StPythiaJetParticle(vec);
    //    array->At(ii) = p;
    //    new(array[ii]) StPythiaJetParticle(vec);
    //    cout<<"ientry:"<<array.GetEntriesFast()<<endl;
    jetpar->setMom(vec->momentum());
    jetpar->setPos(vec->position());
    jetpar->setPdg(vec->pdg());
    jetpar->setStatus(vec->status());
    jetpar->setMother(vec->mother());
    jetpar->setDaughterfirst(vec->daughterfirst());
    jetpar->setDaughtersecond(vec->daughtersecond());
    jetpar->setCharge(vec->charge());
  }
  //  cout<<"entry:"<<particles->GetEntries()<<endl;  
  //  cout<<"entry:"<<array.GetEntries()<<endl;  

}
int StPythiaJetMaker::Finish()
{
  fin->Write();
  fin->Close();

  return 1;
}
