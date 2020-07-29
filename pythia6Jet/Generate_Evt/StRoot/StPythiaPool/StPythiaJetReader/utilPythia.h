#ifndef UTIL_PYTHIA
#define UTIL_PYTHIA

class StPythiaJetEvent;
class StPythiaJetCandidate;
class StPythiaJetParticle;
class TVector3;
class TParticle;

bool chargedHadron(const StPythiaJetParticle *par);
double parFrag(const StPythiaJetParticle *par, const TVector3 &vec);
TVector3 parLocalMomentum(const StPythiaJetParticle *par, const TVector3 &vec);
//TParticle
bool chargedHadron(const TParticle *par);
double parFrag(const TParticle *par, const TVector3 &vec);
TVector3 parLocalMomentum(const TParticle *par, const TVector3 &vec);

void PrintR(StPythiaJetCandidate* jet, StPythiaJetEvent *evntpar);
double FindRmin(StPythiaJetCandidate* jet, StPythiaJetEvent *evntpar);
StPythiaJetCandidate *FindParJet(StPythiaJetCandidate* jet, StPythiaJetEvent *evntpar, double range = 0.9);
double DeltaR(double etaA, float phiA, float etaB, float phiB);

#endif
