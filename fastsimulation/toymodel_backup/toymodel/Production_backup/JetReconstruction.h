#ifndef __JetReconstruction__hh
#define __JetReconstruction__hh

#include "TString.h"
#include "TClonesArray.h"

#include <vector>
#include <fastjet/config.h>            

#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/ClusterSequenceActiveArea.hh>
#include <fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/Subtractor.hh> 
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>

using std::vector;
using fastjet::PseudoJet;

class TFile;
class TTree;
class TClonesArray;

class JetReconstruction
{
 public:
  JetReconstruction(TString wrkdir, Double_t rparam, Double_t pTcut);
  ~JetReconstruction();
  
  void Run();

 private:
  void FillInputVector();
  void RecoJets();
  void RunEmbedding();
  void Reset();

 protected:
  // JET RECO CONFIG
  Double_t frparam;
  Double_t frparam_BG;
  Double_t fpTcut;
  Double_t fetaMax;
  
  // FASTJET
  vector<PseudoJet> kt_jets;
  vector<PseudoJet> akt_jets;
  vector<PseudoJet> input_data;

/*  
  FJWrapper *kt_data;
  FJWrapper *akt_data;
  FJWrapper *akt_data_emb;
*/
  // SAVED INFORMATION
  Double_t frho[3];
  Double_t fsigma[3];
  TClonesArray fktjetsarr;
  TClonesArray faktjetsarr;
  TClonesArray fembeddingarr;
  TClonesArray *fpartarr;

  // I/O
  TFile *finput;
  TFile *foutput;
  TTree *ftreetoy;
  TTree *ftreejets;
};

#endif
