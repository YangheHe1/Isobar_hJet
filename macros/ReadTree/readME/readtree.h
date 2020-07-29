//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 22 09:05:53 2020 by ROOT version 5.34/38
// from TTree METree/METree
// found on file: MCTree_029.root
//////////////////////////////////////////////////////////

#ifndef readtree_h
#define readtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

#include "TTree.h"
#include "TSystem.h"

using namespace std;

#include "TProfile.h"
#include <string>
#include <TH2.h>
#include <TH1.h>

#include <fstream>
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"

#include <iostream>

using std::vector;
using std::string;

const Int_t MAXTRACK=10000;
// Fixed size dimensions of array or collections stored in the TTree if any.

class readtree {
public :

   readtree(string filelist, int fr, int tr, int job);
   void InitHist();
   void TrkLoop();
   

   int from;
   int to;
   int Njob;
   char name[200];

   TFile *fout;
   TH1D *Hntrk;
   TH1D *Hntrk_class;
   TH1D *Hvz;
   TH1D *Hpt;
   TH1D *Heta;
   TH1D *Hphi;
   TH1D *Hcharge;
   TH2D *Hphi_vs_eta;

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        ME_refmult;
   Int_t           ME_numTrk;
   Double_t        ME_Vz;
   Double_t        ME_Vx;
   Double_t        ME_Vy;
   Double_t        ME_Psi2;
   Double_t        ME_Pt[MAXTRACK];   //[ME_numTrk]
   Double_t        ME_Eta[MAXTRACK];   //[ME_numTrk]
   Double_t        ME_Phi[MAXTRACK];   //[ME_numTrk]
   Double_t        ME_Charge[MAXTRACK];   //[ME_numTrk]
   Double_t        ME_Dca[MAXTRACK];   //[ME_numTrk]

   // List of branches
   TBranch        *b_ME_refmult;   //!
   TBranch        *b_ME_numTrk;   //!
   TBranch        *b_ME_Vz;   //!
   TBranch        *b_ME_Vx;   //!
   TBranch        *b_ME_Vy;   //!
   TBranch        *b_ME_Psi2;   //!
   TBranch        *b_ME_Pt;   //!
   TBranch        *b_ME_Eta;   //!
   TBranch        *b_ME_Phi;   //!
   TBranch        *b_ME_Charge;   //!
   TBranch        *b_ME_Dca;   //!

   readtree(TTree *tree=0);
   virtual ~readtree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef readtree_cxx
readtree::readtree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MCTree_029.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MCTree_029.root");
      }
      f->GetObject("METree",tree);

   }
   Init(tree);
}

readtree::~readtree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t readtree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t readtree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void readtree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ME_refmult", &ME_refmult, &b_ME_refmult);
   fChain->SetBranchAddress("ME_numTrk", &ME_numTrk, &b_ME_numTrk);
   fChain->SetBranchAddress("ME_Vz", &ME_Vz, &b_ME_Vz);
   fChain->SetBranchAddress("ME_Vx", &ME_Vx, &b_ME_Vx);
   fChain->SetBranchAddress("ME_Vy", &ME_Vy, &b_ME_Vy);
   fChain->SetBranchAddress("ME_Psi2", &ME_Psi2, &b_ME_Psi2);
   fChain->SetBranchAddress("ME_Pt", ME_Pt, &b_ME_Pt);
   fChain->SetBranchAddress("ME_Eta", ME_Eta, &b_ME_Eta);
   fChain->SetBranchAddress("ME_Phi", ME_Phi, &b_ME_Phi);
   fChain->SetBranchAddress("ME_Charge", ME_Charge, &b_ME_Charge);
   fChain->SetBranchAddress("ME_Dca", ME_Dca, &b_ME_Dca);
   Notify();
}

Bool_t readtree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void readtree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t readtree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef readtree_cxx
