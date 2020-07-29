//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 21 11:43:58 2020 by ROOT version 5.34/38
// from TTree MCTree/MCTree
// found on file: F_mixedevents_mult_0_vz_0_Psi_1.root
//////////////////////////////////////////////////////////

#ifndef readtree_h
#define readtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>


#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TMath.h"
#include "TFile.h"
#include "TH2.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TComplex.h"

#include <algorithm>
#include <TRandom.h>


using namespace std;

const Int_t MAXTRACK=10000;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class readtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        M_refmult;
   Double_t        M_Vz;
   Double_t        M_Vx;
   Double_t        M_Vy;
   Double_t        M_Psi2;
   Int_t           M_numTrk;
   Double_t        M_Pt[MAXTRACK];   //[M_numTrk]
   Double_t        M_Eta[MAXTRACK];   //[M_numTrk]
   Double_t        M_Phi[MAXTRACK];   //[M_numTrk]
   Double_t        M_Charge[MAXTRACK];   //[M_numTrk]
   Double_t        M_Dca[MAXTRACK];   //[M_numTrk]

   // List of branches
   TBranch        *b_M_refmult;   //!
   TBranch        *b_M_Vz;   //!
   TBranch        *b_M_Vx;   //!
   TBranch        *b_M_Vy;   //!
   TBranch        *b_M_Psi2;   //!
   TBranch        *b_M_numTrk;   //!
   TBranch        *b_M_Pt;   //!
   TBranch        *b_M_Eta;   //!
   TBranch        *b_M_Phi;   //!
   TBranch        *b_M_Charge;   //!
   TBranch        *b_M_Dca;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("F_mixedevents_mult_0_vz_0_Psi_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("F_mixedevents_mult_0_vz_0_Psi_1.root");
      }
      f->GetObject("MCTree",tree);

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

   fChain->SetBranchAddress("M_refmult", &M_refmult, &b_M_refmult);
   fChain->SetBranchAddress("M_Vz", &M_Vz, &b_M_Vz);
   fChain->SetBranchAddress("M_Vx", &M_Vx, &b_M_Vx);
   fChain->SetBranchAddress("M_Vy", &M_Vy, &b_M_Vy);
   fChain->SetBranchAddress("M_Psi2", &M_Psi2, &b_M_Psi2);
   fChain->SetBranchAddress("M_numTrk", &M_numTrk, &b_M_numTrk);
   fChain->SetBranchAddress("M_Pt", M_Pt, &b_M_Pt);
   fChain->SetBranchAddress("M_Eta", M_Eta, &b_M_Eta);
   fChain->SetBranchAddress("M_Phi", M_Phi, &b_M_Phi);
   fChain->SetBranchAddress("M_Charge", M_Charge, &b_M_Charge);
   fChain->SetBranchAddress("M_Dca", M_Dca, &b_M_Dca);
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
