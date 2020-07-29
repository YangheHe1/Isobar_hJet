//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 18 10:18:54 2020 by ROOT version 5.34/30
// from TTree PythiaTree/Pythia Record
// found on file: pythia6_pTHat8to30_1.root
//////////////////////////////////////////////////////////

#ifndef readtree_h
#define readtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>
#include <TVector3.h>
#include <TAttLine.h>
#include <TParticle.h>

#include "TH1.h"
#include "TH2.h"
#include <cmath>
#include "TProfile.h"
#include "TMath.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include <string>
#include <fstream>
#include <iostream>
#include "TProfile2D.h"

#include "TSystem.h"

using namespace std;
#include "TVector3.h"
#include "TLorentzVector.h"

#include <iostream>

using std::vector;
using std::string;

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxmParticles = 100000;

class readtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //StPythiaEvent   *PythiaBranch;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Int_t           mRunId;
   Int_t           mEventId;
   Int_t           mProcessId;
   Int_t           mTune;
   UInt_t          mVertex_fUniqueID;
   UInt_t          mVertex_fBits;
   Double_t        mVertex_fX;
   Double_t        mVertex_fY;
   Double_t        mVertex_fZ;
   Float_t         mS;
   Float_t         mT;
   Float_t         mU;
   Float_t         mPt;
   Float_t         mCosTheta;
   Float_t         mX1;
   Float_t         mX2;
   Int_t           mMstu72;
   Int_t           mMstu73;
   Int_t           mMstp111;
   Float_t         mPartonALL;
   Float_t         mDF1[34];
   Float_t         mDF2[34];
   Float_t         mF1[2];
   Float_t         mF2[2];
   Int_t           mParticles_;
   UInt_t          mParticles_fUniqueID[kMaxmParticles];   //[mParticles_]
   UInt_t          mParticles_fBits[kMaxmParticles];   //[mParticles_]
   Short_t         mParticles_fLineColor[kMaxmParticles];   //[mParticles_]
   Short_t         mParticles_fLineStyle[kMaxmParticles];   //[mParticles_]
   Short_t         mParticles_fLineWidth[kMaxmParticles];   //[mParticles_]
   Int_t           mParticles_fPdgCode[kMaxmParticles];   //[mParticles_]
   Int_t           mParticles_fStatusCode[kMaxmParticles];   //[mParticles_]
   Int_t           mParticles_fMother[kMaxmParticles][2];   //[mParticles_]
   Int_t           mParticles_fDaughter[kMaxmParticles][2];   //[mParticles_]
   Float_t         mParticles_fWeight[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fCalcMass[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPx[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPy[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPz[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fE[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVx[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVy[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVz[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fVt[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPolarTheta[kMaxmParticles];   //[mParticles_]
   Double_t        mParticles_fPolarPhi[kMaxmParticles];   //[mParticles_]

   // List of branches
   TBranch        *b_PythiaBranch_fUniqueID;   //!
   TBranch        *b_PythiaBranch_fBits;   //!
   TBranch        *b_PythiaBranch_mRunId;   //!
   TBranch        *b_PythiaBranch_mEventId;   //!
   TBranch        *b_PythiaBranch_mProcessId;   //!
   TBranch        *b_PythiaBranch_mTune;   //!
   TBranch        *b_PythiaBranch_mVertex_fUniqueID;   //!
   TBranch        *b_PythiaBranch_mVertex_fBits;   //!
   TBranch        *b_PythiaBranch_mVertex_fX;   //!
   TBranch        *b_PythiaBranch_mVertex_fY;   //!
   TBranch        *b_PythiaBranch_mVertex_fZ;   //!
   TBranch        *b_PythiaBranch_mS;   //!
   TBranch        *b_PythiaBranch_mT;   //!
   TBranch        *b_PythiaBranch_mU;   //!
   TBranch        *b_PythiaBranch_mPt;   //!
   TBranch        *b_PythiaBranch_mCosTheta;   //!
   TBranch        *b_PythiaBranch_mX1;   //!
   TBranch        *b_PythiaBranch_mX2;   //!
   TBranch        *b_PythiaBranch_mMstu72;   //!
   TBranch        *b_PythiaBranch_mMstu73;   //!
   TBranch        *b_PythiaBranch_mMstp111;   //!
   TBranch        *b_PythiaBranch_mPartonALL;   //!
   TBranch        *b_PythiaBranch_mDF1;   //!
   TBranch        *b_PythiaBranch_mDF2;   //!
   TBranch        *b_PythiaBranch_mF1;   //!
   TBranch        *b_PythiaBranch_mF2;   //!
   TBranch        *b_PythiaBranch_mParticles_;   //!
   TBranch        *b_mParticles_fUniqueID;   //!
   TBranch        *b_mParticles_fBits;   //!
   TBranch        *b_mParticles_fLineColor;   //!
   TBranch        *b_mParticles_fLineStyle;   //!
   TBranch        *b_mParticles_fLineWidth;   //!
   TBranch        *b_mParticles_fPdgCode;   //!
   TBranch        *b_mParticles_fStatusCode;   //!
   TBranch        *b_mParticles_fMother;   //!
   TBranch        *b_mParticles_fDaughter;   //!
   TBranch        *b_mParticles_fWeight;   //!
   TBranch        *b_mParticles_fCalcMass;   //!
   TBranch        *b_mParticles_fPx;   //!
   TBranch        *b_mParticles_fPy;   //!
   TBranch        *b_mParticles_fPz;   //!
   TBranch        *b_mParticles_fE;   //!
   TBranch        *b_mParticles_fVx;   //!
   TBranch        *b_mParticles_fVy;   //!
   TBranch        *b_mParticles_fVz;   //!
   TBranch        *b_mParticles_fVt;   //!
   TBranch        *b_mParticles_fPolarTheta;   //!
   TBranch        *b_mParticles_fPolarPhi;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pythia6_pTHat8to30_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pythia6_pTHat8to30_1.root");
      }
      f->GetObject("PythiaTree",tree);

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

   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_PythiaBranch_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_PythiaBranch_fBits);
   fChain->SetBranchAddress("mRunId", &mRunId, &b_PythiaBranch_mRunId);
   fChain->SetBranchAddress("mEventId", &mEventId, &b_PythiaBranch_mEventId);
   fChain->SetBranchAddress("mProcessId", &mProcessId, &b_PythiaBranch_mProcessId);
   fChain->SetBranchAddress("mTune", &mTune, &b_PythiaBranch_mTune);
   fChain->SetBranchAddress("mVertex.fUniqueID", &mVertex_fUniqueID, &b_PythiaBranch_mVertex_fUniqueID);
   fChain->SetBranchAddress("mVertex.fBits", &mVertex_fBits, &b_PythiaBranch_mVertex_fBits);
   fChain->SetBranchAddress("mVertex.fX", &mVertex_fX, &b_PythiaBranch_mVertex_fX);
   fChain->SetBranchAddress("mVertex.fY", &mVertex_fY, &b_PythiaBranch_mVertex_fY);
   fChain->SetBranchAddress("mVertex.fZ", &mVertex_fZ, &b_PythiaBranch_mVertex_fZ);
   fChain->SetBranchAddress("mS", &mS, &b_PythiaBranch_mS);
   fChain->SetBranchAddress("mT", &mT, &b_PythiaBranch_mT);
   fChain->SetBranchAddress("mU", &mU, &b_PythiaBranch_mU);
   fChain->SetBranchAddress("mPt", &mPt, &b_PythiaBranch_mPt);
   fChain->SetBranchAddress("mCosTheta", &mCosTheta, &b_PythiaBranch_mCosTheta);
   fChain->SetBranchAddress("mX1", &mX1, &b_PythiaBranch_mX1);
   fChain->SetBranchAddress("mX2", &mX2, &b_PythiaBranch_mX2);
   fChain->SetBranchAddress("mMstu72", &mMstu72, &b_PythiaBranch_mMstu72);
   fChain->SetBranchAddress("mMstu73", &mMstu73, &b_PythiaBranch_mMstu73);
   fChain->SetBranchAddress("mMstp111", &mMstp111, &b_PythiaBranch_mMstp111);
   fChain->SetBranchAddress("mPartonALL", &mPartonALL, &b_PythiaBranch_mPartonALL);
   fChain->SetBranchAddress("mDF1[34]", mDF1, &b_PythiaBranch_mDF1);
   fChain->SetBranchAddress("mDF2[34]", mDF2, &b_PythiaBranch_mDF2);
   fChain->SetBranchAddress("mF1[2]", mF1, &b_PythiaBranch_mF1);
   fChain->SetBranchAddress("mF2[2]", mF2, &b_PythiaBranch_mF2);
   fChain->SetBranchAddress("mParticles", &mParticles_, &b_PythiaBranch_mParticles_);
   fChain->SetBranchAddress("mParticles.fUniqueID", mParticles_fUniqueID, &b_mParticles_fUniqueID);
   fChain->SetBranchAddress("mParticles.fBits", mParticles_fBits, &b_mParticles_fBits);
   fChain->SetBranchAddress("mParticles.fLineColor", mParticles_fLineColor, &b_mParticles_fLineColor);
   fChain->SetBranchAddress("mParticles.fLineStyle", mParticles_fLineStyle, &b_mParticles_fLineStyle);
   fChain->SetBranchAddress("mParticles.fLineWidth", mParticles_fLineWidth, &b_mParticles_fLineWidth);
   fChain->SetBranchAddress("mParticles.fPdgCode", mParticles_fPdgCode, &b_mParticles_fPdgCode);
   fChain->SetBranchAddress("mParticles.fStatusCode", mParticles_fStatusCode, &b_mParticles_fStatusCode);
   fChain->SetBranchAddress("mParticles.fMother[2]", mParticles_fMother, &b_mParticles_fMother);
   fChain->SetBranchAddress("mParticles.fDaughter[2]", mParticles_fDaughter, &b_mParticles_fDaughter);
   fChain->SetBranchAddress("mParticles.fWeight", mParticles_fWeight, &b_mParticles_fWeight);
   fChain->SetBranchAddress("mParticles.fCalcMass", mParticles_fCalcMass, &b_mParticles_fCalcMass);
   fChain->SetBranchAddress("mParticles.fPx", mParticles_fPx, &b_mParticles_fPx);
   fChain->SetBranchAddress("mParticles.fPy", mParticles_fPy, &b_mParticles_fPy);
   fChain->SetBranchAddress("mParticles.fPz", mParticles_fPz, &b_mParticles_fPz);
   fChain->SetBranchAddress("mParticles.fE", mParticles_fE, &b_mParticles_fE);
   fChain->SetBranchAddress("mParticles.fVx", mParticles_fVx, &b_mParticles_fVx);
   fChain->SetBranchAddress("mParticles.fVy", mParticles_fVy, &b_mParticles_fVy);
   fChain->SetBranchAddress("mParticles.fVz", mParticles_fVz, &b_mParticles_fVz);
   fChain->SetBranchAddress("mParticles.fVt", mParticles_fVt, &b_mParticles_fVt);
   fChain->SetBranchAddress("mParticles.fPolarTheta", mParticles_fPolarTheta, &b_mParticles_fPolarTheta);
   fChain->SetBranchAddress("mParticles.fPolarPhi", mParticles_fPolarPhi, &b_mParticles_fPolarPhi);
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
