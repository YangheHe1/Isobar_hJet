//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 13 03:07:36 2020 by ROOT version 5.34/38
// from TTree EventTree/EventTree
// found on file: B74F2EA5ADE63E0429C2ECD3ABDDB0DA_449.root
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
const int NCENT=9;

// Fixed size dimensions of array or collections stored in the TTree if any.

class readtree {
public :

   readtree(string filelist, int fr, int tr);
   void InitHist();
   void TrkLoop();
   bool refmult_check(int __nBTOFMatch, int __refMult);
   bool ntrk_check(int __nBTOFMatch, int __ntrk);
   int centrality(int _refMult);

   int from;
   int to;
   int Njob;
   char name[200];

   double NPt;


   TFile *fout;
   TH1D *Hvz;
   TH1D *Hrefmult;
   TH2D *Href_ntrk;
   TH2D *Htof_ref;
   TH2D *Htof_ref_r;
   TH2D *Htof_ntrk;

   TH1D *Hntrk;
   TH1D *Hntrk_C;
   TH1D *Hntrk_P;

   TH1D *Hpt;
   TH1D *Heta_n;
   TH1D *Heta_n_removed;
   TH1D *Heta;
   TH1D *Hphi;
   TH1D *Hcharge;
   TH1D *Hdca;
   TH2D *Hphi_vs_eta;




   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        refmult;
   Double_t        Vz;
   Double_t        Vx;
   Double_t        Vy;
   Double_t        Vz_VPD;
   Double_t        zdc;
   Double_t        bbc;
   Double_t        NBTOFMultfit;
   Int_t           numTrk;
   Double_t        Pt[MAXTRACK];   //[numTrk]
   Double_t        Eta[MAXTRACK];   //[numTrk]
   Double_t        Phi[MAXTRACK];   //[numTrk]
   Double_t        Charge[MAXTRACK];   //[numTrk]
   Double_t        Dca[MAXTRACK];   //[numTrk]
   Int_t           NumEPDHit;
   Double_t        NumMip[MAXTRACK];   //[NumEPDHit]
   Int_t           TileID[MAXTRACK];   //[NumEPDHit]
   Double_t        EPDX[MAXTRACK];   //[NumEPDHit]
   Double_t        EPDY[MAXTRACK];   //[NumEPDHit]
   Double_t        EPD_Psi;
   Double_t        TPC_Psi;

   // List of branches
   TBranch        *b_refmult;   //!
   TBranch        *b_Vz;   //!
   TBranch        *b_Vx;   //!
   TBranch        *b_Vy;   //!
   TBranch        *b_Vz_VPD;   //!
   TBranch        *b_zdc;   //!
   TBranch        *b_bbc;   //!
   TBranch        *b_NBTOFMultfit;   //!
   TBranch        *b_numTrk;   //!
   TBranch        *b_Pt;   //!
   TBranch        *b_Eta;   //!
   TBranch        *b_Phi;   //!
   TBranch        *b_Charge;   //!
   TBranch        *b_Dca;   //!
   TBranch        *b_NumEPDHit;   //!
   TBranch        *b_NumMip;   //!
   TBranch        *b_TileID;   //!
   TBranch        *b_EPDX;   //!
   TBranch        *b_EPDY;   //!
   TBranch        *b_EPD_Psi;   //!
   TBranch        *b_TPC_Psi;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("B74F2EA5ADE63E0429C2ECD3ABDDB0DA_449.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("B74F2EA5ADE63E0429C2ECD3ABDDB0DA_449.root");
      }
      f->GetObject("EventTree",tree);

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

   fChain->SetBranchAddress("refmult", &refmult, &b_refmult);
   fChain->SetBranchAddress("Vz", &Vz, &b_Vz);
   fChain->SetBranchAddress("Vx", &Vx, &b_Vx);
   fChain->SetBranchAddress("Vy", &Vy, &b_Vy);
   fChain->SetBranchAddress("Vz_VPD", &Vz_VPD, &b_Vz_VPD);
   fChain->SetBranchAddress("zdc", &zdc, &b_zdc);
   fChain->SetBranchAddress("bbc", &bbc, &b_bbc);
   fChain->SetBranchAddress("NBTOFMultfit", &NBTOFMultfit, &b_NBTOFMultfit);
   fChain->SetBranchAddress("numTrk", &numTrk, &b_numTrk);
   fChain->SetBranchAddress("Pt", Pt, &b_Pt);
   fChain->SetBranchAddress("Eta", Eta, &b_Eta);
   fChain->SetBranchAddress("Phi", Phi, &b_Phi);
   fChain->SetBranchAddress("Charge", Charge, &b_Charge);
   fChain->SetBranchAddress("Dca", Dca, &b_Dca);
   fChain->SetBranchAddress("NumEPDHit", &NumEPDHit, &b_NumEPDHit);
   fChain->SetBranchAddress("NumMip", NumMip, &b_NumMip);
   fChain->SetBranchAddress("TileID", TileID, &b_TileID);
   fChain->SetBranchAddress("EPDX", EPDX, &b_EPDX);
   fChain->SetBranchAddress("EPDY", EPDY, &b_EPDY);
   fChain->SetBranchAddress("EPD_Psi", &EPD_Psi, &b_EPD_Psi);
   fChain->SetBranchAddress("TPC_Psi", &TPC_Psi, &b_TPC_Psi);
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
