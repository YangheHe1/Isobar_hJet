//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 25 10:17:12 2020 by ROOT version 5.34/38
// from TTree JetTree/JetTree
// found on file: SE_JetTree_619.root
//////////////////////////////////////////////////////////

#ifndef readjet_h
#define readjet_h

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
#include "TProfile.h"
#include <string>
#include <TH2.h>
#include <TH1.h>
#include <TRandom.h>
#include <fstream>
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"

using std::vector;
using std::string;

const Int_t MAXTRACK=10000;
const Int_t NCENT=9;
static const double value_pi = 3.14;
// Fixed size dimensions of array or collections stored in the TTree if any.

class readjet {
public :

   readjet(string filelist, int fr, int tr, int njob, double area_cut);
   void InitHist();
   void JetLoop();
   int centrality( int _refMult );

   int from;
   int to;
   int Njob;
   double AreaCut;
   char name[200];
   TFile *fout;
   int centid;
   //TRandom phi11; 
   //double TriggerPhi;
   TH1D *CjetPt[2];
   TH1D *CjetPt_9_30[2];
   TH1D *CjetPt_3_7[2];
   TH1D *CjetPt_7_9[2];
   TH1D *Ntrigger[2];
   TH1D *Ntrigger9_30[2];
   TH1D *Ntrigger3_7[2];
   TH1D *Ntrigger7_9[2];
   TH1D *triggerPt[2];
   TH1D *Hrho[2];

   TH1D *HArea[2];
   TH1D *HArea9_30[2];
   TH2D *HArea_Pt[2];
   TH2D *HArea_Pt9_30[2];
   
   TH1D *HArea_JPt5[2];
   TH1D *HArea9_30_JPt5[2];
   TH2D *HArea_Pt_JPt5[2];
   TH2D *HArea_Pt9_30_JPt5[2];

   TH2D *Hrho_vs_M[2];
   TH1D *Harea7_30;
   TH1D *Harea9_30;
   TH2D *Harea_Pt;
   TH2D *Harea_Pt9_30;
   TH1D *Harea7_30_Pt5;
   TH1D *Harea9_30_Pt5;
   TH1D *HTrg7_30;
   TH1D *HTrg9_30;
   TH1D *Hrho_all;
   TH1D *Hrho7_10;
   TH1D *Hrho10_20;
   TH1D *Hrho20_30;

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           eventIndex;
   Double_t        refmult;
   Int_t           PrimTrk;
   Double_t        TrgEta;
   Double_t        TrgPhi;
   Double_t        TrgEt;
   Int_t           TrgId;
   Double_t        Rho;
   Double_t        Sigma;
   Int_t           NJets;
   Int_t           JetIndex[MAXTRACK];   //[NJets]
   Int_t           JetNCons[MAXTRACK];   //[NJets]
   Double_t        JetPt[MAXTRACK];   //[NJets]
   Double_t        JetPtCorr[MAXTRACK];   //[NJets]
   Double_t        JetEta[MAXTRACK];   //[NJets]
   Double_t        JetPhi[MAXTRACK];   //[NJets]
   Double_t        JetE[MAXTRACK];   //[NJets]
   Double_t        JetArea[MAXTRACK];   //[NJets]

   // List of branches
   TBranch        *b_EventIndex;   //!
   TBranch        *b_refmult;   //!
   TBranch        *b_PrimTrk;   //!
   TBranch        *b_TrgEta;   //!
   TBranch        *b_TrgPhi;   //!
   TBranch        *b_TrgEt;   //!
   TBranch        *b_TrgId;   //!
   TBranch        *b_Rho;   //!
   TBranch        *b_Sigma;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_JetIndex;   //!
   TBranch        *b_JetNCons;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetPtCorr;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetE;   //!
   TBranch        *b_JetArea;   //!

   readjet(TTree *tree=0);
   virtual ~readjet();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef readjet_cxx
readjet::readjet(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SE_JetTree_619.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("SE_JetTree_619.root");
      }
      f->GetObject("JetTree",tree);

   }
   Init(tree);
}

readjet::~readjet()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t readjet::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t readjet::LoadTree(Long64_t entry)
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

void readjet::Init(TTree *tree)
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

   fChain->SetBranchAddress("eventIndex", &eventIndex, &b_EventIndex);
   fChain->SetBranchAddress("refmult", &refmult, &b_refmult);
   fChain->SetBranchAddress("PrimTrk", &PrimTrk, &b_PrimTrk);
   fChain->SetBranchAddress("TrgEta", &TrgEta, &b_TrgEta);
   fChain->SetBranchAddress("TrgPhi", &TrgPhi, &b_TrgPhi);
   fChain->SetBranchAddress("TrgEt", &TrgEt, &b_TrgEt);
   fChain->SetBranchAddress("TrgId", &TrgId, &b_TrgId);
   fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
   fChain->SetBranchAddress("Sigma", &Sigma, &b_Sigma);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("JetIndex", JetIndex, &b_JetIndex);
   fChain->SetBranchAddress("JetNCons", JetNCons, &b_JetNCons);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetPtCorr", JetPtCorr, &b_JetPtCorr);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetE", JetE, &b_JetE);
   fChain->SetBranchAddress("JetArea", JetArea, &b_JetArea);
   Notify();
}

Bool_t readjet::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void readjet::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t readjet::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef readjet_cxx
