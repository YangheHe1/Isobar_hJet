#define readtree_cxx
#include "readtree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

readtree::readtree(string filelist, int fr, int tr,int job):from(fr),to(tr),Njob(job){

    TChain* chain  = new TChain("METree","");

    char fname[400];
    ifstream lis(filelist.c_str());
    int cnt=0;
    while(!lis.eof()){
        string filename;
        lis >> filename;
        sprintf(fname,"%s",filename.c_str());

        if(!filename.empty()) {
            if(cnt<to && cnt>=from){
                cout << fname << endl; chain->Add(fname);
            }
        }
        cnt++;
        if(cnt>1000000) {cout<<"Too Many Files"<<endl;break;}
    }

    Init(chain);

}

void readtree::InitHist(){

   sprintf(name,"Mixper/MC_Out_from%d_to%d.root", from, to);
   fout = new TFile(name,"RECREATE");

   Hvz = new TH1D("Hvz","vz distribution",160,-40,40);
   Hvz->GetXaxis()->SetTitle("vz");
   Hvz->GetYaxis()->SetTitle("events");

   Hntrk = new TH1D("Hntrk","track distribution for whole central(scaled)",1000,0,1000);
   Hntrk->GetXaxis()->SetTitle("charged track");
   Hntrk->GetYaxis()->SetTitle("events");   

   Hntrk_class = new TH1D("Hntrk_class","track distribution for one class(not scaled)",1000,0,1000);
   Hntrk_class->GetXaxis()->SetTitle("charged track");
   Hntrk_class->GetYaxis()->SetTitle("events"); 

    Hpt = new TH1D("Hpt","Hpt",100,0,50);
    Hpt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    Hpt->GetYaxis()->SetTitle("counts");

    Heta = new TH1D("Heta","Heta",40,-2,2);
    Heta->GetXaxis()->SetTitle("#eta");
    Heta->GetYaxis()->SetTitle("counts");

    Hphi= new TH1D("Hphi","Hphi",80,-4,4);
    Hphi->GetXaxis()->SetTitle("#phi");
    Hphi->GetYaxis()->SetTitle("counts");

    Hphi_vs_eta = new TH2D("Hphi_vs_eta","Hphi_vs_eta",60,-3,3,20,-1,1);
    Hphi_vs_eta->GetXaxis()->SetTitle("#phi");
    Hphi_vs_eta->GetYaxis()->SetTitle("#eta");

    Hcharge = new TH1D("Hcharge","Hcharge ",10,-5,5);
    Hcharge ->GetXaxis()->SetTitle("charge");
    Hcharge ->GetYaxis()->SetTitle("counts");

}

void readtree::TrkLoop(){

   for(int j=0; j < ME_numTrk; j++){
      double _pt = ME_Pt[j];
      double _eta = ME_Eta[j];
      double _phi = ME_Phi[j];
      double _dca = ME_Dca[j];
      double _charge = ME_Charge[j];

      Hpt->Fill(_pt);
      if(_pt<0.5){
      Hphi->Fill(_phi);
      Heta->Fill(_eta);
      Hcharge->Fill(_charge);
      Hphi_vs_eta->Fill(_phi,_eta);}

   }

}


void readtree::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L readtree.C
//      Root > readtree t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   cout<<"total events "<<nentries<<endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Hntrk->Fill(ME_numTrk);
      Hntrk_class->Fill(ME_numTrk);
      Hvz->Fill(ME_Vz);

      TrkLoop();

   }
   double Nev = Hntrk->GetEntries();
   Hntrk->Scale(1./Nev);
   double scaler = 1./Njob;
   Hntrk->Scale(scaler);
   fout -> Write();

}
