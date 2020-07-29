#define readtree_cxx
#include "readtree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

readtree::readtree(string filelist, int fr, int tr):from(fr),to(tr){

    TChain* chain  = new TChain("EventTree","");

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

int readtree::centrality( int _refMult ){
    
    float   CentralityBins  [NCENT] = {213,149,102,68,44,27,15,8,4} ;  //
    Int_t   MiddleBinID     [NCENT] = { 0,1,2,3,4,5,6,7,8};  // ID Number
    
    Int_t   myCentrality=-1;
    
    for(int i=0;i!=NCENT;i++){
        if( _refMult > CentralityBins[i] ){
            myCentrality = MiddleBinID[i];
            break;
        }
        else{
            myCentrality = -1;
        }
    }
    
    return myCentrality ;
    
}

bool readtree::refmult_check (int __nBTOFMatch, int __refMult){

//We do the Y-projection in a nBTOFMatch window, fit by double negative binomial distribution then get these parameters.
///Recommand use max as 3, set min as 4.

double a0=6.44225, a1=1.37398, a2=-0.00446098, a3=2.60578e-05, a4= -6.84022e-08, a5=6.18087e-11;
//double b0=2.52126730672253, b1=0.128066911940844, b2=-0.000538959206681944, b3=1.21531743671716e-06, b4=-1.01886685404478e-09;
//double c0=4.79427731664144, c1=0.187601372159186, c2=-0.000849856673886957, c3=1.9359155975421e-06, c4=-1.61214724626684e-09;



double refmultcutmax=a0+a1*(__nBTOFMatch)+a2*pow(__nBTOFMatch,2)+a3*pow(__nBTOFMatch,3)+a4*pow(__nBTOFMatch,4)+a5*pow(__nBTOFMatch,5);


        //cout<<"refmult,maxcut,mincut= "<<__refMult<<" "<<refmultcutmax<<" "<<refmultcutmin<<endl;

        if( __refMult<refmultcutmax){
                return true;
        }else{
                return false;
        }

}

bool readtree::ntrk_check(int __nBTOFMatch, int __ntrk){

    double b0=-14.7884, b1=1.07711, b2=-0.000186869, b3=1.20443e-06, b4=-1.745e-09;
    double c0=11.1393, c1=1.88776, c2=-0.0024491, c3=8.43354e-06, c4=-1.29845e-08;

    double ntrkcutmin=b0+b1*(__nBTOFMatch)+b2*pow(__nBTOFMatch,2)+b3*pow(__nBTOFMatch,3)+b4*pow(__nBTOFMatch,4);
    double ntrkcutmax=c0+c1*(__nBTOFMatch)+c2*pow(__nBTOFMatch,2)+c3*pow(__nBTOFMatch,3)+c4*pow(__nBTOFMatch,4);

    if( __ntrk<ntrkcutmax && __ntrk>ntrkcutmin){
                return true;
        }else{
                return false;
        }



}


void readtree::InitHist(){

   sprintf(name,"result/MC_Out_from%d_to%d.root", from, to);
   fout = new TFile(name,"RECREATE");

   Hvz = new TH1D("Hvz","vz distribution",160,-40,40);
   Hvz->GetXaxis()->SetTitle("vz");
   Hvz->GetYaxis()->SetTitle("events");

  Heta_n = new TH1D("Heta_n","Heta_n",200,0,2);
  Heta_n_removed = new TH1D("Heta_n_removed","Heta_n_removed",20,0,2);

   Hntrk = new TH1D("Hntrk","track distribution",1000,0,1000);
   Hntrk->GetXaxis()->SetTitle("charged track");
   Hntrk->GetYaxis()->SetTitle("events");

   Hntrk_C = new TH1D("Hntrk_C","track distribution for whole central",1000,0,1000);
   Hntrk_C->GetXaxis()->SetTitle("charged track");
   Hntrk_C->GetYaxis()->SetTitle("events");

   Hntrk_P = new TH1D("Hntrk_P","track distribution for whole peripheral",1000,0,1000);
   Hntrk_P->GetXaxis()->SetTitle("charged track");
   Hntrk_P->GetYaxis()->SetTitle("events");

   Hrefmult = new TH1D ("Hrefmult","Hrefmult",800,0,800);
   Hrefmult->GetXaxis()->SetTitle("refmult");
   Hrefmult->GetYaxis()->SetTitle("events");

   Href_ntrk = new TH2D("Href_ntrk","Href_ntrk",500,0,500,1000,0,1000);
   Href_ntrk->GetXaxis()->SetTitle("refMult");
   Href_ntrk->GetYaxis()->SetTitle("charged track");

   Htof_ref = new TH2D("Htof_ref","Htof_ref",500,0,500,1000,0,1000);
   Htof_ref->GetXaxis()->SetTitle("refMult");
   Htof_ref->GetYaxis()->SetTitle("TOFMatch");

   Htof_ref_r = new TH2D("Htof_ref_r","Htof_ref_r",1000,0,1000,1000,0,1000);
   Htof_ref_r->GetXaxis()->SetTitle("TOFMatch");
   Htof_ref_r->GetYaxis()->SetTitle("refMult");

   Htof_ntrk = new TH2D("Htof_ntrk","Htof_ntrk",500,0,500,1000,0,1000);
   Htof_ntrk->GetYaxis()->SetTitle("charged track");
   Htof_ntrk->GetXaxis()->SetTitle("TOFMatch");

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

    Hdca = new TH1D("Hdca","Hdca",150,0,1.5);
    Hdca ->GetXaxis()->SetTitle("dca");
    Hdca ->GetYaxis()->SetTitle("counts");

}

void readtree::TrkLoop(){


    double Neta=0;
    double Peta=0;
   for(int j=0; j < numTrk; j++){
      double _pt = Pt[j];
      double _eta = Eta[j];
      double _phi = Phi[j];
      double _dca = Dca[j];
      double _charge = Charge[j];

      if(TMath::Abs(_eta)>1) continue;
      if(TMath::Abs(_dca)>1) continue;
      if(_pt<0.2||_pt>30)  continue;

      NPt++;
      if(_eta<0) Neta++;
      if(_eta>0) Peta++;
      Hpt->Fill(_pt);
      
      if(_pt<30){
      Hphi->Fill(_phi);
      Heta->Fill(_eta);
      Hcharge->Fill(_charge);
      Hdca->Fill(_dca);
      Hphi_vs_eta->Fill(_phi,_eta);}

   }
        double ngap=abs(Neta-Peta)/NPt;
        if(refmult<27&&refmult>8) Heta_n->Fill(ngap);

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

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //____cut______________
      if(Vz>25) continue;
      //tof vs ref pile up rejection________________________________________
      //linear
      
      int Ntofmatch = NBTOFMultfit;
     /* double NBTOFMult_fit;
      if(Ntofmatch < 14){ NBTOFMult_fit = 0.59207 + 2.1317*Ntofmatch; }
      else{ NBTOFMult_fit = 18.770 + 1.0699*Ntofmatch; }
      if(refmult >= NBTOFMult_fit) continue;
      */
      //fit
      if(!refmult_check(Ntofmatch,refmult)) continue;

      //_____________________________________
	
	Htof_ref_r->Fill(Ntofmatch,refmult);     
      //___ntrk vs tof pile up_______________
      int NTRK=0;
      for(int j=0; j < numTrk; j++){
      double _pt = Pt[j];
      double _eta = Eta[j];
      double _dca = Dca[j];
      
      if(TMath::Abs(_eta)>1) continue;
      if(TMath::Abs(_dca)>1) continue;
      if(_pt<0.2||_pt>30)  continue;

      NTRK++;

      }

      if(!ntrk_check(Ntofmatch,NTRK)) continue;

      //___________________________________________________


      NPt=0;
      TrkLoop();

      Hntrk->Fill(NPt);
      int centid= centrality(refmult);
      if(centid==0) Hntrk_C->Fill(NPt);
      if(centid==6||centid==7) Hntrk_P->Fill(NPt);
      Hrefmult->Fill(refmult);
      Href_ntrk->Fill(refmult,NPt);
      Htof_ref->Fill(refmult,Ntofmatch);
      Htof_ntrk->Fill(Ntofmatch,NPt);
      Hvz->Fill(Vz);

      


   }

   fout -> Write();

}
