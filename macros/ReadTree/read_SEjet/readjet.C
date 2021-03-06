#define readjet_cxx
#include "readjet.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

readjet::readjet(string filelist, int fr, int tr, int njob, double area_cut):from(fr),to(tr),Njob(njob),AreaCut(area_cut){

    TChain* chain  = new TChain("JetTree","");

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

int readjet::centrality( int _refMult ){
    
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

void readjet::InitHist(){

sprintf(name,"result/Out_from%d_to%d.root", from, to);
fout = new TFile(name,"RECREATE");

for(int icent=0;icent<2;icent++){

sprintf(name,"recoil_jet_cent%d",icent);
CjetPt[icent] = new TH1D(name,"",70,-20,50);
CjetPt[icent] ->Sumw2();

sprintf(name,"recoil_jet_trg9_30_cent%d",icent);
CjetPt_9_30[icent] = new TH1D(name,"",70,-20,50);
CjetPt_9_30[icent] ->Sumw2();

sprintf(name,"recoil_jet_trg3_7_cent%d",icent);
CjetPt_3_7[icent] = new TH1D(name,"",70,-20,50);
CjetPt_3_7[icent] ->Sumw2();

sprintf(name,"recoil_jet_trg7_9_cent%d",icent);
CjetPt_7_9[icent] = new TH1D(name,"",70,-20,50);
CjetPt_7_9[icent] ->Sumw2();

sprintf(name,"number_of_trigger_cent%d",icent);
Ntrigger[icent] = new TH1D(name,"",2,0,2);

sprintf(name,"number_of_trigger9_30_cent%d",icent);
Ntrigger9_30[icent] = new TH1D(name,"",2,0,2);

sprintf(name,"number_of_trigger3_7_cent%d",icent);
Ntrigger3_7[icent] = new TH1D(name,"",2,0,2);

sprintf(name,"number_of_trigger7_9_cent%d",icent);
Ntrigger7_9[icent] = new TH1D(name,"",2,0,2);

sprintf(name,"trigger_pt_dist_cent%d",icent);
triggerPt[icent] = new TH1D(name,"",40,0,40);

sprintf(name,"Rho_cent%d",icent);
Hrho[icent] = new TH1D(name,"rho",100,0,50);

sprintf(name,"Rho_vs_M_cent%d",icent);
Hrho_vs_M[icent] = new TH2D(name,"",820,0,820,100,0,50);

sprintf(name,"Area_trg7_30_cent%d",icent);
HArea[icent] = new TH1D(name,"Area",150,0,1.5);

sprintf(name,"Area_trg9_30_cent%d",icent);
HArea9_30[icent] = new TH1D(name,"Area 9-30",150,0,1.5);

sprintf(name,"Area_vs_Pt_trg7_30_cent%d",icent);
HArea_Pt[icent] = new TH2D(name,"Area vs pt",150,0,1.5,70,-20,50);

sprintf(name,"Area_vs_Pt_trg9_30_cent%d",icent);
HArea_Pt9_30[icent] = new TH2D(name,"Area vs pt trg 9-30",150,0,1.5,70,-20,50);

//JET PT>5
sprintf(name,"Area_trg7_30_JPt5_cent%d",icent);
HArea_JPt5[icent] = new TH1D(name,"Area jet pt>5",150,0,1.5);

sprintf(name,"Area_trg9_30_JPt5_cent%d",icent);
HArea9_30_JPt5[icent] = new TH1D(name,"Area 9-30 pt>5",150,0,1.5);

sprintf(name,"Area_vs_Pt_trg7_30_JPt5_cent%d",icent);
HArea_Pt_JPt5[icent] = new TH2D(name,"Area vs pt jet pt>5",150,0,1.5,70,-20,50);

sprintf(name,"Area_vs_Pt_trg9_30_JPt5_cent%d",icent);
HArea_Pt9_30_JPt5[icent] = new TH2D(name,"Area vs pt trg 9-30 jet pt>5",150,0,1.5,70,-20,50);


}
Hrho_all = new TH1D("Hrho_all","Hrho_all",100,0,50);

Harea7_30 = new TH1D("Harea7_30","Harea trg7-30",150,0,1.5);
Harea7_30_Pt5 = new TH1D("Harea7_30_Pt5","Harea trg7-30 jet pt>5",150,0,1.5);
Harea9_30 = new TH1D("Harea9_30","Harea trg9-30",150,0,1.5);
Harea9_30_Pt5 = new TH1D("Harea9_30_Pt5","Harea trg9-30 jet pt>5",150,0,1.5);
Harea_Pt = new TH2D("Harea_Pt","Area vs pt trg 7-30",150,0,1.5,70,-20,50);
Harea_Pt9_30 = new TH2D("Harea_Pt9_30","Area vs pt trg 9-30",150,0,1.5,70,-20,50);

HTrg7_30 = new TH1D ("Htrg7_30","Htrg7_30",2,0,2);
HTrg9_30 = new TH1D ("Htrg9_30","Htrg9_30",2,0,2);

Hrho7_10 = new TH1D("Hrho7_10","rho TT 7-10",100,0,50);
Hrho10_20 = new TH1D("Hrho10_20","rho TT 10-20",100,0,50);
Hrho20_30 = new TH1D("Hrho20_30","rho TT 20-30",100,0,50);
}

void readjet::JetLoop(){



        for(int j=0; j < NJets; j++)
        {
                double iJetEta = JetEta[j];
                double iJetPt = JetPt[j];
                double iJetPtCorr = JetPtCorr[j];
                double iJetPhi = JetPhi[j];
                double iJetArea = JetArea[j];
                double trig_jet_deltaPhi=0.;
                double deltaPhi= TrgPhi-iJetPhi;
                if( deltaPhi < 0. ){trig_jet_deltaPhi= deltaPhi + 2.*value_pi;}
                else {trig_jet_deltaPhi= deltaPhi;}

                if (trig_jet_deltaPhi < value_pi-(value_pi/4.0)  || trig_jet_deltaPhi > value_pi+(value_pi/4.0) ) continue;

                if(iJetArea<AreaCut) continue;
                double iJetPtC = iJetPt-(iJetArea*Rho);
                //Harea->Fill(iJetArea);

      if(TrgEt>7&&TrgEt<30){

         Harea7_30->Fill(iJetArea);
         Harea_Pt->Fill(iJetArea,iJetPtC);
         if(iJetPtC>5) Harea7_30_Pt5->Fill(iJetArea);

		   if(centid==0){
			   CjetPt[0]->Fill(iJetPtC);

            HArea[0]->Fill(iJetArea);
            HArea_Pt[0]->Fill(iJetArea,iJetPtC);

            if(iJetPtC>5) {
               HArea_JPt5[0]->Fill(iJetArea);
               HArea_Pt_JPt5[0]->Fill(iJetArea,iJetPtC);
            }

			}
         //______________________________________________________
		   if(centid==6||centid==7){
			   CjetPt[1]->Fill(iJetPtC);

            HArea[1]->Fill(iJetArea);
            HArea_Pt[1]->Fill(iJetArea,iJetPtC);

            if(iJetPtC>5) {
               HArea_JPt5[1]->Fill(iJetArea);
               HArea_Pt_JPt5[1]->Fill(iJetArea,iJetPtC);
            }

			}//_________________

      }
      //=================================================================

      if(TrgEt>9&&TrgEt<30){

         Harea9_30->Fill(iJetArea);
         Harea_Pt9_30->Fill(iJetArea,iJetPtC);
         if(iJetPtC>5) Harea9_30_Pt5->Fill(iJetArea);

		   if(centid==0){

			   CjetPt_9_30[0]->Fill(iJetPtC);

            HArea9_30[0]->Fill(iJetArea);
            HArea_Pt9_30[0]->Fill(iJetArea,iJetPtC);

            if(iJetPtC>5) {
               HArea9_30_JPt5[0]->Fill(iJetArea);
               HArea_Pt9_30_JPt5[0]->Fill(iJetArea,iJetPtC);
            }

			}

		   if(centid==6||centid==7){
			   CjetPt_9_30[1]->Fill(iJetPtC);

            HArea9_30[1]->Fill(iJetArea);
            HArea_Pt9_30[1]->Fill(iJetArea,iJetPtC);

            if(iJetPtC>5) {
               HArea9_30_JPt5[1]->Fill(iJetArea);
               HArea_Pt9_30_JPt5[1]->Fill(iJetArea,iJetPtC);
            }

			}


      }

      //==================================================================================
/*      if(TrgEt>3&&TrgEt<7){
         
		   if(centid==0){
			   CjetPt_3_7[0]->Fill(iJetPtC);
			}

		   if(centid==6||centid==7){
			   CjetPt_3_7[1]->Fill(iJetPtC);
			}
      }
*/
      if(TrgEt>7&&TrgEt<9){
		   if(centid==0){
			   CjetPt_7_9[0]->Fill(iJetPtC);
			}

		   if(centid==6||centid==7){
			   CjetPt_7_9[1]->Fill(iJetPtC);
			}
      }



	}

}

void readjet::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L readjet.C
//      Root > readjet t
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
   cout<<"We have "<<nentries<<" events"<<endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      centid= centrality(refmult);
/*
      if(TrgEt>3&&TrgEt<7) {
         if(centid==0)Ntrigger3_7[0]->Fill(1);
         if(centid==6||centid==7) Ntrigger3_7[1]->Fill(1);

      }
*/      
      if(TrgEt>7&&TrgEt<30){
         if(centid==0) {
            Ntrigger[0]->Fill(1);
            if(TrgEt>9&&TrgEt<30) Ntrigger9_30[0]->Fill(1);
            
            if(TrgEt>7&&TrgEt<9) Ntrigger7_9[0]->Fill(1);
            Hrho[0]->Fill(Rho);
            triggerPt[0]->Fill(TrgEt);
            Hrho_vs_M[0]->Fill(PrimTrk,Rho);}
         if(centid==6||centid==7) {
            Ntrigger[1]->Fill(1);
            if(TrgEt>9&&TrgEt<30) Ntrigger9_30[1]->Fill(1);
            
            if(TrgEt>7&&TrgEt<9) Ntrigger7_9[1]->Fill(1);
            Hrho[1]->Fill(Rho);
            triggerPt[1]->Fill(TrgEt);
            Hrho_vs_M[1]->Fill(PrimTrk,Rho);}
      	   
	      Hrho_all->Fill(Rho);   
         HTrg7_30->Fill(1);
         if(TrgEt>9&&TrgEt<30) HTrg9_30->Fill(1);
      }

	if(TrgEt>7&&TrgEt<10)Hrho7_10->Fill(Rho);
	if(TrgEt>10&&TrgEt<20)Hrho10_20->Fill(Rho);
	if(TrgEt>20&&TrgEt<30)Hrho20_30->Fill(Rho);
      JetLoop();

   }
/*
   double scaler=1./Njob;
   for(int i=0;i<2;i++){
	double NC= Hrho[i]->GetEntries();
	if(NC!=0)Hrho[i]->Scale(1./NC);
        Hrho[i]->Scale(scaler);
   }
*/
   fout -> Write();
}
