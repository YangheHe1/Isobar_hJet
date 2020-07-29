#define readtree_cxx
#include "readtree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

readtree::readtree(string filelist, int fr, int tr, int loop):from(fr),to(tr){

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

   myloop     = loop;

}

int readtree::centrality( int _refMult ){
    
    float   CentralityBins  [NCENT] = {235,180,137,103,76,54,24} ;  //
    Int_t   MiddleBinID     [NCENT] = { 0,1,2,3,4,5,6};  // ID Number
    
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


void readtree::InitHist(){
    if(myloop == 0) sprintf(name,"recenter/recenter_loop%d_fr%d.root",myloop,from);
    if(myloop == 1) sprintf(name,"flattening/flattening_loop%d_fr%d.root",myloop,from);
    if(myloop == 2) sprintf(name,"result/result_loop%d_fr%d.root",myloop,from);
    fout = new TFile(name,"recreate");
    
    //===============recentering part=================
    
    for (int ih=0; ih<NHAR; ih++) {
        
      for (int icent=0; icent<NCENT; icent++) {
         if (myloop !=0)  continue;
         sprintf(name,"hqxqy_raw_har%d_cent%d", ih, icent);
                hqxqy_raw[ih][icent] = new TH2F(name,"",400,-1.,1.,400,-1.,1.);
      }
        
    }
    
   
    
    for (int ih=0; ih<NHAR; ih++) {
        
      if (myloop !=0)  continue;
      sprintf(name,"h_profTrk_sumx_har%d", ih);
      h_profTrk_sumx[ih] = new TProfile(name,"",100,0-0.5, 100-0.5);
      h_profTrk_sumx[ih] ->Sumw2();
            
      sprintf(name,"h_profTrk_sumy_har%d", ih);
      h_profTrk_sumy[ih] = new TProfile(name,"",100,0-0.5, 100-0.5);
      h_profTrk_sumy[ih] ->Sumw2();
        
    }
    //end of recentering

    Hepd_psi = new TH1D("Hepd_psi","Hepd_psi",32,-(0.5-1.0/64)*PI,(1.5-1.0/64)*PI);
    Hepd_psi->GetXaxis()->SetTitle("#Psi");
    Hepd_psi->GetYaxis()->SetTitle("events");

    Htpc_psi = new TH1D("Htpc_psi","Htpc_psi",32,-(0.5-1.0/64)*PI,(1.5-1.0/64)*PI);
    Htpc_psi->GetXaxis()->SetTitle("#Psi");
    Htpc_psi->GetYaxis()->SetTitle("events");
    
    //Flattening part
    
    for (int ih=0; ih<NHAR; ih++) {
            for (int icent=0 ; icent<NCENT; icent++) {
                
                if (myloop !=1 )  continue;
                sprintf(name,"hProf_Trkflatc_har%d_cent%d", ih, icent);
                h_prof_TrkflatC[ih][icent] = new TProfile(name, "", NK, -0.5, NK-0.5, -1.1, 1.1);
                h_prof_TrkflatC[ih][icent] ->Sumw2();
                
                sprintf(name,"hProf_Trkflats_har%d_cent%d", ih, icent);
                h_prof_TrkflatS[ih][icent] = new TProfile(name, "", NK, -0.5, NK-0.5, -1.1, 1.1);
                h_prof_TrkflatS[ih][icent] ->Sumw2();
                
            }
    }//end of ic
    //end of Flattening

    for (int ih=0; ih<NHAR; ih++) {
         for (int icent=0; icent<NCENT; icent++) {
               if (myloop !=2)  continue;
               sprintf(name,"hqxqy_har%d_cent%d", ih, icent);
               hqxqy[ih][icent] = new TH2F(name,"",400,-1.,1.,400,-1.,1.);
         }
    }
    
    for (int itype=0; itype<3; itype++) {
        for (int ih=0; ih<NHAR; ih++) {
                for (int icent=0 ; icent<NCENT; icent++) {
                    if (myloop !=2)  continue;
                    sprintf(name,"hqn_type%d_har%d_cent%d", itype, ih, icent);
                    hqn[itype][ih][icent] = new TH1D(name, name, 1000, 0, 1);
                }
        }
    }
    
    for (int itype=0; itype<3; itype++) {
        for (int ih=0; ih<NHAR; ih++) {
         
                for (int icent=0 ; icent<NCENT; icent++) {
                    if (myloop !=2)  continue;
                    sprintf(name,"hPsi_type%d_har%d_cent%d", itype, ih, icent);
                    hPsi[itype][ih][icent] = new TH1D(name, name, 360, -PI, PI);
                }
            
        }
    }


}

void readtree::TrkLoop(){


   memset(QxRaw, 0, sizeof(QxRaw));
   memset(QyRaw, 0, sizeof(QyRaw));
   memset(QwRaw, 0, sizeof(QwRaw));

   for(int j=0; j < numTrk; j++){
      double _pt = Pt[j];
      double _eta = Eta[j];
      double _phi = Phi[j];
      double _dca = Dca[j];
      double _charge = Charge[j];

      if(_pt>5.0) continue;
      if(TMath::Abs(_eta)>1.0) continue;
      if(TMath::Abs(_dca)>1.0) continue;
      double weight = (_pt<2.0)?_pt:2.0;

      for(int ih = 0; ih != 2; ih++){
            
            
                
        QwRaw[ih]+= weight;
        QxRaw[ih] += weight*cos((ih+2)*_phi);
        QyRaw[ih] += weight*sin((ih+2)*_phi);
                
      }//ih

   }   // track loop

   double Psi2 = TMath::ATan2(QyRaw[0],QxRaw[0]);
   Psi2 /= 2.0;

   Htpc_psi->Fill(Psi2);
   Hepd_psi->Fill(EPD_Psi);


}

bool readtree::ForRecenter(){
    
    if(myloop!=0) return true;
    
    for (int ih=0; ih != NHAR; ih++) {
   
            
      hqxqy_raw[ih][CentID] ->Fill(QxRaw[ih]/QwRaw[ih],QyRaw[ih]/QwRaw[ih]);
            
      h_profTrk_sumx[ih] ->Fill(CentID, QxRaw[ih]/QwRaw[ih]);
      h_profTrk_sumy[ih] ->Fill(CentID, QyRaw[ih]/QwRaw[ih]);
            
    }
    
    return true;
    
}

bool readtree::ForFlattening(){
    
    if (myloop!=1) return true;
    
    static int initArray = 0;
    
    if(initArray==0){
        
        memset(Trk_sumxN, 0, sizeof(Trk_sumxN));
        memset(Trk_sumyN, 0, sizeof(Trk_sumyN));
        
        
        TFile* myfin;
        sprintf(name,"recenter/recenter_loop0_fr%d.root",from);
        myfin = TFile::Open(name);
        cout<<"We read in recenter info "<<name<<endl;
        
        for (int ih=0; ih != NHAR; ih++) {
            
                
         sprintf(name,"h_profTrk_sumx_har%d", ih);
         h_profTrk_sumxN[ih] = (TProfile*)myfin->Get(name);
         sprintf(name,"h_profTrk_sumy_har%d", ih);
         h_profTrk_sumyN[ih] = (TProfile*)myfin->Get(name);
                
            
        }
        
        cout<<"We extract the calibTrk information"<<endl;
        for (int ih=0; ih<NHAR; ih++) {
            
            for (int icent=0; icent<NCENT; icent++) {
               Trk_sumxN[ih][icent] = h_profTrk_sumxN[ih]->GetBinContent(icent+1) ;
               Trk_sumyN[ih][icent] = h_profTrk_sumyN[ih]->GetBinContent(icent+1) ;
                    
            }
            
        }
        
        cout<<"finish readin"<<endl;
        initArray = 1;
    }//end of read in
    
    //    cout<<"finish readin"<<endl;
    
    for (int ih=0; ih<NHAR; ih++) {
        
            
      double qx = QxRaw[ih]/QwRaw[ih];
      double qy = QyRaw[ih]/QwRaw[ih];
            
      double QxC  = qx - Trk_sumxN[ih][CentID];
      double QyC  = qy - Trk_sumyN[ih][CentID];
            
      double PsiC = atan2(QyC, QxC)/double(ih+2);
      double psi2pi = PsiC*(ih+2);

      for (int ik=0; ik<NK; ik++) {
         h_prof_TrkflatC[ih][CentID]->Fill(ik, cos( (ik+1)*psi2pi) );
         h_prof_TrkflatS[ih][CentID]->Fill(ik, sin( (ik+1)*psi2pi) );
      }
            
      
    }
    
    return true;
}

bool readtree::CalQn(){
    
    if (myloop!=2) return true;
    
    memset(QxNew, 0, sizeof(QxNew));
    memset(QyNew, 0, sizeof(QyNew));
    memset(QwNew, 0, sizeof(QwNew));
    memset(PsiNew, 0, sizeof(PsiNew));
    
    static int initArrayx = 0;
    
    if(initArrayx==0){
        cout<<"We begin extraction in calQn"<<endl;
        TFile* myfin;
        sprintf(name,"recenter/recenter_loop0_fr%d.root",from);
        myfin = TFile::Open(name);
        cout<<"We read in recenter info "<<name<<endl;
        
        for (int ih=0; ih != NHAR; ih++) {
         
         sprintf(name,"h_profTrk_sumx_har%d", ih);
         h_profTrk_sumxN[ih] = (TProfile*)myfin->Get(name);
         sprintf(name,"h_profTrk_sumy_har%d", ih);
         h_profTrk_sumyN[ih] = (TProfile*)myfin->Get(name);
                
   
        }
        
        TFile* myfin2;
        sprintf(name,"flattening/flattening_loop1_fr%d.root",from);
        myfin2 = TFile::Open(name);
        cout<<"We read in flattening info "<<name<<endl;
        
        //Flattening part
        for (int ih=0; ih<NHAR; ih++) {
            
         for (int icent=0 ; icent<NCENT; icent++) {
                    
            sprintf(name,"hProf_Trkflatc_har%d_cent%d", ih, icent);
            h_prof_TrkflatCN[ih][icent] = (TProfile*)myfin2->Get(name);
            sprintf(name,"hProf_Trkflats_har%d_cent%d", ih, icent);
            h_prof_TrkflatSN[ih][icent] = (TProfile*)myfin2->Get(name);
    
         }
            
        }//end of ic
        
        initArrayx=1;
    }//end of read in
    
    for (int ih=0; ih<NHAR; ih++) {
        
            
            double Qx=QxRaw[ih];
            double Qy=QyRaw[ih];
            double Qw=QwRaw[ih];
            
            double  QxC = Qx -  Qw*h_profTrk_sumxN[ih]->GetBinContent(CentID+1) ;
            double  QyC = Qy -  Qw*h_profTrk_sumyN[ih]->GetBinContent(CentID+1) ;
            
            double Psi  = atan2(Qy, Qx)/double(ih+2);
            double PsiC = atan2(QyC, QxC)/double(ih+2);
            
            double psiC2pi = (ih+2)*PsiC;
            double deltaPsi = 0;
            
            for (int ik=0; ik<NK; ik++) {
                
                double flatcos = 0;
                double flatsin = 0;
                
                flatcos = h_prof_TrkflatCN[ih][CentID]->GetBinContent(ik+1);
                flatsin = h_prof_TrkflatSN[ih][CentID]->GetBinContent(ik+1);
                double cosPsiC = cos( (ik+1)*psiC2pi );
                double sinPsiC = sin( (ik+1)*psiC2pi );
                deltaPsi += (-flatsin*cosPsiC + flatcos*sinPsiC)*2/(ik+1);
                
            }
            
            double PsiF = atan2( sin(psiC2pi+deltaPsi), cos(psiC2pi+deltaPsi) )/double(ih+2);
            //Keep Magnitude, just shift the phase
            double QxF  = sqrt(QxC*QxC + QyC*QyC) * cos( (ih+2)*PsiF );
            double QyF  = sqrt(QxC*QxC + QyC*QyC) * sin( (ih+2)*PsiF );
            
            QxNew[0][ih] = Qx;
            QyNew[0][ih] = Qy;
            
            QxNew[1][ih] = QxC;
            QyNew[1][ih] = QyC;
            
            QxNew[2][ih] = QxF;
            QyNew[2][ih] = QyF;
            
        
    }//ih
    
    for (int ih=0; ih<NHAR; ih++) {
            
      hqxqy[ih][CentID] -> Fill(QxNew[2][ih]/QwRaw[ih],QyNew[2][ih]/QwRaw[ih]);
            
    }
    
    for (int itype=0; itype<3; itype++) {
        for (int ih=0; ih<NHAR; ih++) {

                
         PsiNew[itype][ih] = atan2(QyNew[itype][ih], QxNew[itype][ih])/(ih+2);
         double qn = sqrt( pow(QyNew[itype][ih],2) +pow(QxNew[itype][ih],2)  )/QwRaw[ih];
                
         hqn[itype][ih][CentID] ->Fill(qn);
         hPsi[itype][ih][CentID]->Fill(PsiNew[itype][ih]);
      }
    }//calib
    
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
   cout<<"We have "<<nentries<<" events"<<endl;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%10000==0) cout<<"working on "<<jentry<<endl;
       
       CentID     = centrality(refmult);
       if(CentID<0 || CentID > NCENT ) continue;
       TrkLoop();
       ForRecenter();
       ForFlattening();
       CalQn();
   }

   fout -> Write();

}
