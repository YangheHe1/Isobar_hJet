#include<iostream>
#include<vector>
#include<algorithm>
#include<TH1.h>
#include<fstream>
#include"TFile.h"

void Divide_Nch(){

  TFile *fin=new TFile("Step2/QA_afterallcuts_Sep10.root","read");

  TH1D *hNch_C=(TH1D*)fin->Get("HNCharge_C_1");
  Int_t nentries_C=hNch_C->GetEntries();

  TH1D *hNch_P=(TH1D*)fin->Get("HNCharge_P_1");
  Int_t nentries_P=hNch_P->GetEntries();
  
  
  vector<float> cent_vec;

  
  
  for(int i=0;i<=100;i+=25){
    cent_vec.push_back(0.01*i);

  }

  //cent_vec.push_back(100);

  cout<<"Central Total entries = "<<nentries_C<<endl;
  cout<<"Peripheral Total entries = "<<nentries_P<<endl;
  cout<<"Index = "<<cent_vec.size()<<endl;

  const int NCent=cent_vec.size();
  int ncent = NCent-1;
  //Central
  int max_bin=hNch_C->FindLastBinAbove(0,1);
  
  Float_t Nch_bin_C[NCent][5];
 
  for(int i=0;i<NCent;i++){
    Nch_bin_C[i][0]=cent_vec.at(i);
    
    double a =Nch_bin_C[i][0]*nentries_C;
    Int_t mult=(Int_t) a;
    
    cout<<"sum number"<<mult<<endl;
    if(i==0) {Nch_bin_C[i][1]=mult;Nch_bin_C[i][2]=mult;}
    else{
      Nch_bin_C[i][1]=mult;
      mult=mult-Nch_bin_C[i-1][1];
      Nch_bin_C[i][2]=mult;
    }
    
    
  }
 
  
  float mult_sum=0;
  float bw=hNch_C->GetBinWidth(1);
  for(int i=max_bin;i>0;i--){
    float n=hNch_C->GetBinContent(i);
    
    mult_sum +=n;
    for(int k=1;k<ncent;k++){
	    if(mult_sum<=Nch_bin_C[k][1])
	  {
      Nch_bin_C[k][3]=hNch_C->GetBinCenter(i)+0.5*bw;
      Nch_bin_C[k][4]=mult_sum;
	  }
     
    }
    if(i==1) cout<<mult_sum<<endl;
  }
  Nch_bin_C[0][4]=0;
  Nch_bin_C[0][3]=hNch_C->GetBinCenter(max_bin)+0.5*bw;

  int first_bin=hNch_C->FindFirstBinAbove(0,1);
  Nch_bin_C[ncent][4]=nentries_C;
  Nch_bin_C[ncent][3]=hNch_C->GetBinCenter(first_bin)+0.5*bw;
  

  ofstream fout;
  fout.open("Nch_out_cental.txt");
  for(int i=0;i<NCent;i++){
    fout<<i<<"\t"<<"percent: "<<Nch_bin_C[i][0]<<"\t"<<"sum number: "<<Nch_bin_C[i][1]<<"\t"<<"num per bin: "<<Nch_bin_C[i][2]<<"\t"<<"real sum number: "<<Nch_bin_C[i][4]<<"\t\t"<<"Nch bin: "<<Nch_bin_C[i][3]<<endl;   
  }

  fout.close();


  //Peripheral
  int max_bin_P=hNch_P->FindLastBinAbove(0,1);
  
  Float_t Nch_bin_P[NCent][5];
 
  for(int i=0;i<NCent;i++){
    Nch_bin_P[i][0]=cent_vec.at(i);
    double a = Nch_bin_P[i][0]*nentries_P;
    Int_t mult=(Int_t) a;
    if(i==0) {Nch_bin_P[i][1]=mult;Nch_bin_P[i][2]=mult;}
    else{
      Nch_bin_P[i][1]=mult;
      mult=mult-Nch_bin_P[i-1][1];
      Nch_bin_P[i][2]=mult;
    }
    
    
  }
 
  
  float mult_sum_P=0;
  float bw_P=hNch_P->GetBinWidth(1);
  for(int i=max_bin_P;i>0;i--){
    float n=hNch_P->GetBinContent(i);
    
    mult_sum_P +=n;
    for(int k=1;k<ncent;k++){
	    if(mult_sum_P<=Nch_bin_P[k][1])
	  {
      Nch_bin_P[k][3]=hNch_P->GetBinCenter(i)+0.5*bw_P;
      Nch_bin_P[k][4]=mult_sum_P;
	  }
     
    }
    if(i==1) cout<<mult_sum_P<<endl;
  }
  Nch_bin_P[0][4]=0;
  Nch_bin_P[0][3]=hNch_P->GetBinCenter(max_bin_P)+0.5*bw_P;

  int first_bin=hNch_P->FindFirstBinAbove(0,1);
  Nch_bin_P[ncent][4]=nentries_P;
  Nch_bin_P[ncent][3]=hNch_P->GetBinCenter(first_bin)+0.5*bw_P;
  

  ofstream fout1;
  fout1.open("Nch_out_Peripheral.txt");
  for(int i=0;i<NCent;i++){
    fout1<<i<<"\t"<<"percent: "<<Nch_bin_P[i][0]<<"\t"<<"sum number: "<<Nch_bin_P[i][1]<<"\t"<<"num per bin: "<<Nch_bin_P[i][2]<<"\t"<<"real sum number: "<<Nch_bin_P[i][4]<<"\t\t"<<"Nch bin: "<<Nch_bin_P[i][3]<<endl;   
  }

  fout1.close();

  //hNch_C->Draw();

}
