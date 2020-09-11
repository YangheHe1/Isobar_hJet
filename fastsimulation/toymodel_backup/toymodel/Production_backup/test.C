#include "Efficiency.h"
void test(bool cent=1, TString ratio="pp", int cutset=1)
{

	TH1D* heffi=new TH1D("heffi","efficiency",40,0,10);
	int N=5000;
	for(int i=0; i<N; i++){
		double pt=gRandom->Uniform(0.1,10);
		double epsilon=efficiency11(pt, cent, ratio, cutset,0);
		cout<<"epsilon="<<epsilon<<endl;
		heffi->Fill(pt,epsilon);
	}
	heffi->Scale(10.0/N,"width");
	//heffi->Draw();
}
