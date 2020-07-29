#include <TSystem>
#include <ctime>
#include <string>
#include <iostream>

class StMaker;
class StChain;
class StPicoDstMaker;

StChain *chain;


void runHFullJet(
	Int_t nevents=-1,
	const char* inputfilelist="file.list",
	const char*     out_name ="test",
	const char* bad_run_list ="badrun.list",
	const char* run_number_list="runnumber.list" )
{

	

	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gSystem->Load("StPicoEvent");
	gSystem->Load("StPicoDstMaker");
	gSystem->Load("StDbBroker");
	gSystem->Load("St_db_Maker.so");
	gSystem->Load("StRefMultCorr");	
	gSystem->Load("StEpdUtil");
	gSystem->Load("HFullJetAnalysis");
		
/*
    const char* database = "mysql://db04.star.bnl.gov:3414/RunLog?timeout=60";
    const char* user = "yanghe";
    const char* pass = "";
    TMySQLServer* mysql = TMySQLServer::Connect(database,user,pass);

*/


	string str;
	ifstream file;
	file.open(inputfilelist);
	ofstream fout;
	fout.open("picos.list");
	while (std::getline(file, str)){
		std::istringstream iss(str);
		string nf;
		vector<string> p__;
		while (iss >> nf)
		{
			p__.push_back(nf);
		}
		if(p__.size()<=0) continue;
		fout<<p__.at(0)<<"\n";
		std::cout <<p__.at(0)<<"\n";
	}
	fout.close();	






	chain = new StChain();

	StPicoDstMaker *picoMaker =0X0;
	StPicoDstMaker::PicoIoMode IoMode = 2;
	picoMaker = new StPicoDstMaker(IoMode,"picos.list","picoDst");

	HFullJetAnalysis *myanalysis = new HFullJetAnalysis(
		out_name,
		picoMaker,
		bad_run_list,
		run_number_list   );

	cout << "chain->Init()" <<endl;
	chain->Init();

	int total = picoMaker->chain()->GetEntries();
	cout << " Total entries = " << total << endl;

	int nEvents = nevents;
	if(nEvents>total||nEvents<0) nEvents = total;
	cout<<"Number of Events = "<<nEvents<<endl;

	for (Int_t i=0; i<nEvents; i++){
		if(i%50==0) {
			 cout << "- Finished " << i << " events." << endl;
		}
		chain->Clear();
		int iret = chain->Make(i);
		if (iret) { cout << "Bad return code!" << iret << endl; break;}
	}

	chain->Finish();

	delete chain;

}


























