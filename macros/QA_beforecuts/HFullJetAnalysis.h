#ifndef HFullJetAnalysis_def
#define HFullJetAnalysis_def


#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TRandom.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "St_db_Maker/St_db_Maker.h"
#include "StEmcUtil/database/StBemcTables.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/projection/StEmcPosition.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include <map>


#include "StEvent/StDcaGeometry.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;

class TString;
class TH1F;
class TH1D;
class TProfile;
const Int_t runmin=0;
const Int_t runmax=80;
const Int_t runbins=80;
const Int_t NCENT = 9;
const Int_t MAXTRACK=10000;

class HFullJetAnalysis : public StMaker 
{

	private:
	StPicoDstMaker* picoReader;
	St_db_Maker*    starDb;
	vector<int>     badrunlist;
	vector<int>     runidlist;
	

	string          infile2_;
	TFile*          fout;
	
	TH2D*			hVtxXvsY[2];
	TH1D*			HVz[2];
	TH1D*			refMultHist[3];
	TH2D*			vz_refmult[2];
	TH1D*			etahist[2];
	TH1D*			phihist[2];
	TH1D*			Pt_all_hist[2];
	TH1D*			Hdelta_Vz[2];
	TH2D*			H2d_Vz[2];
	TH2D*			Htof_ref[2];
	TH2D*			Htof_ntrk[2];


	TH1D*			Pthist_C;
	TH1D*			Pthist_P;
	TH1D*			Pthist_Cut0;
	TH1D*			Pthist_Cut1;
	TH1D*			Pthist_Cut2;
	TH1D*			Pthist_Cut3;
	TH1D*			nutrigger;
	TH1D*			triggerdist_C[2];
	TH1D*			triggerdist_P[2];
	TH1D*			eventC;
	TH1D*			eventP;
	TH1D* 			HNCharge_C[2];
	TH1D* 			HNCharge_P[2];
	

	

	public:
	HFullJetAnalysis(
		const char*     Out,
		StPicoDstMaker* picoDstMaker,
		const char*     _bad_run_list,
		const char*     _run_number_list
		);

	virtual ~HFullJetAnalysis();
	virtual Int_t Init();
	virtual Int_t Make();
	virtual void  Clear(Option_t *opt="");
	virtual Int_t Finish();
	bool refmult_check(int __nBTOFMatch, int __refMult);
	bool ntrk_check(int __nBTOFMatch, int __ntrk);
	Int_t centrality( float _refMult );

	ClassDef(HFullJetAnalysis, 1)
};

#endif











































