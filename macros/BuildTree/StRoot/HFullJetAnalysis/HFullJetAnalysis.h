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

const Int_t MAXTRACK=10000;

class HFullJetAnalysis : public StMaker 
{

	private:
	StPicoDstMaker* picoReader;
	St_db_Maker*    starDb;
	vector<int>     badrunlist;
	vector<int>     runidlist;
	TTree          *outTree;
	
	Double_t	QxPos[2][3];
	Double_t	QyPos[2][3];
	Double_t	QwPos[2][3];
	Double_t	QxNeg[2][3];
	Double_t	QyNeg[2][3];
	Double_t	QwNeg[2][3];
	

	Int_t		RunId;
	Int_t		EventId;
	Short_t    refmult;
	Float_t    Vx;
	Float_t    Vy;
   	Float_t    Vz;
	Float_t	zdc;
	Float_t	bbc;
	Float_t	Vz_VPD;
	Float_t 	NBTOFMultfit;

	//TPC information
	Float_t 	NCharge;
	Int_t	 	numTrk;
   	
	
	Float_t		Pt[MAXTRACK];
	Float_t		Px[MAXTRACK];
	Float_t		Py[MAXTRACK];
	Float_t		Pz[MAXTRACK];
	Float_t		Eta[MAXTRACK];
	Float_t		Phi[MAXTRACK];
	Float_t  		Charge[MAXTRACK];
	Float_t		Dca[MAXTRACK];
	Float_t		nHits[MAXTRACK];

	//EPD hit information
    Int_t 		EPDmaxHit;
    Short_t 		NumEPDHit;
    
    Float_t 	NumMip[MAXTRACK];
    Int_t 		TileID[MAXTRACK];
    
    Float_t 	EPDX[MAXTRACK];
    Float_t 	EPDY[MAXTRACK];
    Float_t 	EPDZ[MAXTRACK];
	Float_t		EPD_PP[MAXTRACK];
	Float_t		EPD_TT[MAXTRACK];
	Float_t		EPD_eta[MAXTRACK];
	Float_t		EPD_phi[MAXTRACK];


	//BBC information
	Int_t    mBbcQ[48];
	
	Double_t  	EPD_Psi;
	Double_t  	TPC_Psi;

	string          infile2_;
	TFile*          fout;
	
	

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

	ClassDef(HFullJetAnalysis, 1)
};

#endif











































