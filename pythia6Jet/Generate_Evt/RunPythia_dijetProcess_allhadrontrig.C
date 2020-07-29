void RunPythia_dijetProcess_allhadrontrig(const int nevents = 40000, 
		     const float ptmin = 4, 
			    const float ptmax = 100, //no use of high pThatmas
		     const int seed = 30, 
		     const float val = 370, 
		     const int iout = 0
			    )
{
  TString outputPyth;
  outputPyth = "/gpfs/mnt/gpfs01/star/pwg/yanghe/Pythia6/";
  outputPyth +="hTrig_pTHat4to100/hTrig_pTHat4to100";
  outputPyth +="_";
  outputPyth += iout;
  outputPyth +=".root";
  //  outputPyth +=

  Double_t TrigPTCut = 7.0; // to set trigger event to store in the root files
  
  gSystem->Load("StJetSkimEvent");
  gSystem->Load("StPythiaRecord");
  gSystem->Load("libfastjet.so");
  gSystem->Load("libfastjettools.so");
  gSystem->Load("libsiscone.so");
  gSystem->Load("libsiscone_spherical.so");
  gSystem->Load("libfastjetplugins.so");
  gSystem->Load("StJetFinder");
  
  gSystem->Load("StPythiaJet");
  assert(gSystem->Load("/star/u/zchang/software/local/lib/libLHAPDF.so") == 0);
  assert(gSystem->Load("/star/u/zchang/software/local/lib/libPythia6-6.4.28.so") == 0);
  //assert(gSystem->Load("$HOME/software/local/lib/libPythia6-6.4.28.so") == 0);
  //assert(gSystem->Load("libPythia6_4_26.so") == 0);

  // Create an instance of Pythia
  TPythia6* pythia = new TPythia6;

  printf("Seed = %d\n",seed);

  // Physics processes
   pythia->SetMSEL(1);		// QCD jets
  //  pythia->SetMSEL(10);		// QCD jets
  pythia->SetMRPY(1,seed);      // Seed random number generator
  pythia->SetCKIN(3,ptmin);	// Lower partonic pT bound in GeV
//  pythia->SetCKIN(4,ptmax);	// Higher partonic pT bound in GeV
  //pythia->SetMSTP(5,100);	// CDF tune A
  //pythia->SetMSTP(5,320);	// Perugia 0 tune
  //pythia->SetMSTP(5,326);	// Perugia 6 tune
  //pythia->SetMSTP(5,370);	// Perugia 2012 tune Pythia 6.4.28
  pythia->SetMSTP(5,val);	// Perugia 2012 tune Pythia 6.4.28

  // Make the following stable
  pythia->SetMDCY(102,1,0);  // PI0 111
  //  pythia->SetMDCY(106,1,0);  // PI+ 211
  //  pythia->SetMDCY(109,1,0);  // ETA 221
  //  pythia->SetMDCY(116,1,0);  // K+ 321
  //  pythia->SetMDCY(112,1,0);  // K_SHORT 310
  //  pythia->SetMDCY(105,1,0);  // K_LONG 130
  //  pythia->SetMDCY(164,1,0);  // LAMBDA0 3122
  //  pythia->SetMDCY(167,1,0);  // SIGMA0 3212
  //  pythia->SetMDCY(162,1,0);  // SIGMA- 3112
  //  pythia->SetMDCY(169,1,0);  // SIGMA+ 3222
  //  pythia->SetMDCY(172,1,0);  // Xi- 3312
  //  pythia->SetMDCY(174,1,0);  // Xi0 3322
  //  pythia->SetMDCY(176,1,0);  // OMEGA- 3334

    pythia->Initialize("cms","p","p",200); // p+p collisions at sqrt(s)=200 GeV
  //  pythia->Initialize("cms","p","p",510); // p+p collisions at sqrt(s)=500 GeV
  
    pythia->SetMSTP(5, 0);
  // Choice of PDF
  //  pythia->SetMSTP(51,7); // CTEQ 5L
  //  pythia->SetMSTP(52,1); // Pythia Internal PDF
  //pythia->SetMSTP(52, 2); // libPDF
  //pythia->SetMSTP(51, 263000); //CTEQ 6M 1000*10000 + 0 (1000 * NGroup + NSet)
  //263000 NNPDF 3.0 LO as0130
  //262000 NNPDF 3.0 LO as0118
  //200200 NNPDF 2.1 LO as0119
  //230000 NNPDF 2.3 NLO as 0119
  //10000 cteq6m (not validated)
  //10042 cteq6l1

  //cout<<"PARP(82) from "<<pythia->GetPARP(82)<<" to ";
  //pythia->SetPARP(82, pythia->GetPARP(82)*1.073);
  pythia->SetPARP(90, 0.213); //exponent for pT_0 which controls underlying events
  //pythia->SetPARP(91, 1.0); // set primordial k_T (default for Perugia 0 is 2)
  //pythia->SetMSTP(81, 20); //turn off multiple interactions
  //pythia->SetPARP(82, pt0); // set pT_0 cutoff pT for underlying events
  //pythia->SetMSTP(61, 0); //initial state radiation
  //  pythia->Initialize("cms","p","p",510); // p+p collisions at sqrt(s)=500 GeV
  pythia->Initialize("cms","p","p",200); // p+p collisions at sqrt(s)=500 GeV
   
  //pythia record maker
  StPythiaRecordMaker *recordmaker = new StPythiaRecordMaker;
  //  recordmaker->SetFile(outrecord);
  recordmaker->SetFile(outputPyth.Data());
  //recordmaker->SetSkim(1);
  if(recordmaker->Init()) printf("Pythia record Maker Initialized\n");

  //  if(jetmaker->Init()) printf("Pythia jet maker initialized\n");
   TClonesArray* particles = new TClonesArray("TParticle", 1000);
   
  // Event loop
   Int_t eventcounter=0;
  for (int iEvent = 1; iEvent <= nevents; ++iEvent) {
    pythia->GenerateEvent();
    //    if (iEvent % 10000 <= 10) pythia->Pylist(1);
    //    recordmaker->Make(iEvent, pythia);
    //    jetmaker->Make(iEvent, pythia);

    //__________
    pythia->ImportParticles(particles,"All");
    Int_t np = particles->GetEntriesFast();
    
    //    cout<<"MSTP: "<<pythia->GetMSTP(5)<<endl;
    //    cout<<" np: "<<np<<endl;
    Int_t flagtrigge = 0; //flagging trigger events 
    // Particle loop
    for (Int_t ip = 0; ip < np; ip++) {
      TParticle* part = (TParticle*) particles->At(ip);
      Int_t ist = part->GetStatusCode();
      // Positive codes are final particles.
      if (ist <= 0) continue;
      Int_t pdg = part->GetPdgCode();
      //Float_t charge = TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      //      if (charge == 0.) continue;
      //      Float_t eta = part->Eta();
      Float_t pt  = part->Pt();
      //      mParticles_fStatusCode      
      int stcode = part->GetStatusCode();
      
      //      if(pdg ==111 || pdg == 211 || pdg == -211) 

            if(pdg == 211 || pdg == -211
         || pdg == 321 || pdg == -321
         || pdg == 2212 || pdg == -2212
	       || pdg == 111)      // for all hadrons	      
        {
	  //  	  if(stcode==1){
	  if( pt > TrigPTCut && stcode == 1) {
	    //	  cout<<pdg<<"  "<<pt<<"  "<<stcode<<endl;
	  flagtrigge = 1;
	   }	  
	}      // trigger selection gamma


    } //particle loop

    if(flagtrigge==1){
      //cout<<"flagtrigge :"<<flagtrigge<<endl;
      recordmaker->Make(iEvent, pythia);
      eventcounter++;
      if(eventcounter % 100 == 0)cout<<" trigger event# "<<eventcounter<<endl;
    }

     
  } // End event loop
   
  //write output
  recordmaker->Finish();
  //  jetmaker->Finish();
  
  pythia->Pystat(1);		// Print PYTHIA statistics
  printf("PARI[1] = %e\n",pythia->GetPARI(1));
  printf("Tune: MSTP[5] = %f\n",pythia->GetMSTP(5));
  printf("PDF: MSTP[51] = %f\n",pythia->GetMSTP(51));
  printf("PARP[90]: PARP[90] = %f\n",pythia->GetPARP(90));
  printf("kT: PARP[91] = %f\n",pythia->GetPARP(91));
  printf("pT0: PARP[82] = %f\n",pythia->GetPARP(82));
  printf("ue qq enh: PARP[87] = %f\n",pythia->GetPARP(87));
  printf("ue qq enh scale: PARP[88] = %f\n",pythia->GetPARP(88));
  printf("Q^2 dependence: MSTP[57] = %f\n",pythia->GetMSTP(57));
  printf("ISR: MSTP[61] = %f\n",pythia->GetMSTP(61));
  printf("UE: MSTP[81] = %f\n",pythia->GetMSTP(81));

 }
