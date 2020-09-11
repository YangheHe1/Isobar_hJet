void run_pythiaEmb()
{
  gSystem->Load("$TOYMODELDIR/Production/libThrm.so");
  //gSystem->Load("$TOYMODELDIR/Analysis/libThrmAna.so");
  //gSystem->Load("$TOYMODELDIR/Unfolding/libtoyUnfold.so");
  
  TStopwatch timer;
  timer.Start();

Int_t nevents=atoi(gSystem->Getenv("NEVENTS"));

ThrmPythiaEmb *pyemb=new ThrmPythiaEmb();
pyemb->Run(nevents);

  timer.Stop();
  timer.Print();

}
