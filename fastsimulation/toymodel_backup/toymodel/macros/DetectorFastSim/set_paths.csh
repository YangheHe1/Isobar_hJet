#!/bin/csh
#STARJet
set STARLIB_VER = SL18c_embed
set BASEPATH = /star/u/rusnak/JET_analysis_run11/STARJet

#other software
setenv ROOUNFOLD $BASEPATH/software_SL19c/RooUnfold/v-trunk-custom
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$ROOUNFOLD

setenv FASTJETDIR $BASEPATH/software_SL19c/fastjet3
setenv PATH $PATH\:$FASTJETDIR/bin
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$FASTJETDIR/lib

#Parametrized model
setenv TOYMODELDIR /star/u/yanghe/Isobar_200_jet/toymodel
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$TOYMODELDIR/Production
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$TOYMODELDIR/Analysis
setenv LD_LIBRARY_PATH $LD_LIBRARY_PATH\:$TOYMODELDIR/Unfolding
