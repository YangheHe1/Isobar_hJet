#!/bin/bash
STARLIB_VER="SL18c_embed"
BASEPATH="/star/u/rusnak/JET_analysis_run11/STARJet"

export FASTJETDIR="$BASEPATH/software_SL19c/fastjet3"
export PATH="$PATH:$FASTJETDIR/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$FASTJETDIR/lib"

export TOYMODELDIR="/star/u/yanghe/Isobar_200_jet/toymodel"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TOYMODELDIR/Production"
