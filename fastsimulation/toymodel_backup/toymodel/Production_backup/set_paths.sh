#!/bin/bash
export FASTJETDIR="/star/u/yanghe/pwg/Analysis_software/fastjet3"
export PATH="$PATH:$FASTJETDIR/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$FASTJETDIR/lib"

export TOYMODELDIR="/star/u/yanghe/Isobar_200_jet/toymodel"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$TOYMODELDIR/Production"
