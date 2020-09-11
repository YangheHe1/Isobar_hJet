#!/bin/bash

source set_paths.sh
WPATH="$TOYMODELDIR/DataOut/test"
export OUTPUTDIR=$WPATH
export WRKDIR=$WPATH

export NEVENTS=5000

export PTCUT=0.2
export RPARAM=0.2
export MAXRAP=1.0
export AREACUT=0.09
export NJOBS=1

export CENTRAL=1 #central or peripheral collisions
export CHARGED=1 #charged jets only
export JETFRAG="2u1g" #jet fragmentation: u | g
export EFFICORR=1 #apply tracking efficiency 
export PTSMEAR=1 #track pT smearing
export MOMRES=2 #momentum resolution parametrization 0: global tracks, 1: primary tracks, 2: primary tracks, more accurate
#export CUTSET=1 #set of track cuts which we are using (see ~/jet_analysis/STARJet/Analysis/Embedding/out_cutsetX/cuts.h for details on cuts)
export EFFTYPE="pp" #tracking efficiency model: pp | AuAu - like
#export EFFPATH="$ANALYSISDIR/efficiency" #path to tracking efficiency files
export EFF_INCREMENT=0 #increase/decrease tracking efficiency for systematic studies
export DOSCALE=0 #rescale charge jets to the momentum of created full jet
export TOFEFFI=0 #correc also for TOF+BEMC maching efficiency (for pp study)
export TOFEFF_INCREMENT=0 #increase/decrease TOF+BEMC matching efficiency for systematic studies
export FIXEDSEED=1 #1: use SEED as a value for the random number generator's seed 0: use computer time as a seed
export SEED=42 #seed for random number generator

root -l -q -b run_pythiaEmb.C 
#done
