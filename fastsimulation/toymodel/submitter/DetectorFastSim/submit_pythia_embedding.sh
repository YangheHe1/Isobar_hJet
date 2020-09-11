#!/bin/bash
source set_paths.sh
LOGDIR="$TOYMODELDIR/submitter/log"
ERRDIR="$TOYMODELDIR/submitter/err"
MACRODIR="$TOYMODELDIR/macros/DetectorFastSim"
MACROFILE="run_pythiaEmb.C"

kEVENTS=2000 #how many thousands of events do we want to run

export PTCUT=0.2
export MAXRAP=1.0
export CHARGED=1 #charged jets only
export EFFICORR=1 #apply tracking efficiency 
export PTSMEAR=1 #track pT smearing
export MOMRES=2 #TPC momentum resolution; 0: sigma=0.01*pT^2 (global tracks) 1: sigma=0.005*pT^2 (primary tracks, simple) 2: sigma=a+b*pT+c*pT^2 (primary tracks, more accurate) 3: sigma=0.003*pT^2 (for sys uncertainty)
export EFFTYPE="AuAu" #tracking efficiency model: pp | AuAu - like
#export EFFPATH="$HOME/jet_analysis/STARJet/analysis/efficiency" #path to tracking efficiency files
export EFF_INCREMENT=0 #increase/decrease tracking efficiency for systematic studies, [-100%,100%]
export TOFEFFI=0 #correc also for TOF+BEMC maching efficiency (for pp study)
export TOFEFF_INCREMENT=0 #increase/decrease TOF+BEMC matching efficiency for systematic studies
export JETFRAG="2u1g" #jet fragmentation: u | g | 2u1g
#export CUTSET=1 #set of track cuts which we are using (see ~/jet_analysis/STARJet/Analysis/Embedding/out_cutsetX/cuts.h for details on cuts), cutset=3: same as cutset 1, but for global tracks instead of primary tracks
NFIT=15 #number of TPC fit points used in real data analysis (used only for description)
export DOSCALE=0 #rescale charge jets to the momentum of created full jet
export FIXEDSEED=1 #1: use SEED as a value for the random number generator's seed 0: use computer time as a seed
export SEED=42 #seed for random number generator
export Do3D_Eff=0 #0 pT eff; 1 pt eta phi eff
#SUFFIX="_TOF${TOFEFF_INCREMENT}_pp"
SUFFIX="_GPC" #dataset label

for CENTRAL in 1 0 #central or peripheral collisions
do
	export CENTRAL

for RPARAM in 0.2
do
export RPARAM

   if [ $RPARAM == "0.2" ]; then
      ACUT=0.07
   fi
   if [ $RPARAM == "0.3" ]; then
      ACUT=0.2
   fi
   if [ $RPARAM == "0.4" ]; then
      ACUT=0.4
   fi
   if [ $RPARAM == "0.5" ]; then
      ACUT=0.65
   fi

	export AREACUT=$ACUT
	echo "Area cut: $AREACUT"

	if [ $EFFICORR -eq 1 ]; then
   TYPE1="_effcorr"
   else
   TYPE1=""
   fi

   if [ $PTSMEAR -eq 1 ]; then
   TYPE2="_pTsmear"
   else
   TYPE2=""
   fi

	if [ $CENTRAL -eq 0 ]; then
	CENTRALITY="_peripheral"
	else
	CENTRALITY="_central"
	fi

	TYPE="pyEmb_${kEVENTS}k_charged_R${RPARAM}${TYPE1}${TYPE2}${CENTRALITY}_${JETFRAG}_eff${EFF_INCREMENT}_${EFFTYPE}_nfit${NFIT}_momRes${MOMRES}${SUFFIX}"

	export NEVENTS=100000 #number of events per job
   START=0
   export NJOBS=$(( kEVENTS * 1000 / NEVENTS ))
   #START=75
   #NJOBS=100
   echo "number of jobs: $NJOBS"

	export OUTPUTDIR=$SCRATCH
	OUT_PATH_FINAL="$TOYMODELDIR/DataOut/pythia/jetonly/$TYPE"
	if [ ! -e $OUT_PATH_FINAL ]; then
	  	mkdir -p $OUT_PATH_FINAL
	fi

	if [ ! -e tmp ]; then
		mkdir -p tmp
	fi

	TEMPLATE_NAME="pyEmb_cent${CENTRAL}_R${RPARAM}.xml"

#===========================
#create submission xml file
#===========================
echo "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" > tmp/$TEMPLATE_NAME
echo "<job nProcesses=\"$NJOBS\" simulateSubmission = \"false\" >" >> tmp/$TEMPLATE_NAME
echo "<command>" >> tmp/$TEMPLATE_NAME
echo setenv PTCUT $PTCUT >> tmp/$TEMPLATE_NAME
echo setenv MAXRAP $MAXRAP >> tmp/$TEMPLATE_NAME
echo setenv CHARGED $CHARGED >> tmp/$TEMPLATE_NAME
echo setenv EFFICORR $EFFICORR >> tmp/$TEMPLATE_NAME
echo setenv PTSMEAR $PTSMEAR >> tmp/$TEMPLATE_NAME
echo setenv MOMRES $MOMRES >> tmp/$TEMPLATE_NAME
echo setenv EFFTYPE $EFFTYPE >> tmp/$TEMPLATE_NAME
echo setenv EFF_INCREMENT $EFF_INCREMENT >> tmp/$TEMPLATE_NAME
echo setenv TOFEFFI $TOFEFFI >> tmp/$TEMPLATE_NAME
echo setenv TOFEFF_INCREMENT $TOFEFF_INCREMENT >> tmp/$TEMPLATE_NAME
echo setenv FIXEDSEED $FIXEDSEED >> tmp/$TEMPLATE_NAME
echo setenv SEED $SEED >> tmp/$TEMPLATE_NAME
echo setenv JETFRAG $JETFRAG >> tmp/$TEMPLATE_NAME
#echo setenv CUTSET $CUTSET >> tmp/$TEMPLATE_NAME
echo setenv DOSCALE $DOSCALE >> tmp/$TEMPLATE_NAME
echo setenv Do3D_Eff $Do3D_Eff >> tmp/$TEMPLATE_NAME
echo setenv CENTRAL $CENTRAL >> tmp/$TEMPLATE_NAME
echo setenv RPARAM $RPARAM >> tmp/$TEMPLATE_NAME
echo setenv AREACUT $AREACUT >> tmp/$TEMPLATE_NAME
echo setenv NEVENTS $NEVENTS >> tmp/$TEMPLATE_NAME
echo setenv NJOBS $NJOBS >> tmp/$TEMPLATE_NAME
echo setenv OUTPUTDIR $OUTPUTDIR >> tmp/$TEMPLATE_NAME
echo "  starver $STARLIB_VER" >> tmp/$TEMPLATE_NAME
echo "  source $TOYMODELDIR/set_paths.csh" >> tmp/$TEMPLATE_NAME
echo "  cd $MACRODIR"
echo "  pwd" >> tmp/$TEMPLATE_NAME
echo "  root4star -q -b -l $MACROFILE" >> tmp/$TEMPLATE_NAME
echo "  </command>" >> tmp/$TEMPLATE_NAME
echo "  <stdout URL=\"file:$LOGDIR/\$JOBID.log\"/>" >> tmp/$TEMPLATE_NAME
echo "  <stderr URL=\"file:$ERRDIR/\$JOBID.err\"/>" >> tmp/$TEMPLATE_NAME
echo "  <output fromScratch=\"*.root\" toURL=\"file:$OUT_PATH_FINAL/\"/>" >> tmp/$TEMPLATE_NAME
echo "  <SandBox>" >> tmp/$TEMPLATE_NAME
echo "    <Package>" >> tmp/$TEMPLATE_NAME
echo "	      <File>file:$MACRODIR/$MACROFILE</File>" >> tmp/$TEMPLATE_NAME
echo "			  </Package>" >> tmp/$TEMPLATE_NAME
echo "			  </SandBox>" >> tmp/$TEMPLATE_NAME
echo "			  </job>" >> tmp/$TEMPLATE_NAME

#let's submit
	cd tmp
	star-submit $TEMPLATE_NAME 
	cd ..


done #loop over R
done #loop ove centralities
