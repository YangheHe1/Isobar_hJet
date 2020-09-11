#!/bin/sh
echo "RUNNING PROGRAM"

if [ -e log1 ];then
    rm -rf log1
    echo "We now delete log"
fi
mkdir log1

N0=425

A_cut_R2=0.05
A_cut_R3=0.2
A_cut_R4=0.35


for((I1=0;I1<N0;I1+=1))
do
    I2=`expr  $I1 +  1`
    echo $I1 $I2
    Filenm01=FBasy_$I1.job
    Filenm02=FBasy_$I1.sh

    cp condor_script $Filenm01
    echo "Output          = log1/Corrjob_$I1.out">>$Filenm01
    echo "Error           = log1/Corrjob_$I1.err">>$Filenm01
    echo "Log             = log1/Corrjob_$I1.log">>$Filenm01
    echo "Executable      = $Filenm02">>$Filenm01
    echo "queue"                        >>$Filenm01

    echo "#!/bin/sh    ">$Filenm02
    echo 'echo $0'>>$Filenm02
    echo "sh run.sh $I1 $I2 $A_cut_R2 >log1/job_$I1.log">>$Filenm02
    echo "rm $Filenm01">>$Filenm02
    echo "rm $Filenm02">>$Filenm02

    chmod 711 $Filenm02
    condor_submit $Filenm01 
done
