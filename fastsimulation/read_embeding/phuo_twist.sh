#!/bin/sh
echo "RUNNING PROGRAM"

if [ -e log ];then
    rm -rf log
    echo "We now delete log"
fi
mkdir log


N0=36

for((I1=0;I1<N0;I1+=1))
do
    
    echo $I1 
    Filenm01=FBasy_$I1.job
    Filenm02=FBasy_$I1.sh

    cp condor_script $Filenm01
    echo "Output          = log/Corrjob_$I1.out">>$Filenm01
    echo "Error           = log/Corrjob_$I1.err">>$Filenm01
    echo "Log             = log/Corrjob_$I1.log">>$Filenm01
    echo "Executable      = $Filenm02">>$Filenm01
    echo "queue"                        >>$Filenm01

    echo "#!/bin/sh    ">$Filenm02
    echo 'echo $0'>>$Filenm02
    echo "sh run.sh $I1 $N0 >log/job_$I1.log">>$Filenm02
    echo "rm $Filenm01">>$Filenm02
    echo "rm $Filenm02">>$Filenm02

    chmod 711 $Filenm02
    condor_submit $Filenm01 
done
