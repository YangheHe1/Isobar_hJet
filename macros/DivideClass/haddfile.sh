#!/bin/sh

cd /star/data05/scratch/yanghe/mixclass/central/




N1=4
N3=2


    for((I1=0;I1<N1;I1+=1))
    do
        for((I3=0;I3<N3;I3+=1))
        do

            Filename="F_mixedevents_mult_${I1}_vz_$1_Psi_${I3}"
            File="${Filename}.root"
            echo "$File"
            cd $Filename
            hadd $File MC*.root
            mv $File ../RootFile
            cd ..

        done



    done


