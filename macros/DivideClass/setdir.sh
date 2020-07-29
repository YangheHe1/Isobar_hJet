#!/bin/sh

cd /star/data05/scratch/yanghe/mixclass/central/

N=(4 2 2)

for((I1=0;I1<N[0];I1+=1))
do
    for((I2=0;I2<N[1];I2+=1))
    do
        for((I3=0;I3<N[2];I3+=1))
        do

            Filename="F_mixedevents_mult_${I1}_vz_${I2}_Psi_${I3}"
            
            echo "$Filename"
            mkdir $Filename

        done



    done


done
