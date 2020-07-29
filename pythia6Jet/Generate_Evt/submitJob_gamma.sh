#!/bin/bash

# check first your previous" val " and change it for next run otherwise it would overwrite your previously generated root files 
val=400

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do


n=`expr $val + $i`

  echo "submitting job -- $n"
      ./Execute.sh $n
  echo "Fined job -- $n"
#       -q
#  root4star -b -q 'RunPythia_GammaProcess.C(100,2,30,$i,370,$i)';

done
