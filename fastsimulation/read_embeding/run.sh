#!/bin/sh
order=$1
export PARTICLE="Pi"
export FILELIST="Pifilelist/Pi_file_$order.list"
export CENTRAL=2
export GLOBAL=0
export OUTPUT=$order
export NJOBS=$2

starver SL18h
pwd
#ls -latr
root4star -l <<EOF
.L StMiniMcTree.C
StMiniMcTree uu
uu.Loop()
EOF
#.q
