# n=$1
echo $1
 root4star -q -b "RunPythia_dijetProcess_allhadrontrig.C(1000000,4,100,$1,370,$1)"
#(Event, pThatmin, pThatmax,seed,val,rootfilenumber)
