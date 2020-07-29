#!/bin/bash
root -b <<EOF

gSystem->Load("readtree_C.so")
readtree* myana = new readtree("fileP.list", $1, $2, $3);
myana->InitHist();
myana->Loop();
.q
EOF

