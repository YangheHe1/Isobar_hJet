#!/bin/bash
root -b <<EOF

gSystem->Load("readjet_C.so")
readjet* myana = new readjet("fileR2.list", $1, $2, $3);
myana->InitHist();
myana->Loop();
.q
EOF

