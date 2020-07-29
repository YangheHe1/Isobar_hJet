#!/bin/bash
root -b <<EOF

gSystem->Load("readjet_C.so")
readjet* myana = new readjet("fileR4.list", $1, $2);
myana->InitHist();
myana->Loop();
.q
EOF

