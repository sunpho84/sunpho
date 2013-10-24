#!/bin/bash

wget http://www.netlib.org/lapack/lapack-3.4.2.tgz -O -|tar xzf -
cd lapack*/
mv make.inc.example make.inc
make all
make
make -j 8 man
