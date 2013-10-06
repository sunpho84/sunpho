#!/bin/bash

wget http://www.fftw.org/fftw-3.3.3.tar.gz -O -|tar xzvf -
cd fftw*/
./configure --prefix=$HOME
make -j 8
make install
cd ..
rm -fr fftw*/