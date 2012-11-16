#pragma once 

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE

#include <TMinuit.h>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

extern "C" { void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
			 double* w, double* work, int* lwork, int* info );}

#include "tools.cpp"

#include "endianess.cpp"
#include "functions.cpp"
#include "file.cpp"

#include "jack.cpp"
#include "boot.cpp"

#include "jvec.cpp"
#include "bvec.cpp"

#include "grace.cpp"
#include "plot.cpp"
#include "fits.cpp"

#include "rand.cpp"
#include "interpolate.cpp"

#include "gevp.cpp"

