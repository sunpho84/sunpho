#pragma once 

#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE

#include <TMinuit.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>

using namespace std;

bool debug_load=true;

#include "tools.cpp"
#include "auto.cpp"

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

