#pragma once

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define nmom 140

typedef double complex[2];
typedef complex cmom[140];

typedef complex spin[4];
typedef spin spinspin[4];
typedef spinspin ss_propagator[140];

typedef complex color[3];
typedef color su3[3];
typedef su3 cc_propagator[140];

typedef spinspin colorspinspin[3];
typedef colorspinspin su3spinspin[3];

typedef su3spinspin ccss_propagator[nmom];

typedef spin colorspin[3];
typedef colorspin spin_colorspin[3];
typedef spin_colorspin colorspin_colorspin[4];
typedef colorspin_colorspin ape_propagator[nmom];

const double PI=3.14159265358979323846;
const double bc[4]={1,0,0,0};

int big_endian;
int T,L;

double P[4][nmom],SinP[4][nmom];
double P2[nmom],SinP2[nmom],SinP4[nmom];

#include "complex.c"
#include "routines.c"
#include "invert_matrix.c"
#include "propagator.c"
#include "read.c"
#include "momentum.c"
#include "init.c"

