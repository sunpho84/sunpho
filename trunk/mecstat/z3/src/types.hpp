#ifndef _TYPES_HPP
#define _TYPES_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "macros.hpp"

typedef int coords[NDIMS];

typedef unsigned short int z3_t;

//The structure for the random generator
class rnd_gen
{
public:
  int idum;
  int idum2;
  int iv[RAN2_NTAB];
  int iy;
};

struct N_t
{
  int N;
  int N0;
};

#endif
