#ifndef _TYPES_HPP
#define _TYPES_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <complex>

#include "macros.hpp"

using namespace std;

typedef complex<double> dcomplex;
typedef int coords[NDIMS];

//The structure for the random generator
class rnd_gen
{
public:
  int idum;
  int idum2;
  int iv[RAN2_NTAB];
  int iy;
};

#endif
