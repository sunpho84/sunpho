#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#ifndef EXTERN_RANDOM
 #define EXTERN_RANDOM extern
#endif

#include <random>

using namespace std;

#include "types.hpp"

#ifndef M_SQRT_2
 #define M_SQRT_2 0.707106781186547524401
#endif

EXTERN_RANDOM mt19937_64 *gen;

inline double get_unif_double(double max,int site)
{
  uniform_real_distribution<double> dis(0,max);
  return dis(gen[site]);
}
inline double get_gauss_double(int site)
{
  normal_distribution<double> dis(0,1);
  return dis(gen[site]);
}

#endif
