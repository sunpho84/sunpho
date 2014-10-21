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

EXTERN_RANDOM mt19937_64 gen;

inline double get_unif_double(double max)
{
  uniform_real_distribution<double> dis(0,max);
  return dis(gen);
}
inline double get_gauss_double()
{
  normal_distribution<double> dis(0,1);
  return dis(gen);
}
inline dcomplex get_gauss_complex()
{
  normal_distribution<double> dis(0,M_SQRT_2);
  return dcomplex(dis(gen),dis(gen));
}
void set_U1_to_rnd(dcomplex &U);
double get_theta(double a,int k);
double get_theta_1(double a);
void set_ON_to_rnd(dcomplex *O);

#endif
