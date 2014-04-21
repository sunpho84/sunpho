#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#ifndef EXTERN_RANDOM
 #define EXTERN_RANDOM extern
#endif

#include "types.hpp"

//random number generators
#ifdef GOOD_GENERATOR
 EXTERN_RANDOM random_device *rd;
 EXTERN_RANDOM mt19937_64 *gen;
 EXTERN_RANDOM uniform_real_distribution<double> *dis;
#else
 EXTERN_RANDOM rnd_gen gen;
#endif

double get_unif_double(double max,bool incl=false);
void set_U1_to_rnd(dcomplex &U);
double get_theta(double a,int k);
double get_theta_1(double a);
void set_ON_to_rnd(dcomplex *O);
void init_system_to(int cond);
#endif
