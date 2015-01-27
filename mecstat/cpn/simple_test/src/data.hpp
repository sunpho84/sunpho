#ifndef _DATA_HPP
#define _DATA_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <vector>

#include "parameters.hpp"
#include "types.hpp"

using namespace std;

#ifndef EXTERN_DATA
 #define EXTERN_DATA extern
#endif

EXTERN_DATA dcomplex *zeta,*lambda;
EXTERN_DATA dcomplex *zeta_old,*lambda_old;
EXTERN_DATA dcomplex *pi,*fpi;
EXTERN_DATA double *omega,*fomega;
EXTERN_DATA dcomplex *topo_staples_data,*topo_staples_supp_data;
EXTERN_DATA dcomplex **lambda_stout;

void copy_zeta_conf(dcomplex *dest,dcomplex *source);
void copy_lambda_conf(dcomplex *dest,dcomplex *source);

#endif
