#ifndef _DATA_HPP
#define _DATA_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "parameters.hpp"
#include "types.hpp"

using namespace std;

#ifndef EXTERN_DATA
 #define EXTERN_DATA extern
#endif

EXTERN_DATA z3_t *phi;
EXTERN_DATA double act_contr_tab[(2*NDIMS+1)*(2*NDIMS+1)][3];

#endif
