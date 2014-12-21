#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef EXTERN_PARAMETERS
#define EXTERN_PARAMETERS extern
#endif

#include "types.hpp"

//parameters
EXTERN_PARAMETERS int N;
EXTERN_PARAMETERS int L;
EXTERN_PARAMETERS int seed;
EXTERN_PARAMETERS int nsweep;
EXTERN_PARAMETERS int nterm;
EXTERN_PARAMETERS start_cond_t start_cond;
EXTERN_PARAMETERS int nhmc_steps;

#endif
