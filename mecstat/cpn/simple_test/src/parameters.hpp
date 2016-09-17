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
EXTERN_PARAMETERS double beta;
EXTERN_PARAMETERS int L;
EXTERN_PARAMETERS double g;
EXTERN_PARAMETERS int nstout_lev;
EXTERN_PARAMETERS double stout_rho;
EXTERN_PARAMETERS int nhmc_steps;
EXTERN_PARAMETERS int use_topo_pot;
EXTERN_PARAMETERS double th_top;

EXTERN_PARAMETERS double ch_pot;
const int ch_pot_dir=0;
const int ch_pot_n=0;
EXTERN_PARAMETERS int use_charge_pot;

EXTERN_PARAMETERS int compute_corr_each;

#endif
