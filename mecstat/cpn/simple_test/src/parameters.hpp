#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef EXTERN_PARAMETERS
#define EXTERN_PARAMETERS extern
#endif

//parameters
const int nhmc_steps=18;
EXTERN_PARAMETERS int N;
EXTERN_PARAMETERS double beta;
EXTERN_PARAMETERS int L;
EXTERN_PARAMETERS double g;
const double th_top=0;
const int use_topo_pot=0;

EXTERN_PARAMETERS int nstout_lev;
EXTERN_PARAMETERS double stout_rho;

const int chrono_topo_after=300;
const double chrono_topo_coeff=0.2;
const double chrono_topo_width=0.3;
const double chrono_topo_barr=5;

#endif
