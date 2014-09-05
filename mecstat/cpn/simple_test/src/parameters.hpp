#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef EXTERN_PARAMETERS
#define EXTERN_PARAMETERS extern
#endif

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
EXTERN_PARAMETERS int chrono_topo_after;//=300;
EXTERN_PARAMETERS double chrono_topo_coeff;//=0.2;
EXTERN_PARAMETERS double chrono_topo_width;//=0.3;
EXTERN_PARAMETERS double chrono_topo_barr;//=5;
EXTERN_PARAMETERS double chrono_topo_force_out;
EXTERN_PARAMETERS double chrono_topo_well_tempering;
EXTERN_PARAMETERS double chrono_topo_bend;
EXTERN_PARAMETERS int compute_corr_each;

#endif
