#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

//parameters
const int N=21;
//const int N=10;
const int nhmc_steps=15;
const double beta=0.7;
const int L=72;
//const int L=20;
const double g=1/(N*beta);
const double th_top=2;
const int use_topo_pot=1;

const double chrono_topo_coeff=0.1;
const double chrono_topo_width=0.1;

#endif
