#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

const double tau=0.17;
const double kappa=0.05;
const int L0=24,L1=24;
const int L[3]={L0,L0,L1};

const int ntraj=100000;
const int nterm=1000;
const double meta_sigma_N=10;
const double meta_sigma_N0=10;
const double meta_coeff=0;//1000;

#endif
