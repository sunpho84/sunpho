#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

const double tau=0.17;
const double kappa=0.05;
const int L0=24,L1=24;
const int L[3]={L0,L0,L1};

/*
const double tau=0.5490000000;
const double kappa=0;
const int L0=40,L1=20;
const int L[3]={L0,L0,L1};
*/

/*
const double tau=0.45;
const double kappa=-2/tau; //in facts this is h=H*tau
const int L0=30,L1=30;
const int L[3]={L0,L0,L1};
*/

const int ntraj=10000;
const int nterm=1000;
const double meta_sigma_N=10;
const double meta_sigma_N0=10;
const double meta_coeff=0;//1000;

#endif
