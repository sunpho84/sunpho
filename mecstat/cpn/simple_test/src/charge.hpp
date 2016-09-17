#ifndef _CHARGE_HPP
#define _CHARGE_HPP

#include "parameters.hpp"
#include "types.hpp"

extern meta_pars_t chrono_charge;

dcomplex charge(dcomplex *z,dcomplex *l,int n=ch_pot_n,int mu=ch_pot_dir);
void sum_charge_lambda_force(dcomplex *z,dcomplex *l);
void sum_charge_zeta_force(dcomplex *z,dcomplex *l);

#endif

