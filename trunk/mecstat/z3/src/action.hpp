#ifndef _ACTION_HPP
#define _ACTION_HPP

#include "parameters.hpp"
#include "types.hpp"

#ifndef EXTERN_ACTION
 #define EXTERN_ACTION extern
#endif

#include <vector>

EXTERN_ACTION std::vector<N_t> data_N;
EXTERN_ACTION int contr_N[3];
EXTERN_ACTION int contr_N0[3];
EXTERN_ACTION int glb_N,glb_N0;
EXTERN_ACTION int glb_ntypes[3],glb_nequals;
EXTERN_ACTION double contr_act_re[3];
EXTERN_ACTION double contr_act_equals,contr_act_diff;

double compute_action_internal();
double compute_energy_internal();
double compute_magnetization_internal();
double compute_action(z3_t *phi);
int compute_N(z3_t *phi);
int compute_N0(z3_t *phi);

#endif
