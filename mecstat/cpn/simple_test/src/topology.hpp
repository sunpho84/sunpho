#ifndef _TOPOPOLOGY_HPP
#define _TOPOPOLOGY_HPP

extern meta_pars_t chrono_topo;

void compute_unsmeared_topological_force(double *fpi,dcomplex *l);
void sum_topological_force(double *fpi,double rho,int nlev,dcomplex *l);
double compute_theta_pot(dcomplex *l);
double geometric_topology(dcomplex *z);
double geometric_topology_simplified(dcomplex *z);
double topology(dcomplex *l);

#endif
