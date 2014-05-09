#ifndef _TOPOPOLOGY_HPP
#define _TOPOPOLOGY_HPP

#include <vector>

extern std::vector<double> chrono_topo_past_values;

void draw_chrono_topo_potential();
void compute_unsmeared_topological_force(double *fpi,dcomplex *l);
void compute_topological_force(double *fpi,double rho,int nlev,dcomplex *l);
double compute_theta_pot(dcomplex *l);
double geometric_topology(dcomplex *z);
double geometric_topology_simplified(dcomplex *z);
double topology(dcomplex *l);
void finish_topological_force(double *f,dcomplex *,dcomplex *l);

#endif
