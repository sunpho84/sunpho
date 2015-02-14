#ifndef _TOPOPOLOGY_HPP
#define _TOPOPOLOGY_HPP

#include <vector>

extern int ngrid;
extern vector<double> topo_grid;

void load_chrono_topo_potential();
void draw_chrono_topo_potential();
void draw_chrono_topo_force();
void compute_unsmeared_topological_force(double *fpi,dcomplex *l);
void compute_topological_force(double *fpi,double rho,int nlev,dcomplex *l);
double compute_theta_pot(dcomplex *l);
double compute_theta_pot(double Q);
double geometric_topology(dcomplex *z);
double geometric_topology_simplified(dcomplex *z);
double topology(dcomplex *l);
void update_chrono_potential(double Q);
void finish_topological_force(double *f,dcomplex *,dcomplex *l);

#endif
