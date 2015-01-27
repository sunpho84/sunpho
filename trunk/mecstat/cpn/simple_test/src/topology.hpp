#ifndef _TOPOPOLOGY_HPP
#define _TOPOPOLOGY_HPP

#include <vector>

extern std::vector<double> chrono_topo_past_values;
extern std::vector<double> chrono_topo_past_weight;
extern int ngrid;
extern double dx_grid;
extern vector<double> topo_grid,topo_grid_ave;

void draw_chrono_topo_potential(bool ave=false);
void draw_chrono_topo_force(int isweep);
void compute_unsmeared_topological_force(double *fpi,dcomplex *l);
void compute_topological_force(double *fpi,double rho,int nlev,dcomplex *l);
double compute_theta_pot(dcomplex *l,int isweep=0);
double compute_theta_pot(double Q,int isweep=0,bool ave=false);
double geometric_topology(dcomplex *z);
double geometric_topology_simplified(dcomplex *z);
double topology(dcomplex *l);
void update_chrono_potential(double Q,int isweep);
void finish_topological_force(double *f,dcomplex *,dcomplex *l);

#endif
