#ifndef _ACTION_HPP
#define _ACTION_HPP
#include "types.hpp"
double action(dcomplex *z,dcomplex *l);
double topo_action(double rho,int nlev,dcomplex *l);
double energy(dcomplex *z,dcomplex *l);
double link_action(dcomplex *z,dcomplex *l,int site,int mu);
double link_energy(dcomplex *z,dcomplex *l,int site,int mu);
double link_staple_action(dcomplex *z,dcomplex *l,int site,int mu);
double link_staple_energy(dcomplex *z,dcomplex *l,int site,int mu);
double site_action(dcomplex *z,dcomplex *l,int site);
double site_energy(dcomplex *z,dcomplex *l,int site);
double site_staple_action(int site,dcomplex *staple,dcomplex *z);
double site_staple_energy(int site,dcomplex *staple,dcomplex *z);
double topo_action(dcomplex *l);
#endif
