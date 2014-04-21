#ifndef _ACTION_HPP
#define _ACTION_HPP

#include "types.hpp"

double energy();
double action();
double site_energy(int site);
double site_action(int site);
double link_energy(int site,int mu);
double link_action(int site,int mu);
double site_staple_energy(int site,dcomplex *staple);
double site_staple_action(int site,dcomplex *staple);
double link_staple_energy(int site,int mu);
double link_staple_action(int site,int mu);

#endif
