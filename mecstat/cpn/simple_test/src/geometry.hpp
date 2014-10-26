#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifndef _GEOMETRY_HPP
#define _GEOMETRY_HPP

#ifndef EXTERN_GEOMETRY
 #define EXTERN_GEOMETRY extern
#endif

#include "macros.hpp"
#include "parameters.hpp"
#include "types.hpp"

//geometry
EXTERN_GEOMETRY int V,V_per_par;
EXTERN_GEOMETRY int *neigh_data;
EXTERN_GEOMETRY int npar;
EXTERN_GEOMETRY int *lx_of_par;

//return the neighbors
inline int &neighdw(int site,int mu)
{return neigh_data[0+2*(mu+NDIMS*site)];}
inline int &neighup(int site,int mu)
{return neigh_data[1+2*(mu+NDIMS*site)];}

int site_of_coords(coords c);
void coords_of_site(coords c,int site);

#endif
