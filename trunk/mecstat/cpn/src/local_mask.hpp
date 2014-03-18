#ifndef _LOCAL_MASK_HPP
#define _LOCAL_MASK_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <vector>

#include "all_to_all_comm.hpp"
#include "geometry.hpp"

//wrapper
struct coords_t_wrap
{
  coords_t data;
  coords_t_wrap(coords_t in){for(int mu=0;mu<NMU;mu++) data[mu]=in[mu];}
  int operator[](int i){return data[i];}
};

//structure to define required neighbors
struct req_neighs_per_site_t: std::vector<coords_t_wrap>
{void add_neighbor(coords_t disp){push_back(disp);}};

//structure to hold the local view of a field
//and movements on the lattice
struct local_mask_t
{
  int nsites;
  neighs_t *neighs;  
};

#endif
