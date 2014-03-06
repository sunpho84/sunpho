#ifndef _GEOMETRY_HPP
#define _GEOMETRY_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>

#define NMU 2
#define BW 0
#define FW 1

typedef int coords_t[NMU];
typedef int neighs_t[2*NMU];

//structure to hold geometry
struct geometry_t
{
  MPI_Comm rank_comm;      //communicator
  int nranks;              //total number of ranks
  coords_t nranks_per_dir; //ranks per direction
  coords_t rank_coords;    //coordinates in the proc grid
  neighs_t neigh_ranks;    //neighbors ranks
  
  coords_t glb_sizes,loc_sizes; //local dimensions
  int glb_size,loc_size;        //number of sites
  
  void loc_coords_of_loc_site(coords_t loc_coords,int loc_site){coords_of_site(loc_coords,loc_site,loc_sizes);}
  int loc_site_of_loc_coords(coords_t loc_coords){return site_of_coords(loc_coords,loc_sizes);}
  void glb_coords_of_glb_site(coords_t glb_coords,int glb_site){coords_of_site(glb_coords,glb_site,glb_sizes);}
  int glb_site_of_glb_coords(coords_t glb_coords){return site_of_coords(glb_coords,glb_sizes);}
  void start(coords_t ext_glb_sizes);
  
  geometry_t(coords_t ext_glb_sizes){start(ext_glb_sizes);};
  geometry_t(int glb_comm_size)
  {
    coords_t ext_glb_sizes;
    for(int mu=0;mu<NMU;mu++) ext_glb_sizes[mu]=glb_comm_size;
    start(ext_glb_sizes);
  }
private:
  void coords_of_site(coords_t coords,int loc_site,coords_t sizes);
  int site_of_coords(coords_t coords,coords_t sizes);
  geometry_t();
};

#endif
