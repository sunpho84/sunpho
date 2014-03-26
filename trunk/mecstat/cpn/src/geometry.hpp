#ifndef _GEOMETRY_HPP
#define _GEOMETRY_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <algorithm>
#include <functional>
#include <mpi.h>
#include <vector>

#include "debug.hpp"
#include "random.hpp"

class geometry_t;
class per_site_neighs_t;

//coordinate of a point
class coords_t : public std::vector<int>
{
public:
  coords_t(size_t ndims=0) : std::vector<int>(ndims) {}
  
  size_t ndims(){return size();}
  void check_compatibility(const coords_t &a,const coords_t &b)
  {if(a.size()!=b.size()) CRASH_SOFTLY("incompatible vectors");}
  int total_product(){int p=1;for(std::vector<int>::iterator it=begin();it!=end();it++)p*=*it;return p;}
  
  coords_t operator+(coords_t a)
  {coords_t b(a.ndims());transform(this->begin(),this->end(),a.begin(),b.begin(),std::plus<int>());return b;}
  coords_t operator-(coords_t a)
  {coords_t b(a.ndims());transform(this->begin(),this->end(),a.begin(),b.begin(),std::minus<int>());return b;}
  coords_t operator*(coords_t a)
  {coords_t b(a.ndims());transform(this->begin(),this->end(),a.begin(),b.begin(),std::multiplies<int>());return b;}
  coords_t operator/(coords_t a)
  {coords_t b(a.ndims());transform(this->begin(),this->end(),a.begin(),b.begin(),std::divides<int>());return b;}
  coords_t operator%(coords_t a)
  {coords_t b(a.ndims());transform(this->begin(),this->end(),a.begin(),b.begin(),std::modulus<int>());return b;}
};

#include "per_site_neighs.hpp"
#include "neighs.hpp"

//structure to hold geometry
struct geometry_t
{
  size_t ndims;
  int cart_rank;
  MPI_Comm cart_comm;        //communicator
  coords_t nranks_per_dir;   //ranks per direction
  coords_t rank_coords;      //coordinates in the proc grid
  
  bool homogeneous_partitioning;         //store if each node has the same number of sites
  coords_t glb_sizes,loc_sizes;          //local dimensions
  coords_t loc_comm_sizes;               //dimensions an all ranks apart last one
  int nglb_sites,nloc_sites;             //number of sites
  
  coords_t glb_coords_of_loc_origin;                     //global coords of local origin
  std::vector<coords_t> loc_coords_of_loc_site_table;    //lookup-table for loc site loc coords
  std::vector<coords_t> glb_coords_of_loc_site_table;    //lookup-table for loc site glb coords
  std::vector<int> glb_site_of_loc_site_table;           //lookup-table for loc->glb conversion
  
  bool loc_rnd_gen_inited;                    //flag for remembering if local random generator inited
  std::vector<rnd_gen_t> loc_rnd_gen;         //local random generators
  
  neighs_t *first_neighbors;      //connection at first neighbors
  
  void print();
  int bulk_volume(coords_t L);
  int bulk_recip_lat_volume(coords_t R,coords_t L);
  int compute_border_variance(coords_t L,coords_t P,int factorize_processor);

  coords_t glb_coords_of_origin_of_rank_coords(coords_t in_coords){return std::min(loc_comm_sizes*in_coords,glb_sizes);}
  coords_t loc_sizes_of_rank_coords(coords_t in_coords);
  void rank_and_loc_site_of_rel_glb_coords(int &rank,int &site,coords_t glb_coords)
  {rank_and_loc_site_of_glb_coords(rank,site,(glb_coords+glb_sizes)%glb_sizes);}
  void rank_and_loc_site_of_glb_coords(int &rank,int &site,coords_t coords);

  coords_t loc_coords_of_loc_site(int loc_site){return loc_coords_of_loc_site_table[loc_site];}
  int loc_site_of_loc_coords(coords_t loc_coords){return site_of_coords(loc_coords,loc_sizes);}
  int loc_site_of_glb_coords(coords_t glb_coords){return site_of_coords(glb_coords-glb_coords_of_loc_origin,loc_sizes);}
  int glb_site_of_loc_coords(coords_t loc_coords){return glb_site_of_glb_coords(loc_coords+glb_coords_of_loc_origin);}
  int glb_site_of_loc_rel_coords(coords_t loc_coords)
  {return glb_site_of_glb_rel_coords(loc_coords+glb_coords_of_loc_origin);}
  coords_t glb_coords_of_glb_site(int glb_site){return coords_of_site(glb_site,glb_sizes);}
  coords_t glb_coords_of_loc_site(int loc_site){return glb_coords_of_loc_site_table[loc_site];}
  int glb_site_of_glb_coords(coords_t glb_coords){return site_of_coords(glb_coords,glb_sizes);}
  int glb_site_of_loc_site(int loc_site){return glb_site_of_loc_site_table[loc_site];}
  int glb_site_of_glb_rel_coords(coords_t glb_coords)
  {return site_of_coords((glb_coords+glb_sizes)%glb_sizes,glb_sizes);}
  void start(coords_t ext_glb_sizes);
  void init_loc_rnd_gen();
  
  ~geometry_t();
  geometry_t(coords_t ext_glb_sizes) : ndims(ext_glb_sizes.size()) {start(ext_glb_sizes);};
  geometry_t(size_t ndims,int glb_comm_size) : ndims(ndims)
  {
    coords_t ext_glb_sizes(ndims);
    for(size_t dim=0;dim<ndims;dim++) ext_glb_sizes[dim]=glb_comm_size;
    start(ext_glb_sizes);
  }
private:
  void partition_ranks_homogeneously();
  void partition_ranks_inhomogeneously();
  coords_t coords_of_site(int loc_site,coords_t sizes);
  int site_of_coords(coords_t coords,coords_t sizes);
  geometry_t();
};


#endif
