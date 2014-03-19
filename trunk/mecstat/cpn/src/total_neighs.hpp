#ifndef _TOTAL_NEIGHS_HPP
#define _TOTAL_NEIGHS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <map>
#include <vector>

#include "geometry.hpp"
#include "per_site_neighs.hpp"

//hold the site per each rank
typedef std::map<int,std::map<int,int> > site_list_per_rank_t;

//holds temporarily info on neighs
typedef std::map<int,std::map<int,std::vector<int> > > temp_total_neighs_t;

//hold info on ranks asking id, size ans position
class rank_to_ask_t
{
public:
  int rank;
  int size;
  int dest;
};

//hold info on ranks asking id, size ans position
class rank_asking_t
{
public:
  int rank;
  int size;
  int dest;
  int *list_from;
};

//holds the neighbors
class total_neighs_t
{
  total_neighs_t();
public:
  geometry_t *geometry;    //geometry used to create

  size_t nneighs_per_site; //neighbors per site
  size_t nouter_sites;     //total number of sites non local
  size_t nsites_to_send;   //number of sites to send away
  size_t ntotal_sites;     //total number of sites, including outer ones
  int *neighs;             //movements
  
  int nranks_to_ask;    //number of ranks to ask to
  int nranks_asking;    //number of ranks asking data
  
  rank_to_ask_t *ranks_to_ask;    //info on the ranks to ask
  rank_asking_t *ranks_asking;    //info on asking ranks
  
  int *operator[](int isite){return neighs+isite*nneighs_per_site;}
  void mark_all_neighbors(int loc_site,coords_t glb_site_coords,per_site_neighs_t &per_site_neighs,
			  site_list_per_rank_t &outer_sites_per_rank,bool recursive=false);

  total_neighs_t(geometry_t *geometry,per_site_neighs_t &per_site_neighs);
  ~total_neighs_t();
};

#endif
