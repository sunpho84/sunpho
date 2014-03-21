#include <iostream>

#include "debug.hpp"
#include "geometry.hpp"
#include "simul.hpp"
#include "per_site_neighs.hpp"
#include "total_neighs.hpp"

using namespace std;

//internal main
void in_main(int narg,char **arg)
{
  GET_THREAD_ID();
  
  //test allocate and deallocate  
  double *v=NEW_COMMON("v") double[20];
  DELETE_COMMON(v);


  geometry_t *geometry=NEW_COMMON("geometry") geometry_t(2,10);

  //define first neighbors
  per_site_neighs_t *first_neighbors_per_site=NEW_COMMON("first_neigh_per_site") per_site_neighs_t;

  if(IS_MASTER_THREAD)
    {
      coords_t site(geometry->ndims);
      for(size_t dim=0;dim<geometry->ndims;dim++) site[dim]=0;
      for(int du=-1;du<=1;du+=2)
	for(size_t dim=0;dim<geometry->ndims;dim++)
	  {
	    site[dim]=du;
	    first_neighbors_per_site->add_neighbor(site);
	    site[dim]=0;
	  }
    }
  THREAD_BARRIER();

  
  total_neighs_t *first_neighbors=NEW_COMMON("first_neighbors") total_neighs_t(geometry,first_neighbors_per_site);
  
  if(IS_MASTER_THREAD)
    {
      //
    }
  
  DELETE_COMMON(first_neighbors_per_site);
  DELETE_COMMON(first_neighbors);
  DELETE_COMMON(geometry);
  
  THREAD_BARRIER();
}

int main(int narg,char**arg)
{
  simul=new simul_t(narg,arg,in_main);
  
  return 0;
}
