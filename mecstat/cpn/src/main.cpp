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
  double *v=ALLOCATE("v",20,double);
  FREE(v);

  //geometry_t *geometry=CAST_PTR_FROM_MASTER_THREAD(new geometry_t(2,10));
  geometry_t *geometry=CAST_PTR_FROM_MASTER_THREAD(new geometry_t(2,10));
  
  //define first neighbors
  per_site_neighs_t first_neighbors;

  coords_t site(geometry->ndims);
  for(size_t dim=0;dim<geometry->ndims;dim++) site[dim]=0;
  for(int du=-1;du<=1;du+=2)
    for(size_t dim=0;dim<geometry->ndims;dim++)
      {
	site[dim]=du;
	first_neighbors.add_neighbor(site);
	site[dim]=0;
      }
  
  total_neighs_t(geometry,first_neighbors);
  
  if(IS_MASTER_THREAD)
    {
      //
    }
  THREAD_BARRIER();
}

int main(int narg,char**arg)
{
  simul=new simul_t(narg,arg,in_main);
  
  return 0;
}
