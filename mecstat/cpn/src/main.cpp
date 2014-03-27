#include <iostream>

#include "debug.hpp"
#include "field.hpp"
#include "geometry.hpp"
#include "neighs.hpp"
#include "per_site_neighs.hpp"
#include "simul.hpp"

using namespace std;

int compute_energy(field_t<int> &t)
{
  GET_THREAD_ID();
  
  t.sync_outer_sites();

  int E=0;
  PARALLEL_FOR(iel,0,t.neighs_ptr->geometry->nloc_sites)
    {
      int s=t[iel];
      int a=0;
      for(size_t dir=0;dir<t.neighs_ptr->nneighs_per_site;dir++)
	{
	  a+=t[t.get_neigh(iel,dir)];
	  //MASTER_PRINTF("%d %d %d\n",iel,(int)dir,t[t.get_neigh(iel,dir)]);
	}
      E+=s*a;
    }
  
  return rank_threads_reduce(E);
}

//internal main
void in_main(int narg,char **arg)
{
  simul->init_glb_rnd_gen(101);
  
  GET_THREAD_ID();
  
  //test allocate and deallocate  
  double *v=NEW_BLOCKING("v") double[20];
  DELETE_BLOCKING(v);

  geometry_t *geometry=NEW_BLOCKING("geometry") geometry_t(2,10);
  geometry->init_loc_rnd_gen();
  
  //fill a field
  field_t<int> t("t",geometry->first_neighbors);
  PARALLEL_FOR_SITES_OF_FIELD(iel,t)
    t[iel]=geometry->loc_rnd_gen[iel].get_pm_one();
  THREAD_BARRIER();
  
  MASTER_PRINTF("Energy: %d\n",compute_energy(t));
  
  DELETE_BLOCKING(geometry);
  
  THREAD_BARRIER();
}

int main(int narg,char**arg)
{
  simul=new simul_t(narg,arg,in_main);
  
  return 0;
}
