#include <iostream>

#include "threads.hpp"
#include "debug.hpp"
#include "field.hpp"
#include "geometry.hpp"
#include "neighs.hpp"
#include "per_site_neighs.hpp"
#include "simul.hpp"

#include "cp_n.hpp"

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
  GET_THREAD_ID();
  
  simul->init_glb_rnd_gen(101);
  
  O_n_t<3> r;
  U1_t c(3);
  if(IS_MASTER_THREAD) r.set_to_rnd(simul->glb_rnd_gen);
  //r.set_to_one();
  MASTER_PRINTF("%lg %lg\n",r[0].real(),r.get_norm());
  
  CP_N_t<3> u(10,3.0,CP_N_t<3>::HOT);
  u.Zeta.normalize();
  //u.Lambda.normalize();
  MASTER_PRINTF("Zeta norm: %lg, Lambda norm: %lg\n",u.Zeta.get_norm2(),u.Lambda.get_norm2());
  MASTER_PRINTF("Zeta norm: %lg, Lambda norm: %lg\n",u.Zeta.get_norm2(),u.Lambda.get_norm2());
  
  //MASTER_PRINTF("Energy: %d\n",compute_energy(t));
  
  THREAD_BARRIER();
}

int main(int narg,char**arg)
{
  simul=new simul_t(narg,arg,in_main);
  
  return 0;
}
