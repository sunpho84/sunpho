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

#define N 2
#define beta 1.1
#define L 12

//internal main
void in_main(int narg,char **arg)
{
  GET_THREAD_ID();
  
  simul->init_glb_rnd_gen(101);
  
  //initialize a CPN simulation with N, as HOT
  CP_N_t<N> u(L,beta,CP_N_t<N>::HOT);
  //u.Zeta.normalize();
  //u.Lambda.normalize();
  //MASTER_PRINTF("Zeta norm: %lg, Lambda norm: %lg\n",u.Zeta.get_norm2(),u.Lambda.get_norm2());
  
  for(int iter=0;iter<1000;iter++)
    {
      for(int i=0;i<u.geometry->nglb_sites;i++)
	{
	  int rank,loc;
	  u.geometry->rank_and_loc_site_of_glb_coords(rank,loc,u.geometry->glb_coords_of_glb_site(i));
	  
	  //change Zeta
	  {	 
	    double ori_ac=u.get_action();
	    O_n_t<N> ori_val;
	    if(simul->rank_id==rank && IS_MASTER_THREAD)
	      {
		ori_val=u.Zeta[loc];
		u.Zeta[loc].set_to_rnd(u.geometry->loc_rnd_gen[loc]);
	      }
	    u.Zeta.mark_touched();
	    double new_ac=u.get_action();
	    double diff=new_ac-ori_ac;

	    double p=exp(-diff);
	    double e;
	    MASTER_THREAD_BROADCAST(e,simul->glb_rnd_gen.get_unif(0,1));
	    bool acc=(e<p);
	    if(!acc)
	      {
		if(simul->rank_id==rank && IS_MASTER_THREAD) u.Zeta[loc]=ori_val;
		u.Zeta.mark_touched();
	      }
	    
	    //change Lambda
	    for(int mu=0;mu<2;mu++)
	      {
		//change Zeta
		double ori_ac=u.get_action();
		U1_t ori_val;
		if(simul->rank_id==rank && IS_MASTER_THREAD)
		  {
		    ori_val=u.Lambda[loc][mu];
		    u.Lambda[loc][mu].set_to_rnd(u.geometry->loc_rnd_gen[loc]);
		  }
		u.Lambda.mark_touched();
		double new_ac=u.get_action();
		double diff=new_ac-ori_ac;
		double p=exp(-diff);
		double e;
		MASTER_THREAD_BROADCAST(e,simul->glb_rnd_gen.get_unif(0,1));
		bool acc=(e<p);
		if(!acc)
		  {
		    if(simul->rank_id==rank && IS_MASTER_THREAD) u.Lambda[loc][mu]=ori_val;
		    u.Lambda.mark_touched();
		  }
		
		//MASTER_PRINTF("rank %d, thread %d, Energy: %lg %lg %lg %lg, %lg %lg, %d\n",simul->rank_id,thread_id,ori_en,new_en,diff,exp(-diff),p,e,acc);
		THREAD_BARRIER();
	      }
	  }
	}
      MASTER_PRINTF("Energy %lg\n",u.get_energy()/u.geometry->nglb_sites/N);
    }
}

int main(int narg,char**arg)
{
  simul=new simul_t(narg,arg,in_main);
  
  return 0;
}
