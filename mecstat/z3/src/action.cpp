#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define EXTERN_ACTION

#include "action.hpp"
#include "data.hpp"
#include "geometry.hpp"
#include "macros.hpp"
#include "parameters.hpp"

//compute N
int compute_N(z3_t *phi)
{
  int N=0;
  for(int s=0;s<V;s++)
    {
      //get value on s
      int sv=phi[s];
      N+=contr_N[sv];
    }
  
  return N;
}

//compute N0
int compute_N0(z3_t *phi)
{
  int N0=0;
  for(int s=0;s<V;s++)
    {
      //get value on s
      int sv=phi[s];
      N0+=contr_N0[sv];
    }
  
  return N0;
}

//compute the energy
double compute_energy_internal()
{
  //compute hopping part of the action
  int ndiff=V-glb_nequals;
  double E=contr_act_equals*glb_nequals+contr_act_diff*ndiff;
  
  return E;
}
  
//use global info
double compute_action_internal()
{
  //compute hopping part of the action
  double E=compute_energy_internal();
  
  //compute time part
  double M=0;
  for(int i=0;i<3;i++) M+=contr_act_re[i]*glb_ntypes[i];
  
  return -tau*(E+kappa*M);
} 

//return the action
double compute_action(z3_t *phi)
{
  //reset
  glb_nequals=0;
  for(int d=0;d<3;d++) glb_ntypes[d]=0;
  
  //compute ingredients
  for(int s=0;s<V;s++)
    {
      //get site value
      z3_t sv=phi[s];
      
      //count each type
      glb_ntypes[sv]++;
      
      //scan all neighbors
      for(int mu=0;mu<NDIMS;mu++)
	{
	  //take neighbors value
	  int sup=neighup(s,mu);
	  z3_t supv=phi[sup];
	  
	  //if they are equal, add it
	  glb_nequals+=(supv==sv);
	}
    }
  
  return compute_action_internal();
}
