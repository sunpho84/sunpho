#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cmath>

#include "action.hpp"
#include "data.hpp"
#include "geometry.hpp"
#include "macros.hpp"
#include "parameters.hpp"
#include "random.hpp"

//update zeta with metropolis
void metro_update_site(int site)
{
  //copy zeta
  dcomplex ori[N];
  for(int n=0;n<N;n++) ori[n]=zeta[site*N+n];
  
  //change of action
  double ori_ac=site_action(zeta,lambda,site);
  set_ON_to_rnd(zeta+site*N,site);
  double fin_ac=site_action(zeta,lambda,site);
  
  //accept or not?
  double diff_ac=fin_ac-ori_ac;
  double t=exp(-diff_ac);
  double p=get_unif_double(1,site);
  if(p>t) for(int n=0;n<N;n++) zeta[site*N+n]=ori[n];
}

//update lambda with metropolis
void metro_update_link(int site,int mu)
{
  //copy lambdaa
  dcomplex ori=lambda[site*NDIMS+mu];
  
  //change of action
  double ori_ac=link_action(zeta,lambda,site,mu);
  set_U1_to_rnd(lambda[site*NDIMS+mu],site);
  double fin_ac=link_action(zeta,lambda,site,mu);
  
  //accept or not?
  double diff_ac=fin_ac-ori_ac;
  double t=exp(-diff_ac);
  double p=get_unif_double(1,site);
  if(p>t) lambda[site*NDIMS+mu]=ori;
}

//sweep all the lattice
void metro_sweep()
{
  //loop over sites
  for(int site=0;site<V;site++)
    {
      metro_update_site(site);
      for(int mu=0;mu<NDIMS;mu++) metro_update_link(site,mu);
    }
}
