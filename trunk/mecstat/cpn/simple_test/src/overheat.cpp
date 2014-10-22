#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cmath>
#include <iostream>

#include "data.hpp"
#include "geometry.hpp"
#include "lambda.hpp"
#include "random.hpp"
#include "staples.hpp"
#include "zeta.hpp"

using namespace std;

//update a site using overrelaxion/heatbath
void overheat_update_site(int site)
{
  //get the staple, its norm and energy
  dcomplex staple[N];
  site_staple(staple,zeta,lambda,site);
  double staple_norm=get_zeta_norm(staple);
  
  //if the staple is 0 we can set zeta to random
  if(staple_norm==0) set_ON_to_rnd(zeta+site*N);
  else
    {
      //generate the new theta according to the microscopic distriution of prob
      double theta_new=get_theta(beta*N*staple_norm,N);
      
      //draw a random vetor orthogonal to the staple until it's not too small
      dcomplex R[N];
      double R_norm;
      do
	{
	  set_ON_to_rnd(R);
	  zeta_orthogonalize_with(R,staple);
	  R_norm=get_zeta_norm(R);
	}
      while(R_norm<1e-3);
      
      //set the components parallel and perpendicular to the staple
      double par_comp=cos(theta_new)/staple_norm;
      double perp_comp=sin(theta_new)/R_norm;
      
      //compose the new vector
      for(int n=0;n<N;n++) zeta[site*N+n]=staple[n]*par_comp+R[n]*perp_comp;
    }
}

//update a link using overrelaxion/heatbath
void overheat_update_link(int site,int mu)
{
  //get the staple
  dcomplex staple;
  link_staple(staple,zeta,site,mu);

  //compute the staple norm
  double staple_norm=sqrt(norm(staple));
  
  //if the staple is 0 we can set zeta to random
  if(staple_norm==0) set_U1_to_rnd(lambda[site*NDIMS]);
  else
    {
      //generate the new theta according to the microscopic distriution of prob
      double theta_new=get_theta_1(beta*N*staple_norm);
      
      //draw a vector orthogonal to the staple until it is not too small
      dcomplex R;
      double R_norm;
      do
	{
	  set_U1_to_rnd(R);
	  lambda_orthogonalize_with(R,staple);
	  R_norm=get_lambda_norm(R);
	}
      while(R_norm<1e-3);
      
      //set the components parallel and perpendicular to the staple
      double par_comp=cos(theta_new)/staple_norm;
      double perp_comp=sin(theta_new)/R_norm;
      lambda[site*NDIMS+mu]=staple*par_comp+R*perp_comp;
    }
}

//sweep all the lattice with overrelaxation/heatbath
void overheat_sweep()
{
  //loop over sites
  for(int site=0;site<V;site++)
    {
      overheat_update_site(site);
      for(int mu=0;mu<NDIMS;mu++) overheat_update_link(site,mu);
    }
}
