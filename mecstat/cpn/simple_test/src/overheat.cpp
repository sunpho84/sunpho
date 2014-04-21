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
  site_staple(staple,site);
  double staple_norm=get_zeta_norm(staple);
  double staple_energy=get_zeta_scalprod(zeta(site),staple);
  
  //compute theta in the simple way
  double ctheta_old=staple_energy/staple_norm;
  double theta_old=acos(ctheta_old);
  
  //if theta is too small we switch to alternative method
  if(fabs(theta_old)<1.e-4) 
    {
      dcomplex diff[N];
      for(int n=0;n<N;n++) diff[n]=zeta(site)[n]-staple[n]/staple_norm;
      theta_old=asin(get_zeta_scalprod(diff,diff));
    }
  
  //it theta is too small the algorithm is undefined
  if(fabs(theta_old)>5.e-8)
    {
      //extract theta and compute its cos
      double theta_new=get_theta(beta*N*staple_norm,N);
      double ctheta_new=cos(theta_new);
      
      //extract remaining components
      double a=ctheta_new/staple_norm,b=ctheta_old/staple_norm,c=sin(theta_new)/sin(theta_old);
      for(int n=0;n<N;n++) zeta(site)[n]=a*staple[n]-(zeta(site)[n]-b*staple[n])*c;
      
      //reunitarize
      zeta_unitarize(zeta(site));
    }
  else cout<<"skipping site "<<site<<": "<<theta_old<<endl;
}

//update a link using overrelaxion/heatbath
void overheat_update_link(int site,int mu)
{
  //get the staple
  dcomplex staple;
  link_staple(staple,site,mu);

  //compute the staple norm and energy
  double staple_norm=sqrt(norm(staple));
  double staple_energy=get_lambda_scalprod(lambda(site)[mu],staple);
  
  //compute theta in the simple way
  double ctheta_old=staple_energy/staple_norm;
  double theta_old=acos(ctheta_old);
  
  //if theta is too small we switch to alternative method
  if(fabs(theta_old)<1.e-4) 
    {
      dcomplex diff=lambda(site)[mu]-staple/staple_norm;
      theta_old=asin(get_lambda_scalprod(diff,diff));
    }
  
  //it theta is too small the algorithm is undefined
  if(fabs(theta_old)>5.e-8)
    {
      //extract theta and compute its cos
      double theta_new=get_theta_1(beta*N*staple_norm);
      double ctheta_new=cos(theta_new);
      
      //extract remaining components  
      double a=ctheta_new/staple_norm,b=ctheta_old/staple_norm,c=sin(theta_new)/sin(theta_old);
      lambda(site)[mu]=a*staple-(lambda(site)[mu]-b*staple)*c;
      
      //reunitarize
      lambda_unitarize(lambda(site)[mu]);
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
