#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "lambda.hpp"
#include "geometry.hpp"
#include "routines.hpp"
#include "staples.hpp"
#include "zeta.hpp"

//update a site using microcanonical
void micro_update_site(int site)
{
  //get the staple, its norm and energy
  dcomplex staple[N];
  site_staple(staple,zeta,lambda,site);
  double staple_norm=get_zeta_norm(staple);
  double staple_energy=get_zeta_real_scalprod(zeta+site*N,staple);
  
  //extract site
  for(int n=0;n<N;n++) zeta[site*N+n]=2*staple_energy/sqr(staple_norm)*staple[n]-zeta[site*N+n];

  //reunitarize
  zeta_unitarize(zeta+site*N);
}

//update a link using microcanonical
void micro_update_link(int site,int mu)
{
  //get the staple
  dcomplex staple;
  link_staple(staple,zeta,site,mu);

  //compute the staple norm and energy
  double staple_norm=sqrt(norm(staple));
  double staple_energy=get_lambda_real_scalprod(lambda[site*NDIMS+mu],staple);
  
  //extract link
  lambda[site*NDIMS+mu]=2*staple_energy/sqr(staple_norm)*staple-lambda[site*NDIMS+mu];

  //reunitarize
  lambda_unitarize(lambda[site*NDIMS+mu]);
}

//sweep all the lattice with microcanonical
void micro_sweep()
{
  //loop over sites
  for(int site=0;site<V;site++)
    {
      micro_update_site(site);
      for(int mu=0;mu<NDIMS;mu++) micro_update_link(site,mu);
    }
}
