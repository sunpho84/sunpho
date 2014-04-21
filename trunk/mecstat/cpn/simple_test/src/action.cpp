#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "parameters.hpp"
#include "staples.hpp"
#include "types.hpp"

//compute the total energy or action
double energy()
{
  double res=0;
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      {
	int site_up=neighup(site,mu);
	for(int n=0;n<N;n++) res+=(conj(zeta(site_up)[n])*zeta(site)[n]*lambda(site)[mu]).real();
	}
  
  return -(2*res-2*V*NDIMS);
}
double action()
{return energy()/g;}

//return the energy/action of a single site
//NB: the total action will be half the sum of the energy of all sites!
double site_energy(int site)
{
  double res=0;
  for(int mu=0;mu<NDIMS;mu++)
    {
      int site_up=neighup(site,mu);
      for(int n=0;n<N;n++) res+=(conj(zeta(site_up)[n])*zeta(site)[n]*lambda(site)[mu]).real();
      int site_dw=neighdw(site,mu);
      for(int n=0;n<N;n++) res+=(conj(zeta(site_dw)[n])*zeta(site)[n]*conj(lambda(site_dw)[mu])).real();
    }
  
  return -(2*res-4*NDIMS);
}
double site_action(int site)
{return site_energy(site)/g;}

//return the energy/action of a single link
double link_energy(int site,int mu)
{
  double res=0;
  
  int site_up=neighup(site,mu);
  for(int n=0;n<N;n++) res+=(conj(zeta(site_up)[n])*zeta(site)[n]*lambda(site)[mu]).real();
  
  return -(2*res-2);
}
double link_action(int site,int mu)
{return link_energy(site,mu)/g;}

double site_staple_energy(int site,dcomplex *staple)
{
  double res=0;
  for(int n=0;n<N;n++) res+=(staple[n]*conj(zeta(site)[n])).real();
  return -(res-4*NDIMS);
}
double site_staple_action(int site,dcomplex *staple)
{return site_staple_energy(site,staple)/g;}

//energy of a link computed from staples
double link_staple_energy(int site,int mu)
{
  dcomplex staple;
  link_staple(staple,site,mu);
  return -((staple*conj(lambda(site)[mu])).real()-2);
}
double link_staple_action(int site,int mu)
{return link_staple_energy(site,mu)/g;}

