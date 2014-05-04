#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "parameters.hpp"
#include "staples.hpp"
#include "topology.hpp"
#include "types.hpp"

//compute the total energy or action
double energy(dcomplex *z,dcomplex *l)
{
  double res=0;
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      {
	int site_up=neighup(site,mu);
	for(int n=0;n<N;n++) res+=(conj(z[site_up*N+n])*z[site*N+n]*l[site*NDIMS+mu]).real();
      }
  
  return -(2*res-2*V*NDIMS);
}
double action(dcomplex *z,dcomplex *l)
{return energy(z,l)/g;}

//return the topological action
double topo_action(dcomplex *l)
{return topology(l)*th_top;}

//return the energy/action of a single site
//NB: the total action will be half the sum of the energy of all sites!
double site_energy(dcomplex *z,dcomplex *l,int site)
{
  double res=0;
  for(int mu=0;mu<NDIMS;mu++)
    {
      int site_up=neighup(site,mu);
      for(int n=0;n<N;n++) res+=(conj(z[site_up*N+n])*z[site*N+n]*l[site*NDIMS+mu]).real();
      int site_dw=neighdw(site,mu);
      for(int n=0;n<N;n++) res+=(conj(z[site_dw*N+n])*z[site*N+n]*conj(l[site_dw*NDIMS+mu])).real();
    }
  
  return -(2*res-4*NDIMS);
}
double site_action(dcomplex *z,dcomplex *l,int site)
{return site_energy(z,l,site)/g;}

//return the energy/action of a single link
double link_energy(dcomplex *z,dcomplex *l,int site,int mu)
{
  double res=0;
  
  int site_up=neighup(site,mu);
  for(int n=0;n<N;n++) res+=(conj(z[site_up*N+n])*z[site*N+n]*l[site*NDIMS+mu]).real();
  
  return -(2*res-2);
}
double link_action(dcomplex *z,dcomplex *l,int site,int mu)
{return link_energy(z,l,site,mu)/g;}

double site_staple_energy(int site,dcomplex *staple,dcomplex *z)
{
  double res=0;
  for(int n=0;n<N;n++) res+=(staple[n]*conj(z[site*N+n])).real();
  return -(res-4*NDIMS);
}
double site_staple_action(int site,dcomplex *staple,dcomplex *z)
{return site_staple_energy(site,staple,z)/g;}

//energy of a link computed from staples
double link_staple_energy(dcomplex *z,dcomplex *l,int site,int mu)
{
  dcomplex staple;
  link_staple(staple,z,site,mu);
  return -((staple*conj(l[site*NDIMS+mu])).real()-2);
}
double link_staple_action(dcomplex *z,dcomplex *l,int site,int mu)
{return link_staple_energy(z,l,site,mu)/g;}

