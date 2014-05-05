#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "types.hpp"

//compute the staple of zeta
void site_staple(dcomplex *staple,dcomplex *z,dcomplex *l,int site)
{
  for(int n=0;n<N;n++) staple[n]=0;

  for(int mu=0;mu<NDIMS;mu++)
    {
      int site_up=neighup(site,mu);
      for(int n=0;n<N;n++) staple[n]+=(double)2.0*z[site_up*N+n]*conj(l[site*NDIMS+mu]);
      int site_dw=neighdw(site,mu);
      for(int n=0;n<N;n++) staple[n]+=(double)2.0*z[site_dw*N+n]*l[site_dw*NDIMS+mu];
    }
}

//compute the staple of lambda
void link_staple(dcomplex &staple,dcomplex *z,int site,int mu)
{
  staple=0;
  int site_up=neighup(site,mu);
  for(int n=0;n<N;n++) staple+=(double)2.0*conj(z[site*N+n])*z[site_up*N+n];
}

//compute the staples of the topological term
void topo_staple(dcomplex &staple,dcomplex *l,int s,int mu)
{
  int nu=!mu;
  staple=l[neighup(s,mu)*NDIMS+nu]*conj(l[neighup(s,nu)*NDIMS+mu]*l[s*NDIMS+nu])-
    conj(l[neighup(neighdw(s,nu),mu)*NDIMS+nu]*l[neighdw(s,nu)*NDIMS+mu])*l[neighdw(s,nu)*NDIMS+nu];
}
