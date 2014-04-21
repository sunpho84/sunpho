#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "types.hpp"

//compute the staple of zeta
void site_staple(dcomplex *staple,int site)
{
  for(int n=0;n<N;n++) staple[n]=0;

  for(int mu=0;mu<NDIMS;mu++)
    {
      int site_up=neighup(site,mu);
      for(int n=0;n<N;n++) staple[n]+=(double)2.0*zeta(site_up)[n]*conj(lambda(site)[mu]);
      int site_dw=neighdw(site,mu);
      for(int n=0;n<N;n++) staple[n]+=(double)2.0*zeta(site_dw)[n]*lambda(site_dw)[mu];
    }
}

//compute the staple of lambda
void link_staple(dcomplex &staple,int site,int mu)
{
  staple=0;
  int site_up=neighup(site,mu);
  for(int n=0;n<N;n++) staple+=(double)2.0*conj(zeta(site)[n])*zeta(site_up)[n];
}
