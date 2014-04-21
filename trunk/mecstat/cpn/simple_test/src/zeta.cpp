#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "parameters.hpp"
#include "types.hpp"

//return a scalar product between two zetas
void get_zeta_scalprod(dcomplex &res,dcomplex *a,dcomplex *b)
{
  res=0;
  for(int n=0;n<N;n++) res+=conj(a[n])*b[n];
}
double get_zeta_scalprod(dcomplex *a,dcomplex *b)
{
  dcomplex sc;
  get_zeta_scalprod(sc,a,b);
  return sc.real();
}

//return inner product of zeta
void get_zeta_P(dcomplex *P,dcomplex *z)
{
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      P[i*N+j]=conj(z[i])*z[j];
}

//return the norm of a zeta
double get_zeta_norm(dcomplex *z)
{
  double norm2=0;
  for(int n=0;n<N;n++) norm2+=norm(z[n]);
  return sqrt(norm2);
}

//reunitarize a zeta
void zeta_unitarize(dcomplex *z)
{
  double zeta_norm_inv=1/get_zeta_norm(z);
  for(int n=0;n<N;n++) z[n]*=zeta_norm_inv;
}

//return the deviation from unitarity of a zeta
double check_zeta_unitarity(dcomplex *z)
{return fabs(get_zeta_norm(z)-1);}
