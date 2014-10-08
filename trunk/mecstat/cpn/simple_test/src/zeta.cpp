#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "parameters.hpp"
#include "types.hpp"

#include <iostream>

using namespace std;

//return a scalar product between two zetas
dcomplex get_zeta_compl_scalprod(dcomplex *a,dcomplex *b)
{
  dcomplex res=0;
  for(int n=0;n<N;n++) res+=conj(a[n])*b[n];
  return res;
}

//only real part
double get_zeta_real_scalprod(dcomplex *a,dcomplex *b)
{
  double res=0;
  for(int n=0;n<N;n++) res+=a[n].real()*b[n].real()+a[n].imag()*b[n].imag();
  return res;
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

//return the norm2 of a zeta
double get_zeta_norm2(dcomplex *z)
{
  double norm2=0;
  for(int n=0;n<N;n++) norm2+=norm(z[n]);
  return norm2;
}

//orthogonalize
void zeta_orthogonalize_with(dcomplex *z,dcomplex *w)
{
  double norm_with=get_zeta_real_scalprod(w,z)/get_zeta_norm2(w);
  for(int n=0;n<N;n++) z[n]-=norm_with*w[n];
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
