#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "geometry.hpp"
#include "types.hpp"

#define EXTERN_DATA
#include "data.hpp"

//copy all the zeta
void copy_zeta_conf(dcomplex *dest,dcomplex *source)
{
#pragma omp parallel for
  for(int site=0;site<V;site++)
    for(int n=0;n<N;n++)
      dest[site*N+n]=source[site*N+n];
}

//copy all the lambda
void copy_lambda_conf(dcomplex *dest,dcomplex *source)
{
#pragma omp parallel for
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      dest[site*NDIMS+mu]=source[site*NDIMS+mu];
}
