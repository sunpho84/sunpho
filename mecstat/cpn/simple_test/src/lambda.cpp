#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cmath>
#include <iostream>

#include "geometry.hpp"
#include "types.hpp"

using namespace std;

//return the result of the scalar product of two lambda
double get_lambda_real_scalprod(dcomplex a,dcomplex b)
{return (conj(a)*b).real();}

//return the norm of a lambda
double get_lambda_norm(dcomplex &l)
{return sqrt(norm(l));}

//return the norm of a lambda
double get_lambda_norm2(dcomplex &l)
{return norm(l);}

//reunitarize a lambda
double lambda_unitarize(dcomplex &l)
{
  double n=get_lambda_norm(l);
  l*=1/n;
  return n;
}

//return the deviation from unitarity of a lambda
double check_lambda_unitarity(dcomplex &l)
{return fabs(get_lambda_norm(l)-1);}

//orthogonalize
void lambda_orthogonalize_with(dcomplex &l,dcomplex w)
{
  double norm_with=get_lambda_real_scalprod(w,l)/get_lambda_norm2(w);
  l-=norm_with*w;
}

//unitarity of the whole conf
void check_lambda_conf_unitarity(dcomplex *l)
{
  double lambda_nonun=0;
#pragma omp parallel for reduction(+:lambda_nonun)
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      lambda_nonun+=get_lambda_norm(l[site*NDIMS+mu])-1;
  
  cout<<"Lambda nonun: "<<lambda_nonun/sqrt(V*NDIMS)<<endl;
}
