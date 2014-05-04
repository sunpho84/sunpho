#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cmath>

#include "types.hpp"

//return the result of the scalar product of two lambda
double get_lambda_real_scalprod(dcomplex a,dcomplex b)
{return (conj(a)*b).real();}

//return the norm of a lambda
double get_lambda_norm(dcomplex &l)
{return sqrt(norm(l));}

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
