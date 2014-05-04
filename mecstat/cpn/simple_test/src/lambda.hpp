#ifndef _LAMBDA_HPP
#define _LAMBDA_HPP

#include "types.hpp"

double get_lambda_real_scalprod(dcomplex a,dcomplex b);
double get_lambda_norm(dcomplex &l);
double lambda_unitarize(dcomplex &l);
double check_lambda_unitarity(dcomplex &l);

#endif
