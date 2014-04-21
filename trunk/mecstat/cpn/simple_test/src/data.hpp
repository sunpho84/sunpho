#ifndef _DATA_HPP
#define _DATA_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "parameters.hpp"
#include "types.hpp"

extern dcomplex *zeta_data,*lambda_data;

//return zeta and lambda
inline dcomplex *zeta(int site)
{return zeta_data+site*N;}
inline dcomplex *lambda(int site)
{return lambda_data+site*NDIMS;}

#endif
