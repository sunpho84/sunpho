#ifndef _ZETA_HPP
#define _ZETA_HPP

double check_zeta_unitarity(dcomplex *z);
void get_zeta_scalprod(dcomplex &res,dcomplex *a,dcomplex *b);
double get_zeta_scalprod(dcomplex *a,dcomplex *b);
double get_zeta_norm(dcomplex *z);
void get_zeta_P(dcomplex *P,dcomplex *z);
void zeta_unitarize(dcomplex *z);

#endif
