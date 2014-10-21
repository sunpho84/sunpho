#ifndef _ZETA_HPP
#define _ZETA_HPP

double check_zeta_unitarity(dcomplex *z);
void check_zeta_conf_unitarity(dcomplex *z);
dcomplex get_zeta_compl_scalprod(dcomplex *a,dcomplex *b);
double get_zeta_real_scalprod(dcomplex *a,dcomplex *b);
double get_zeta_norm(dcomplex *z);
double get_zeta_norm2(dcomplex *z);
void get_zeta_P(dcomplex *P,dcomplex *z);
void zeta_unitarize(dcomplex *z);
void zeta_orthogonalize_with(dcomplex *z,dcomplex *w);

#endif
