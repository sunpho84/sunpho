#ifndef _STOUT_HPP
#define _STOUT_HPP

void stout_lambda(dcomplex *ext_dest,double rho,dcomplex *source);
void stout_lambda_whole_stack(dcomplex **out,double rho,int nstout_lev,dcomplex *in);
void stout_remap_force(dcomplex *&f_out,dcomplex *&f_in,double rho,dcomplex *l_unsm,dcomplex *l_sm);
void stout_remap_force(dcomplex *&f,double rho,int nlevls,dcomplex **l);

#endif
