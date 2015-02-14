#ifndef _LINALGS_H
#define _LINALGS_H

#include "new_types/new_types_definitions.hpp"

namespace bissa
{
  void double_vector_glb_scalar_prod(double *res,double *a,double *b,int n);
  void double_vector_glb_collapse(double *res,double *a,int n);
  void double_vector_copy(double *a,double *b,int n);
  void double_vector_init_to_zero(double *a,int n);
  void double_vector_linear_comb(double *a,double *b,double c,double *d,double e,int n,int OPT=0);
  void double_vector_prod_double(double *out,double *in,double r,int n);
  void double_vector_normalize(double *rat,double *out,double *in,double fact,int n);
  void double_vector_prod_the_summ_double(double *out,double r,double *in1,double *in2,int n);
  void double_vector_summassign(double *out,double *in,int n);
  void double_vector_subt(double *out,double *in1,double *i2,int n);
  void double_vector_summ_double_vector_prod_double(double *a,double *b,double *c,double d,int n,int OPT=0);
  void single_vector_summ_single_vector_prod_single(float *a,float *b,float *c,float d,int n,int OPT=0);
  void parallel_memcpy(void *out,void *in,int n);
}

#endif
