#ifndef _VECTORS_HPP
#define _VECTORS_HPP

#include <immintrin.h>
#include <stdio.h>

#include "utils.hpp"

#define VECT_ALIGNEMENT 32
#define ALIGNED  __attribute__((aligned(VECT_ALIGNEMENT)))

//double vector registries
typedef __m256d dregv;
const size_t dvect_nel=sizeof(dregv)/sizeof(double);
typedef double dvect[dvect_nel] ALIGNED;

//! check that alignement is proper
template <class T> void check_alignement(const T addr)
{
  long int res=(long int)addr%VECT_ALIGNEMENT;
  if(res) CRASH("pointer is %ld offset",res);
}

//! load a vector into a register
inline dregv dvect_load(dvect const addr)
{
  check_alignement(addr);
  return _mm256_load_pd(addr);
}

//! store a register into a vector
inline void dregv_store(dvect addr,dregv v)
{
  check_alignement(addr);
  _mm256_store_pd(addr,v);
}

//! init to zero
inline void dregv_zero(dregv &d)
{d=_mm256_setzero_pd();}

//! summ two register
inline dregv dregv_summ(dregv a,dregv b)
{return _mm256_add_pd(a,b);}

//! subt two register
inline dregv dregv_subt(dregv a,dregv b)
{return _mm256_sub_pd(a,b);}

//! mult two register
inline dregv dregv_prod(dregv a,dregv b)
{return _mm256_mul_pd(a,b);}

//! mult two register and summ
inline dregv dregv_summ_the_prod(dregv a,dregv b,dregv c)
{return _mm256_fmadd_pd(b,c,a);}

//! mult two register and subt
inline dregv dregv_subt_the_prod(dregv a,dregv b,dregv c)
{return _mm256_fnmadd_pd(b,c,a);}

//! swap within 128-lines
inline dregv dregv_shift_to_1032(dregv in)
{return _mm256_permute_pd(in,5);}

//! swap the 128-lines
inline dregv dregv_shift_to_2301(dregv in)
{return _mm256_permute2f128_pd(in,in,7);}

//! makes a full shift dw
inline dregv dregv_shift_to_1230(dregv in)
{
  dregv temp;
  temp=_mm256_permute_pd(in,5);
  return _mm256_blend_pd(temp,_mm256_permute2f128_pd(temp,temp,1),10);
}

//! makes a full shift up
inline dregv dregv_shift_to_3012(dregv in)
{
  dregv temp;
  temp=_mm256_permute_pd(in,5);
  return _mm256_blend_pd(temp,_mm256_permute2f128_pd(temp,temp,1),5);
}

//! shift 0123 into 1032,2301,1230
inline dregv dregv_shift_dw(dregv in,int sh)
{
  switch(sh)
    {
    case 1: return dregv_shift_to_1032(in);
    case 2: return dregv_shift_to_2301(in);
    case 3: return dregv_shift_to_1230(in);
    default: return in;
    }
}

//! shift 0123 into 1032,2301,3012
inline dregv dregv_shift_up(dregv in,int sh)
{
  dregv temp;
  switch(sh)
    {
    case 1: return dregv_shift_to_1032(in);
    case 2: return dregv_shift_to_2301(in);
    case 3: return dregv_shift_to_3012(in);
    default: return in;
    }
}

inline void dvect_fprintf(FILE *fout,dvect c)
{fprintf(fout,"%lg %lg %lg %lg\n",c[0],c[1],c[2],c[3]);}

inline void dvect_printf(dvect c)
{dvect_fprintf(stdout,c);}

///////////////////////////////////////////////////////////////////////////

typedef dvect cdvect[2];
typedef dregv cdregv[2];

#define FOR_RI(a) for(int a=0;a<2;a++)
#define RE 0
#define IM 1

//! load a complex vector into a register
inline void cdvect_load(cdregv out,cdvect const addr)
{FOR_RI(ri) out[ri]=dvect_load(addr[ri]);}

//! store a complex register into a vector
inline void cdregv_store(cdvect addr,cdregv v)
{FOR_RI(ri) dregv_store(addr[ri],v[ri]);}

//! summ two register complex vectors
inline void cdregv_summ(cdregv out,cdregv in1,cdregv in2)
{FOR_RI(ri) out[ri]=dregv_summ(in1[ri],in2[ri]);}

//! subt two register complex vectors
inline void cdregv_subt(cdregv out,cdregv in1,cdregv in2)
{FOR_RI(ri) out[ri]=dregv_subt(in1[ri],in2[ri]);}

//! summ two register complex vectors, multiplying with i the second
inline void cdregv_isumm(cdregv out,cdregv in1,cdregv in2)
{
  out[RE]=dregv_subt(in1[RE],in2[IM]);
  out[IM]=dregv_summ(in1[IM],in2[RE]);
}

//! summ two register complex vectors, multiplying with i the second
inline void cdregv_isubt(cdregv out,cdregv in1,cdregv in2)
{
  out[RE]=dregv_summ(in1[RE],in2[IM]);
  out[IM]=dregv_subt(in1[IM],in2[RE]);
}

//! prod two register complex vectors
inline void cdregv_prod(cdregv out,cdregv in1,cdregv in2)
{
  out[RE]=dregv_subt_the_prod(dregv_prod(in1[RE],in2[RE]),in1[IM],in2[IM]);
  out[IM]=dregv_summ_the_prod(dregv_prod(in1[RE],in2[IM]),in1[IM],in2[RE]);
}

//! prod two register complex vectors, conjugating the first
inline void cdregv_conj1_prod(cdregv out,cdregv in1,cdregv in2)
{
  out[RE]=dregv_summ_the_prod(dregv_prod(in1[RE],in2[RE]),in1[IM],in2[IM]);
  out[IM]=dregv_subt_the_prod(dregv_prod(in1[RE],in2[IM]),in1[IM],in2[RE]);
}

//! prod two register complex vectors, conjugating the second
inline void cdregv_conj2_prod(cdregv out,cdregv in1,cdregv in2)
{cdregv_conj1_prod(out,in2,in1);}

//! summ the prod of two register complex vectors
inline void cdregv_summ_the_prod(cdregv out,cdregv in1,cdregv in2)
{
  out[RE]=dregv_subt_the_prod(dregv_summ_the_prod(out[RE],in1[RE],in2[RE]),in1[IM],in2[IM]);
  out[IM]=dregv_summ_the_prod(dregv_summ_the_prod(out[IM],in1[RE],in2[IM]),in1[IM],in2[RE]);
}

//! summ the prod of two register complex vectors
inline void cdregv_subt_the_prod(cdregv out,cdregv in1,cdregv in2)
{
  out[RE]=dregv_summ_the_prod(dregv_subt_the_prod(out[RE],in1[RE],in2[RE]),in1[IM],in2[IM]);
  out[IM]=dregv_subt_the_prod(dregv_subt_the_prod(out[IM],in1[RE],in2[IM]),in1[IM],in2[RE]);
}

//! summ the prod of two register complex vectors, conjugating the first
inline void cdregv_summ_the_conj1_prod(cdregv out,cdregv in1,cdregv in2)
{
  out[RE]=dregv_subt_the_prod(dregv_subt_the_prod(out[RE],in1[RE],in2[RE]),in1[IM],in2[IM]);
  out[IM]=dregv_summ_the_prod(dregv_subt_the_prod(out[IM],in1[RE],in2[IM]),in1[IM],in2[RE]);
}

//! summ the prod of two register complex vectors, conjugating the second
inline void cdregv_summ_the_conj2_prod(cdregv out,cdregv in1,cdregv in2)
{cdregv_summ_the_conj1_prod(out,in2,in1);}

//! subt the prod of two register complex vectors, conjugating the first
inline void cdregv_subt_the_conj1_prod(cdregv out,cdregv in1,cdregv in2)
{
  out[RE]=dregv_summ_the_prod(dregv_summ_the_prod(out[RE],in1[RE],in2[RE]),in1[IM],in2[IM]);
  out[IM]=dregv_subt_the_prod(dregv_summ_the_prod(out[IM],in1[RE],in2[IM]),in1[IM],in2[RE]);
}

//! subt the prod of two register complex vectors, conjugating the second
inline void cdregv_subt_the_conj2_prod(cdregv out,cdregv in1,cdregv in2)
{cdregv_subt_the_conj1_prod(out,in2,in1);}

//! print a complex vector
inline void cdvect_fprintf(FILE *fout,cdvect in)
{
  for(int i=0;i<4;i++)
    FOR_RI(ri)
      fprintf(fout,"%lg ",in[ri][i]);
  fprintf(fout,"\n");
}
inline void cdvect_printf(cdvect in)
{cdvect_fprintf(stdout,in);}

#endif
