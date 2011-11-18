#pragma once

void ss_set_to_zero(spinspin a)
{memset(a,0,sizeof(spinspin));}

void double_vect_summ(double *a,double *b,double *c,int n)
{for(int i=0;i<n;i++) a[i]=b[i]+c[i];}
void double_vect_subt(double *a,double *b,double *c,int n)
{for(int i=0;i<n;i++) a[i]=b[i]-c[i];}
void double_vect_prod_real(double *a,double *b,double c,int n)
{for(int i=0;i<n;i++) a[i]=b[i]*c;}

void ccss_propagator_summ(ccss_propagator a,ccss_propagator b,ccss_propagator c)
{double_vect_summ((double*)a,(double*)b,(double*)c,sizeof(ccss_propagator)/sizeof(double));}
void ccss_propagator_subt(ccss_propagator a,ccss_propagator b,ccss_propagator c)
{double_vect_subt((double*)a,(double*)b,(double*)c,sizeof(ccss_propagator)/sizeof(double));}
void ccss_propagator_prod_real(ccss_propagator a,ccss_propagator b,double c)
{double_vect_prod_real((double*)a,(double*)b,c,sizeof(ccss_propagator)/sizeof(double));}

void ccss_propagator_summassign(ccss_propagator a,ccss_propagator b)
{ccss_propagator_summ(a,a,b);}
void ccss_propagator_subtassign(ccss_propagator a,ccss_propagator b)
{ccss_propagator_subt(a,a,b);}
void ccss_propagator_prodassign_real(ccss_propagator a,double b)
{ccss_propagator_prod_real(a,a,b);}

void cc_trace_ccss(spinspin a,su3spinspin b)
{
  ss_set_to_zero(a);
  for(int id1=0;id1<4;id1++)
    for(int id2=0;id2<4;id2++)
      for(int ic=0;ic<3;ic++)
	complex_summassign(a[id1][id2],b[ic][ic][id1][id2]);
}

void ccss_trace_ccss(complex a,su3spinspin b)
{
  complex_set_to_zero(a);
  for(int id=0;id<4;id++)
    for(int ic=0;ic<3;ic++)
      complex_summassign(a,b[ic][ic][id][id]);
}

void cc_trace_ccss_propagator(ss_propagator a,ccss_propagator b)
{for(int imom=0;imom<nmom;imom++) cc_trace_ccss(a[imom],b[imom]);}

void ccss_trace_ccss_propagator(cmom a,ccss_propagator b)
{for(int imom=0;imom<nmom;imom++) ccss_trace_ccss(a[imom],b[imom]);}
