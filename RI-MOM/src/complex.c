#pragma once

//set to zero
void complex_set_to_zero(complex a)
{a[0]=a[1]=0;}

//copy
void complex_copy(complex a,complex b)
{a[0]=b[0];a[1]=b[1];}

//abs
double complex_fabs(complex c)
{return sqrt(c[0]*c[0]+c[1]*c[1]);}

//The summ of two complexs
void complex_summ(complex a,complex b,complex c)
{
  a[0]=b[0]+c[0];
  a[1]=b[1]+c[1];
}
void complex_summassign(complex a,complex b)
{complex_summ(a,a,b);}

//The product of two complex number
void complex_prod(complex a,complex b,complex c)
{
  double tmp=b[0]*c[0]-b[1]*c[1];
  a[1]=b[0]*c[1]+b[1]*c[0];
  a[0]=tmp;
}
void complex_prodassign(complex a,complex b)
{complex_prod(a,a,b);}

void complex_subtassign_the_prod(complex a,complex b,complex c)
{
  double tmp=b[0]*c[0]-b[1]*c[1];
  a[1]-=b[0]*c[1]+b[1]*c[0];
  a[0]-=tmp;
}
void complex_summassign_the_prod(complex a,complex b,complex c)
{
  double tmp=b[0]*c[0]-b[1]*c[1];
  a[1]+=b[0]*c[1]+b[1]*c[0];
  a[0]+=tmp;
}

//reciprocal of a complex
void complex_reciprocal(complex rec,complex c)
{
  double module=c[0]*c[0]+c[1]*c[1];
  
  rec[0]=c[0]/module;
  rec[1]=-c[1]/module;
}
