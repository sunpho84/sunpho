#pragma once

double unary_minus(const double a){return -a;}
double unary_plus(const double a){return a;}
double sqr(const double a){return a*a;}

double double_summ(const double a,const double b){return a+b;}
double double_subt(const double a,const double b){return a-b;}
double double_prod(const double a,const double b){return a*b;}
double double_frac(const double a,const double b){return a/b;}

double min(double a,double b)
{
  if(a>=b) return b;
  else return a;
}

double max(double a,double b)
{
  if(a>=b) return a;
  else return b;
}

/////////////////////////////////////////////////////////////

double lin_fun(double const x,double*p)
{return p[0]*x+p[1];}

double const_fun(double const x,double*p)
{return p[0];}

template <class T> T cont_e(T m,double p)
{return sqrt(m*m+3*p*p);}

template <class T> T latt_e(T m,double p)
{return 2*asinh(sqrt(3*sqr(sin(p/2))+sqr(sinh(m/2))));}
