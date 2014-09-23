#include <cmath>
#include <iostream>

using namespace std;

template <class T> T sqr(T in)
{return in*in;}

double cont_e(double m,double p)
{return sqrt(m*m+3*p*p);}

double latt_e(double m,double p)
{return 2*asinh(sqrt(3*sqr(sin(p/2))+sqr(sinh(m/2))));}

double Q2_fun(double m1,double p1,double m2,double p2,double(*fun_e)(double,double)=cont_e)
{
  double e1=fun_e(m1,p1);
  double e2=fun_e(m2,p2);
  
  return (e2-e1)*(e2-e1)-3*(p2-p1)*(p2-p1);
}

double p2_fun_bf(double p1)
{return -p1;}
double p2_fun_2rf(double p1)
{return 0;}

double find_p1(double m1,double m2,double(fun_p2)(double),double(*fun_e)(double,double)=cont_e)
{
  double eps=0.1;
  double p1=0,p2=0;
  double Q2=Q2_fun(m1,p1,m2,p2,fun_e),Q2_old;
  do
    {
      p1+=eps;
      p2=fun_p2(p1);
      
      Q2_old=Q2;
      Q2=Q2_fun(m1,p1,m2,p2,fun_e);
      
      if(Q2*Q2_old<0) eps*=-0.5;
      else
	if(fabs(Q2_old)<fabs(Q2)) eps*=-1;
    }
  while(fabs(Q2)>1.e-14);

  return p1;
}

jack find_p1(jack m1,jack m2,double(*fun_p2)(double),double(*fun_e)(double,double)=cont_e)
{
  int njacks=m1.njack;
  jack out(njacks);
  for(int ijack=0;ijack<=njacks;ijack++) out[ijack]=find_p1(m1[ijack],m2[ijack],fun_p2,fun_e);
  
  return out;
}

template <class T> T find_p1_bf(T m1,T m2,double(*fun_e)(double,double)=latt_e)
{return find_p1(m1,m2,p2_fun_bf,fun_e);}

template <class T> T find_p1_2rf(T m1,T m2,double(*fun_e)(double,double)=latt_e)
{return find_p1(m1,m2,p2_fun_2rf,fun_e);}
