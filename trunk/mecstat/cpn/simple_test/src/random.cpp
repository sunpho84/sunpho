#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <algorithm>
#include <iostream>
#include <cmath>

using namespace std;

#include "data.hpp"
#include "debug.hpp"
#include "geometry.hpp"
#include "init.hpp"
#include "parameters.hpp"
#include "routines.hpp"
#include "tools.hpp"

#define EXTERN_RANDOM
#include "random.hpp"

//compute p(theta)
double fun_ptheta(double theta,double a,int k)
{
  //compute the factor
  double f=sin(theta);
  f*=f;

  //compute first piece
  double p=1;
  for(int i=0;i<k-1;i++) p*=f;
  
  return p*exp(a*cos(theta));
}

#ifdef SMART_EXTRACTION
//obtain theta
double get_theta(double a,int k,int site)
{
  //compute parameters
  double zita=(k-1)/a;
  double theta0=acos(sqrt(1+zita*zita)-zita);
  double ctheta0=cos(theta0);
  double c=sqrt(2*(k-1)*(1-zita*ctheta0)/sqr(sin(theta0)));
  double ptheta0=fun_ptheta(theta0,a,k);
  
  //extract theta
  double theta;
  double eta=0.99;
  bool acc;
  int no=0;
  do
    {
      double chi=get_unif_double(1,site);
      theta=theta0+tan(chi*atan(c*(M_PI-theta0))+(chi-1)*atan(c*theta0))/c;
      
      //reweighting
      double ptheta=fun_ptheta(theta,a,k);
      double pacc=ptheta/ptheta0*(1+sqr(c*(theta-theta0)))*eta;
      double extr=get_unif_double(1,site);
      if(pacc>1) CRASH("pacc: %lg",pacc);
      //cerr<<extr<<" "<<theta<<" "<<theta0<<" "<<zita<<endl;
      acc=(extr<pacc);
      if(!acc) no++;
      if(no>1000) CRASH("k: %lg, a: %d, zita: %lg, ptheta: %lg, ptheta0: %lg, c: %lg, theta: %lg, theta0: %lg, pacc: %lg",
			k,a,zita,ptheta,ptheta0,c,theta,theta0,pacc);
    }
  while(!acc);
  
  return theta;
}

//taken by appendix C of hep-lat/9210016
//correction done: h is returned in place of (wrong) g
double eget_theta_1(double a,int site)
{
  double eps=0.001;
  double as=0.798953686083986;
  double dap=max((double)0.0,a-as);
  double del=0.35*dap+1.03*sqrt(dap);
  double alp=min(sqrt(a*(2-eps)),max(sqrt(eps*a),del));
  double bet=max(alp*alp/a,(double)((cosh(M_PI*alp)-1)/(exp(2*a)-1)))-1;
  double bt1=sqrt((1+bet)/(1-bet));
  
  double h; //result
  bool acc; //accepted or not
  int no=0;
  do
    {
      double r=get_unif_double(1,site);
      double h1=bt1*tan((2*r-1)*atan(tanh(M_PI*alp/2)/bt1));
      h=log((1+h1)/(1-h1))/alp;
      
      //decide if accept or reject
      double g=exp(-a*(1-cos(h)))*(cosh(alp*h)+bet)/(1+bet);
      if(g>1+TINY) CRASH("%lg",g-1);
      double p=get_unif_double(1,site);
      acc=(p<g);
      if(!acc) no++;
      if(no>1000) CRASH("%d",no);
    }
  while(!acc);
  
  cerr<<h<<endl;
  
  return h;
}

double get_theta_1(double a,int site)
{
  int eff=0;
  
  double max=exp(fabs(a));
  double ret,pacc,extr;
  do
    {
      ret=get_unif_double(2*M_PI,site)-M_PI;
      pacc=exp(a*cos(ret));
      extr=get_unif_double(max,site);
      
      if(pacc>max) crash("ahm");
      eff++;
    }
  while(extr>pacc);

  //cerr<<eff<<endl;
  return ret;
}
#else
double get_theta(double a,int k,int site)
{
  double th,no=exp(-fabs(a));
  
  int non=0;
  bool acc;
  do
    {
      th=get_unif_double(2*M_PI,site);
      double pacc=no*pow(sin(th),2*(k-1))*exp(a*cos(th));
      if(pacc>1) CRASH("a%lg",pacc);
      double ext=get_unif_double(1,site);
      acc=(ext<pacc);
      if(!acc) non++;
      if(non>10000) CRASH("non: %d, a: %lg",non,a);
    }
  while(!acc);
  
  return th;
}

double get_theta_1(const double a,int site)
{
  double th,no=exp(-fabs(a));
  
  int non=0;
  bool acc;
  do
    {
      th=get_unif_double(2*M_PI,site);
      double pacc=no*exp(a*cos(th));
      if(pacc>1+TINY) CRASH("a%lg",pacc-1);
      double ext=get_unif_double(1,site);
      acc=(ext<pacc);
      if(!acc) non++;
      if(non>100) CRASH("non: %d",non);
    }
  while(!acc);
  
  return th;
}
#endif

//set an U1 to random
void set_U1_to_rnd(dcomplex &U,int site)
{
  //extract a phase
  double ph=get_unif_double(2*M_PI,site);
  U=dcomplex(cos(ph),sin(ph));
}

//set an O(N) to random respecting the bound of unitarity
void set_ON_to_rnd(dcomplex *O,int site)
{
  //first of all extract their norm
  double w[N+1];
  w[0]=0;
  for(int i=1;i<N;i++) w[i]=get_unif_double(1,site);
  w[N]=1;
  sort(w,w+N+1);
  
  //extract a random complex number for each site using the extracted norm
  for(int i=0;i<N;i++)
    {
      double nor=sqrt(w[i+1]-w[i]);
      double the=get_unif_double(2*M_PI,site);
      O[i]=nor*dcomplex(cos(the),sin(the));
    }
}
