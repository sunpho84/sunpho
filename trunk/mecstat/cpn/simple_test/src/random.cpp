#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <algorithm>
#include <math.h>

using namespace std;

#include "data.hpp"
#include "debug.hpp"
#include "geometry.hpp"
#include "parameters.hpp"
#include "routines.hpp"

#define EXTERN_RANDOM
#include "random.hpp"

//return a double between [0,max)
double get_unif_double(double max,bool incl)
{
  double res;
#ifdef GOOD_GENERATOR
  do res=max*(*dis)(*gen)/dis->b();
  while((!incl)&&(res==max));
#else
  const int im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014,ia2=40692;
  const int iq1=53668,iq2=52774,ir1=12211,ir2=3791,ndiv=1+imm1/RAN2_NTAB;
  const double am=1.0/im1,eps=1.2e-7,rnmx=1-eps;
  int j,k;
    
  k=gen.idum/iq1;
  gen.idum=ia1*(gen.idum-k*iq1)-k*ir1;
  if(gen.idum<0) gen.idum+=im1;
    
  k=gen.idum2/iq2;
  gen.idum2=ia2*(gen.idum2-k*iq2)-k*ir2;
  if(gen.idum2<0) gen.idum2+=im2;
    
  j=gen.iy/ndiv;
  gen.iy=gen.iv[j]-gen.idum2;
  gen.iv[j]=gen.idum;
  if(gen.iy<0) gen.iy+=imm1;
    
  res=max*std::min(am*gen.iy,rnmx);
#endif
  
  return res;
}

//extract gaussianly with *total* standard deviation 1
dcomplex get_gauss_complex()
{
#ifndef M_SQRT_2
 #define M_SQRT_2 0.707106781186547524401
#endif

  double r=M_SQRT_2*sqrt(-2*log(1-get_unif_double(1)));
  double q=2*M_PI*get_unif_double(1);
    
  return dcomplex(r*cos(q),r*sin(q));
}

//return a gaussian variable with standard deviation 1
double get_gauss_double()
{
  static bool flag=false;
  static dcomplex ref;
  
  //at odd turns extract
  if(flag==false) ref=get_gauss_complex();
  flag=!flag;
  
  if(flag) return ref.real()*M_SQRT2;
  else     return ref.imag()*M_SQRT2;
}

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
double get_theta(double a,int k)
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
      double chi=get_unif_double(1);
      theta=theta0+tan(chi*atan(c*(M_PI-theta0))+(chi-1)*atan(c*theta0))/c;
      
      //reweighting
      double ptheta=fun_ptheta(theta,a,k);
      double pacc=ptheta/ptheta0*(1+sqr(c*(theta-theta0)))*eta;
      double extr=get_unif_double(1);
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
double get_theta_1(double a)
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
      double r=get_unif_double(1);
      double h1=bt1*tan((2*r-1)*atan(tanh(M_PI*alp/2)/bt1));
      h=log((1+h1)/(1-h1))/alp;
      
      //decide if accept or reject
      double g=exp(-a*(1-cos(h)))*(cosh(alp*h)+bet)/(1+bet);
      if(g>1+TINY) CRASH("%lg",g-1);
      double p=get_unif_double(1);
      acc=(p<g);
      if(!acc) no++;
      if(no>1000) CRASH("%d",no);
    }
  while(!acc);
  
  return h;
}
#else
double get_theta(double a,int k)
{
  double th,no=exp(-fabs(a));
  
  int non=0;
  bool acc;
  do
    {
      th=get_unif_double(2*M_PI);
      double pacc=no*pow(sin(th),2*(k-1))*exp(a*cos(th));
      if(pacc>1) CRASH("a%lg",pacc);
      double ext=get_unif_double(1);
      acc=(ext<pacc);
      if(!acc) non++;
      if(non>10000) CRASH("non: %d, a: %lg",non,a);
    }
  while(!acc);
  
  return th;
}

double get_theta_1(const double a)
{
  double th,no=exp(-fabs(a));
  
  int non=0;
  bool acc;
  do
    {
      th=get_unif_double(2*M_PI);
      double pacc=no*exp(a*cos(th));
      if(pacc>1+TINY) CRASH("a%lg",pacc-1);
      double ext=get_unif_double(1);
      acc=(ext<pacc);
      if(!acc) non++;
      if(non>100) CRASH("non: %d",non);
    }
  while(!acc);
  
  return th;
}
#endif

//set an U1 to random
void set_U1_to_rnd(dcomplex &U)
{
  //extract a phase
  double ph=get_unif_double(2*M_PI);
  U=dcomplex(cos(ph),sin(ph));
}

//set an O(N) to random respecting the bound of unitarity
void set_ON_to_rnd(dcomplex *O)
{
  //first of all extract their norm in such: ordering
  double w[N+1];
  w[0]=0;
  for(int i=1;i<N;i++) w[i]=get_unif_double(1);
  w[N]=1;
  sort(w,w+N+1);
  
  //extract a random complex number for each site using the extracted norm
  for(int i=0;i<N;i++)
    {
      double nor=sqrt(w[i+1]-w[i]);
      double the=get_unif_double(2*M_PI);
      O[i]=nor*dcomplex(cos(the),sin(the));
    }
}

//initialize the system to hot
void init_system_to_hot()
{
  for(int site=0;site<V;site++)
    {
      //fill the lambda
      for(int mu=0;mu<NDIMS;mu++)
	set_U1_to_rnd(lambda[site*NDIMS+mu]);
      
      //fill the Zeta
      set_ON_to_rnd(zeta+site*N);
    }
}

//initialize to cold
void init_system_to_cold()
{
  for(int site=0;site<V;site++)
    {
      //fill the lambda
      for(int mu=0;mu<NDIMS;mu++) lambda[site*NDIMS+mu]=1;
      
      //fill the Zeta
      for(int n=0;n<N;n++) zeta[site*N+n]=(n==0);
    }
}

//switch
void init_system_to(int cond)
{
  if(cond==HOT) init_system_to_hot();
  else          init_system_to_cold();
}
