#include "include.h"

const int T=48,TH=T/2;
const int njacks=16;

inline double ph2(int t)
{return sqr(sin(2*M_PI*t/T/2));}

int main()
{
  jvec a=jvec_load("time_momentum_prop",T,njacks,2);
  jvec b=jvec_load("time_momentum_prop",T,njacks,3);

  a.print_to_file("prop_re.xmg");
  b.print_to_file("prop_im.xmg");
  
  jvec pr(T,njacks);
  jvec pi(T,njacks);
  for(int t=0;t<T;t++)
    {
      double c=cos(0.5*t*2*M_PI/T);
      double s=sin(0.5*t*2*M_PI/T);
      
      pr[t]=a[t]*c-b[t]*s;
      pi[t]=a[t]*s+b[t]*c;
    }
  
  pr=0.5*(pr-pr.simmetric());
  pi=0.5*(pi+pi.simmetric());
  
  pr.print_to_file("prop_re.xmg");
  pi.print_to_file("prop_im.xmg");
  
  jvec xr(T,njacks);
  jvec xi(T,njacks);
  for(int t=0;t<T;t++)
    {
      double c=cos(t*2*M_PI/T);
      double s=sin(t*2*M_PI/T);
      
      xr[t]=pr[t]*c-pi[t]*s;
      xi[t]=pr[t]*s+pi[t]*c;
    }
  
  xr.print_to_file("prop_x.xmg");
  xi.print_to_file("tooth.xmg");
  
  //out<<ph2(t)<<" "<<a[t]<<endl;
  
  return 0;
}
