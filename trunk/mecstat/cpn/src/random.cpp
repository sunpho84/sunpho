#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <algorithm>

#include "debug.hpp"
#include "random.hpp"

//initialize
void rnd_gen_t::init(int seed)
{
  const int im1=2147483563,ia1=40014;
  const int iq1=53668,ir1=12211;
  int j,k;
  
  //initialization
  idum=seed;
  idum=std::max(idum+1,1);
  idum2=idum;
  for(j=RAN2_NTAB+7;j>=0;j--)
    {
      k=idum/iq1;
      idum=ia1*(idum-k*iq1)-k*ir1;
      if(idum<0) idum+=im1;
      if(j<RAN2_NTAB) iv[j]=idum;
    }
  iy=iv[0];
  
  //mark to be inited
  inited=true;
}

//standard ran2 from numerical recipes
double rnd_gen_t::get_unif(double min,double max)
{
  if(!inited) CRASH_SOFTLY("cannot work if not properly inited");
  
  const int im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014,ia2=40692;
  const int iq1=53668,iq2=52774,ir1=12211,ir2=3791,ndiv=1+imm1/RAN2_NTAB;
  const double am=1.0/im1,eps=1.2e-7,rnmx=1-eps;
  int j,k;
  double out;
  
  k=idum/iq1;
  idum=ia1*(idum-k*iq1)-k*ir1;
  if(idum<0) idum+=im1;
  
  k=idum2/iq2;
  idum2=ia2*(idum2-k*iq2)-k*ir2;
  if(idum2<0) idum2+=im2;
  
  j=iy/ndiv;
  iy=iv[j]-idum2;
  iv[j]=idum;
  if(iy<0) iy+=imm1;
  
  out=std::min(am*iy,rnmx);
  
  return out*(max-min)+min;
}

//return a numer between 0 and 1
int rnd_gen_t::get_pm_one()
{
  double r=get_unif(0,1);
  if(r>=0.5) return 1;
  else return -1;
}

