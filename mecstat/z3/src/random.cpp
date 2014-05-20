#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <algorithm>
#include <cmath>

using namespace std;

#include "parameters.hpp"

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

//return a z3 element
z3_t get_random_z3()
{
  return (char)trunc(get_unif_double(3));
}
