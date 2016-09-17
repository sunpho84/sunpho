#include "../../src/include.h"

const int njacks=16;
const int nmass=19;
const int ncombo=nmass*(nmass+1)/2;
const double a=0.086/0.197;

int icombo(int im1,int im2)
{
  int imc=max(im1,im2);
  int ims=min(im1,im2);
  
  return ims*nmass-(ims*(ims-1))/2+(imc-ims);
}

int main()
{
  jvec M=jvec_load("results",190,16,0);
  jvec Z2=jvec_load("results",190,16,1);

  int ic=icombo(2,5);
  cout<<"M: "<<M[ic]<<endl;
  cout<<"Z: "<<sqrt(Z2[ic])<<endl;
  cout<<"Mcorr: "<<sqrt(sinh(M[ic])/M[ic])<<endl;
  cout<<ic<<endl;
  cout<<M[ic]<<endl;
  cout<<sqrt(Z2[ic])/(1*0+1*M[ic])/a*0.746<<endl;
  
  return 0;
}
