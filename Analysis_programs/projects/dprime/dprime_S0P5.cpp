#include <math.h>
#include "include.h"

const int njacks=16;
int T=48,TH=T/2,L=T/2;
int nml=3;

jvec load_S0P5(const char *path,int im1,int im2)
{
  int ri=1;
  int r=0;
  jvec a(T,njacks),b(T,njacks);
  
  a.load(path,ri+2*( r+2*(im1+nml*( r+2*im2))));
  b.load(path,ri+2*(!r+2*(im1+nml*(!r+2*im2))));
  
  return ((b-a)/2).simmetrized(1);
}

jvec load_P5P5(const char *path,int im1,int im2)
{
  int ri=0;
  int r=0;
  jvec a(T,njacks),b(T,njacks);
  
  a.load(path,ri+2*( r+2*(im1+nml*( r+2*im2))));
  b.load(path,ri+2*(!r+2*(im1+nml*(!r+2*im2))));
  
  return ((a+b)/2).simmetrized(1);
}

int main(int narg,char **arg)
{
  jvec corr_S0P5=load_S0P5("corrs/2pts_S0P5_02_00",2,0);
  jvec corr_S0S0=-load_P5P5("corrs/2pts_S0S0_02_00",2,0);
  jvec corr_P5P5=load_P5P5("corrs/2pts_P5P5_02_00",2,0);
  
  ofstream out("test.xmg");
  out<<"@type xydy"<<endl;
  out<<effective_mass(corr_S0P5)<<"&"<<endl;
  out<<effective_mass(corr_P5P5)<<"&"<<endl;
  out<<effective_mass(corr_S0S0)<<"&"<<endl;
  
  return 0;
}
