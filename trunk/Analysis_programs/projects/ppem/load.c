#include "include.h"

int njack=16;
int L=24,T=48,TH=T/2;
int spat_vol=L*L*L;

int main()
{
  for(int i=0;i<1;i++)
    {
      jvec a(T,njack),b(T,njack);
      
      a.load("oPPo-ss_conf.1.dat",2*i);
      a/=-spat_vol;
      b.load("oPPo-sd_conf.1.dat",2*i);
      b/=spat_vol;
      a=effective_mass(a.simmetrized(1));
      ofstream tmp(combine("/tmp/garb/%03d",i).c_str());
      tmp<<a<<endl;
      cout<<constant_fit(a,15,TH,"/tmp/y")<<endl;
    }
  return 0;
}
