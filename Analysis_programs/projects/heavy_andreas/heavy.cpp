#include "include.h"
#include "../jpsiprime/eff_bis.cpp"

const int T=96,TH=48;
int nconfs=36,njack=nconfs;

jvec load(const char *path)
{
  ifstream fin(path);
  if(!fin.good()) crash("opening");
  
  jvec f(T,nconfs);
    
  int i;
  double dum;
  for(int iconf=0;iconf<nconfs;iconf++)
    for(int t=0;t<T;t++)
      fin>>i>>i>>i>>i>>i>>i>>f[t][iconf]>>dum;
  if(!fin.good()) crash("something went wrong");
  
  f.clusterize();
  
  return f.simmetrized(1);
}

int main(int narg,char **arg)
{
  //jvec corr=load("/Users/francesco/Downloads/meson_ms_0.0133_LL_0.0_t0_mh_0.28_LL_0.0_t0_pp");
  jvec corrA=load("/Users/francesco/Downloads/meson_mh_0.28_LL_0.0_t0_mh_0.28_LL_0.0_t0_pp");
  
  jack MA,ZA;
  int tA=47;
  two_pts_fit(MA,ZA,corrA,tA,TH,"/tmp/MA.xmg","/tmp/ZA.xmg");
  
  jvec corrB=corrA;
  for(int t=0;t<=TH;t++) corrB[t]-=ZA*exp(-MA*TH)*cosh(MA*(TH-t))/MA;
  
  jack MB,ZB;
  int tB=30;
  two_pts_fit(MB,ZB,corrB,tB,tA-2,"/tmp/MB.xmg","/tmp/ZB.xmg");
  
  two_states_fit(MA,MB,ZA,ZB,tB,TH,corrA,"/tmp/MAB.xmg");
    
  cout<<MA<<endl;
  cout<<MB<<endl;

  return 0;
}
