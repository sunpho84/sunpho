#include "include.h"

int T=48,njacks=16;

jvec load(const char *path,int im1,int im2)
{
  jvec out(T,njacks);
  
  int ri=0;
  out.load(path,ri+2*(im1+2*im2));
  
  return out.simmetrized(1);
}

int main()
{
  const char path[]="corrs/2pts_P5P5_00_00";
  
  jvec corr_pi=load(path,0,0);
  jvec eff_pi=effective_mass(corr_pi,T/2);

  jack M_pi=constant_fit(eff_pi,13,23,"Pion_fit.xmg");
  
  cout<<M_pi<<endl;
 
  //cout<<corr_pi;
  
  ////
  
  jvec corr_D=load(path,1,0);
  jvec eff_D=effective_mass(corr_D,T/2);

  jack M_D=constant_fit(eff_D,15,23,"D_fit.xmg");
  
  cout<<M_D<<endl;
  
  //cout<<corr_D;
  
  ////
  
  jvec corr_Eta=load(path,1,1);
  jvec eff_Eta=effective_mass(corr_Eta,T/2);

  jack M_Eta=constant_fit(eff_Eta,15,23,"Eta_fit.xmg");
  
  cout<<M_Eta<<endl;
  
  //cout<<corr_Eta;
  
  ////
  
  return 0;
}
