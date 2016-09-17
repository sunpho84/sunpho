#include "include.h"

int T,TH;
int njacks=16;

//load
jvec load(const char *name,int r,int ri,int par)
{
  jvec a(T,njacks);
  
  a.load(combine("corrs/2pts_%s_00_00",name).c_str(),ri+2*r);
  
  return a.simmetrized(par);
}

jvec load(const char *name,int par=-1)
{return (load(name,0,0,par)+load(name,1,0,par))/2;}

//multiply by 2*t
jvec fun(jvec corr)
{
  jvec f(corr);
  
  for(int t=0;t<=TH;t++) f[t]*=2*t;
  
  return f;
}

//integrate the correlation function
jvec integrate(jvec corr)
{
  jvec summ(TH+1,njacks);
  
  summ[0]=0;
  for(int t=1;t<=TH;t++)
    summ[t]=summ[t-1]+(corr[t]+corr[t-1])/2;
  
  return summ;
}

jvec spline_integrate(jvec corr)
{
  TMatrixD M(TH+1,TH+1);
  for(int t=0;t<TH+1;t++)
    for(int d=0;d<TH+1;d++)
      if(t==0 && d==0) M(t,d)=0;
      else M(t,d)=pow(t,d);
  
  TMatrixD Minv=M.Invert();
  
  jvec A(TH+1,njacks);
  A=0;
  for(int d=0;d<TH+1;d++)
    {
      jack o(njacks);
      o=0;
      for(int t=0;t<TH+1;t++) o+=Minv(d,t)*corr[t];
      for(int t=0;t<TH+1;t++) A[t]+=o*pow(t,d+1)/(d+1);
    }
  
  return A;
}

int main(int narg,char **arg)
{
  if(narg<2) crash("use %s T",arg[0]);
  T=atoi(arg[1]);
  TH=T/2;
  
  jvec corr_VKTK=+(load("V1T1")+load("V2T2")+load("V3T3"))/3;
  jvec corr_TKVK=-(load("T1V1")+load("T2V2")+load("T3V3"))/3;

  jvec sum_TKVK=integrate(fun(corr_TKVK));
  jvec sum_VKTK=integrate(fun(corr_VKTK));
  
  sum_TKVK.print_to_file("/tmp/TKVK_inte.xmg");
  sum_VKTK.print_to_file("/tmp/VKTK_inte.xmg");
  
  (sum_VKTK-sum_TKVK).print_to_file("/tmp/diff_inte.xmg");
  
  /*
  jvec sp_sum_TKVK=spline_integrate(fun(corr_TKVK));
  jvec sp_sum_VKTK=spline_integrate(fun(corr_VKTK));
  
  sp_sum_TKVK.print_to_file("/tmp/TKVK_sp_inte.xmg");
  sp_sum_VKTK.print_to_file("/tmp/VKTK_sp_inte.xmg");
  
  (sp_sum_VKTK-sp_sum_TKVK).print_to_file("/tmp/diff_sp_inte.xmg");
  */
  
  return 0;
}
