#include "../../src/include.h"

void read_pseudo(int T,int njacks)
{
  FILE *fin=open_file("/tmp/p","r");
  jvec corr(T,njacks);
  corr=0;
  
  for(int iconf=0;iconf<njacks;iconf++)
    for(int t=0;t<T;t++)
      if(fscanf(fin,"%lg",&corr[t][iconf])!=1) crash("reading ijack %d t %d",iconf,t);
  fclose(fin);
  
  corr.clusterize();
  
  corr.simmetrized(1).print_to_file("/tmp/pseudo_corr.xmg");
  effective_mass(corr.simmetrized(1)).print_to_file("/tmp/pseudo_corr_eff.xmg");
}

int T,TH;
const int njacks=16;

const int dt_fit=4;
int t_osc_fit;
double *c_osc_two_pts_fit,*e_osc_two_pts_fit;

template <class Ti> Ti fun_osc_two_pts_migrad_fit(Ti Z2A,Ti Z2B,Ti MA,Ti MB,int t)
{
  int s=pow(-1.0,t+1); // this +1 is in line with Massimo thesis where A- is negative
  return Z2A*(exp(-MA*t)-s*exp(-MA*(T-t)))+s*Z2B*(exp(-MB*t)-s*exp(-MB*(T-t)));
}

void ch2_osc_two_pts_migrad_fit(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  int t0=t_osc_fit;
  double Z2A=p[0];
  double Z2B=p[1];
  double MA=p[2];
  double MBA=p[3];
  double MB=MA+MBA;
  
  for(int t=t0;t<t0+dt_fit*2;t+=2)
    {
      double num=c_osc_two_pts_fit[t];
      double teo=fun_osc_two_pts_migrad_fit(Z2A,Z2B,MA,MB,t);
      double diff=num-teo;
      double err=e_osc_two_pts_fit[t];
      double cont=sqr(diff/err);
      ch+=cont;
      if(flag==3)
	cout<<" t="<<t<<", diff=("<<num<<"-"<<teo<<")="<<diff<<" err="<<err<<" cont="<<cont<<endl;
    }
}

void oscill_eff_mass(jvec &Z2A,jvec &Z2B,jvec &MA,jvec &MB,jvec corr)
{
  Z2A=jvec(TH-dt_fit+1,njacks);
  Z2B=jvec(TH-dt_fit+1,njacks);
  MA=jvec(TH-dt_fit+1,njacks);
  MB=jvec(TH-dt_fit+1,njacks);
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(ch2_osc_two_pts_migrad_fit);

  //cout<<corr[0][0]<<endl;
  minu.DefineParameter(0,"Z2A",0.01,0.005,0,1);
  minu.DefineParameter(1,"Z2B",0.01,0.005,0,1);
  
  c_osc_two_pts_fit=new double[TH+1];
  e_osc_two_pts_fit=new double[TH+1];
  for(int iel=0;iel<=TH;iel++) e_osc_two_pts_fit[iel]=corr[iel].err();
  
  for(int ijack_fit=0;ijack_fit<=njacks;ijack_fit++)
    {
      for(int iel=0;iel<=TH;iel++) c_osc_two_pts_fit[iel]=corr[iel][ijack_fit];
      
      for(t_osc_fit=0;t_osc_fit<=TH-dt_fit;t_osc_fit++)
	{
	  double a=0.8,ba=0.5;
	  minu.DefineParameter(2,"MA",a,0.1,0,5);
	  minu.DefineParameter(3,"MBA",ba,0.2,0,5);
	  
	  minu.Migrad();
	  double dum;
	  minu.GetParameter(0,Z2A[t_osc_fit][ijack_fit],dum);
	  minu.GetParameter(1,Z2B[t_osc_fit][ijack_fit],dum);
	  minu.GetParameter(2,a,dum);
	  minu.GetParameter(3,ba,dum);
	  double b=a+ba;
	  
	  MA[t_osc_fit][ijack_fit]=a;
	  MB[t_osc_fit][ijack_fit]=b;
	  double ch2,grad[4],par[4];
	  for(int i=0;i<4;i++) minu.GetParameter(i,par[i],dum);
	  minu.Eval(4,grad,ch2,par,3);
	  cout<<"Z2A: "<<par[0]<<", Z2B: "<<par[1]<<", MA: "<<a<<", MB: "<<b<<", ch2: "<<ch2<<endl;
  	}
    }
}

int main(int narg,char **arg)
{
  //read_pseudo(48,700);
  
  if(narg<4) crash("use %s file nconf T",arg[0]);
  
  const char *path=arg[1];
  int nconfs=atoi(arg[2]);
  cout<<"nconfs possible: "<<nconfs<<endl;
  int clust_size=nconfs/njacks;
  nconfs=clust_size*njacks;
  cout<<"nconfs: "<<nconfs<<endl;
  T=atoi(arg[3]);
  TH=T/2;
  
  FILE *fin=open_file(path,"r");
  
  jvec corr(T,njacks);
  corr=0;
  
  for(int iconf=0;iconf<nconfs;iconf++)
    for(int t=0;t<T;t++)
      {
	double temp;
	if(fscanf(fin,"%lg",&temp)!=1) crash("reading ijack %d t %d",iconf,t);
	corr[t][iconf/clust_size]+=temp;
      }
  
  corr.clusterize(clust_size);
  //corr=corr.simmetrized(+1);
  
  jvec Z2A,Z2B,MA,MB;
  oscill_eff_mass(Z2A,Z2B,MA,MB,corr);
  
  Z2A.print_to_file("/tmp/Z2A.xmg");
  Z2B.print_to_file("/tmp/Z2B.xmg");
  MA.print_to_file("/tmp/MA.xmg");
  MB.print_to_file("/tmp/MB.xmg");
  
  ofstream outA("/tmp/MA_jack");
  ofstream outB("/tmp/MB_jack");
  for(int t=0;t<=TH-dt_fit;t++)
    {
      for(int ijack=0;ijack<=njacks;ijack++) outA<<MA[t][ijack]<<endl;
      for(int ijack=0;ijack<=njacks;ijack++) outB<<MB[t][ijack]<<endl;
      outA<<"&"<<endl;
      outB<<"&"<<endl;
    }
  
  
  corr.print_to_file("/tmp/nucleon_corr.xmg");
  effective_mass(corr,-1,1,2).print_to_file("/tmp/nucleon_corr_eff_mass.xmg");
  
  return 0;
}

