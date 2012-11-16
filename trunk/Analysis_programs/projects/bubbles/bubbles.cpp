#include <include.h>
#include <iostream>
#include <sstream>

using namespace std;

int T=48,TH=T/2;
int njack=16;

jvec compute_prd_tra(jvec evn,jvec odd)
{
  jvec tra=evn*0;
  for(int dt=0;dt<T;dt++)
    {
      tra[dt]=0;
      for(int t1=0;t1<T;t1++)
	{
	  int t2=(t1+dt)%T;
	  tra[dt]+=evn[t1]*odd[t2];
	}
      tra[dt]/=T;
    }
  
  return tra;
}

int tmin=3,tmax=6;
double *corr_med,*corr_err;

//function to fit
double fun_fit(double Z2,double M,double C,int t)
{
  return Z2*exp(-M*TH)*cosh(M*(TH-t))/sinh(M)+C;
}

//chi2 calculation
void chi2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  for(int t=tmin;t<=tmax;t++)
    ch+=sqr((corr_med[t]-fun_fit(p[0],p[1],p[2],t))/corr_err[t]);
}


void compute()
{
  double V=sqrt(TH*TH*TH);
  jvec evn=(jvec_load("bubble_evn_P5",T,njack,1)-jvec_load("bubble_evn_P5",T,njack,3))/2*V;
  jvec odd=(jvec_load("bubble_odd_P5",T,njack,1)-jvec_load("bubble_odd_P5",T,njack,3))/2*V;
  jvec prd=(jvec_load("bubble_prd_P5",T,njack,1)+jvec_load("bubble_prd_P5",T,njack,3)).simmetrized(1)/2*V*V;

  jvec tra=compute_prd_tra(evn,odd).simmetrized(1);
  jvec bub=(prd-tra);
  
  jvec conn=(jvec_load("2pts_00_00_P5P5",T,njack,0)+jvec_load("2pts_00_00_P5P5",T,njack,3)).simmetrized(1)/2;
  
  TMinuit minu;
  minu.SetPrintLevel(-1);
  minu.SetFCN(chi2);
  corr_med=new double[TH+1];
  corr_err=new double[TH+1];
  jvec Mcor=effective_mass(bub),Z2cor(TH+1,njack);
  jack Meff=constant_fit(Mcor,tmin,tmax);
  for(int t=0;t<=TH;t++)
    for(int ijack=0;ijack<=njack;ijack++)
      Z2cor[t].data[ijack]=bub[t].data[ijack]/fun_fit(1,Meff[ijack],0,t);
  jack Z2eff=constant_fit(Z2cor,tmin,tmax);
        
  minu.DefineParameter(0,"Z2",Z2eff[0],Z2eff.err(),0,2*Z2eff[0]);
  minu.DefineParameter(1,"M",Meff[0],Meff.err(),0,2*Meff[0]);
  minu.DefineParameter(2,"C",0,0.01,0,0);
  for(int t=tmin;t<=tmax;t++) corr_err[t]=bub.data[t].err();
  //jacknife analysis
  jack Z2(njack),M(njack),C(njack);
  for(int ijack=0;ijack<njack+1;ijack++)
    {
      //copy data so that glob function may access it
      for(int t=tmin;t<=tmax;t++) corr_med[t]=bub.data[t].data[ijack];
      
      //fit
      double dum;
      minu.Migrad();          
      minu.GetParameter(0,Z2.data[ijack],dum);
      minu.GetParameter(1,M.data[ijack],dum);
      minu.GetParameter(2,C.data[ijack],dum);
    }
  
  cout<<"M: "<<M<<endl;
  cout<<"C: "<<C<<endl;
  bub-=C;
  
  ofstream out("out.xmg");
  out.precision(16);
  out<<"@type xydy"<<endl;
  //out<<effective_mass(bub)<<endl;
  out<<effective_mass(conn)<<"&\n"<<effective_mass(conn-bub);
}

void chris_compute()
{
  jvec evn=jvec_load("/tmp/bub/chris_bubble_evn_CHRIS-P5",T,njack,0);
  jvec odd=jvec_load("/tmp/bub/chris_bubble_odd_CHRIS-P5",T,njack,0);
  jvec prd=jvec_load("/tmp/bub/chris_bubble_prd_CHRIS-P5",T,njack,0);
  
  jvec tra=compute_prd_tra(evn,odd);
  
  cout<<(prd-tra).simmetrized(1)/4<<endl;
}

int main(int narg,char **arg)
{
  compute();
  
  return 0;
}
