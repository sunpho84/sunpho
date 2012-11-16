#include "common.cpp"

int T=48,TH=24;
int tmin=12,tmax=23;
double *corr_fit,*corr_err;

//function to fit
double fun_fit(double Z2,double M,int t)
{
  return Z2*exp(-M*TH)*cosh(M*(TH-t))/sinh(M);
}

//chi2 calculation
void chi2(int &npar,double *fuf,double &ch,double *p,int flag)
{
  ch=0;
  for(int t=tmin;t<=tmax;t++)
    ch+=sqr((corr_fit[t]-fun_fit(p[0],p[1],t))/corr_err[t]);
}

int main()
{
  jvec corr(48,16);
  corr.load("/Users/francesco/QCD/LAVORI/NF2/RUN_HIGH/DATA/3.90/24/0.0100/P5P5",0);
  corr=corr.simmetrized(1);
  
  //define minuit staff
  TMinuit minu(2);
  minu.SetFCN(chi2);
  corr_fit=new double[TH+1];
  corr_err=new double[TH+1];
  
  jack E(njack);
  jack Z2(njack);
  two_pts_fit(E,Z2,corr,tmin,tmax);
  
  minu.DefineParameter(0,"Z2",Z2[0],Z2.err(),0,2*Z2[0]);
  minu.DefineParameter(1,"E",E[0],E.err(),0,2*E[0]);
  for(int t=tmin;t<=tmax;t++) corr_err[t]=corr.data[t].err();
  
  //jacknife analysis
  for(int ijack=0;ijack<njack+1;ijack++)
    {
      //copy data so that glob function may access it
      for(int t=tmin;t<=tmax;t++) corr_fit[t]=corr.data[t].data[ijack];
          
      //fit
      double dum;
      minu.Migrad();          
      minu.GetParameter(0,Z2.data[ijack],dum);
      minu.GetParameter(1,E.data[ijack],dum);
    }
  
  cout.precision(7);
  cout<<E<<endl;
  
  jvec data(13,16);
  
  data.load("/Users/francesco/QCD/LAVORI/NF2/RUN_HIGH/ANALYSIS/P5P5/Pi/data/M_Pi",0);
  cout<<data;
  return 0;
}
