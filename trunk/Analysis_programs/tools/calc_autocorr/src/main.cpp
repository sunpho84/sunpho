#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>

#include "driver.hpp"

using namespace std;

const int clust_size=1;
int jackniffed_size;
int njacks;
int size;
double *data;

//read allocating
void read(const char *path)
{
  //scan
  vector<double> buf;
  ifstream fin(path);
  double t;
  while(fin>>t) buf.push_back(t);
  fin.close();
  
  //alloc
  njacks=buf.size()/clust_size;
  size=njacks*clust_size;
  data=(double*)malloc(sizeof(double)*size);
  for(int i=0;i<size;i++) data[i]=buf[i];
  
  cout<<"Finished reading"<<endl;
}

#define EXCLUDING_LOOP(i,j,ijack)					\
  for(int imin[2]={0,(ijack+1)*clust_size},imax[2]={ijack*clust_size,size},iter=0;iter<2;iter++) \
    for(int j=imin[iter],i=j-iter*clust_size;j<imax[iter];j++,i++)

void autocorr(double *ave_corr,double *jck_corr,double *err_corr)
{
  fftw_complex *in=(fftw_complex*)fftw_malloc(jackniffed_size*sizeof(fftw_complex));
  fftw_complex *out=(fftw_complex*)fftw_malloc(jackniffed_size*sizeof(fftw_complex));
  
  for(int i=0;i<jackniffed_size;i++) ave_corr[i]=err_corr[i]=0;
  
  for(int ijack=0;ijack<njacks;ijack++)
    {
      //compute ave
      double ave=0;
      EXCLUDING_LOOP(i,j,ijack) ave+=data[j];
      ave/=jackniffed_size;
      
      //copy in
      EXCLUDING_LOOP(i,j,ijack)
	{
	  in[i][0]=data[j]-ave;
	  in[i][1]=0;
	}
      
      //take fftw
      fftw_plan pf=fftw_plan_dft_1d(jackniffed_size,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
      fftw_execute(pf);
      fftw_destroy_plan(pf);
  
      //take module
      EXCLUDING_LOOP(i,j,ijack)
	{
	  in[i][0]=out[i][0]*out[i][0]+out[i][1]*out[i][1];
	  in[i][1]=0;
	}

      //take fftw
      fftw_plan pb=fftw_plan_dft_1d(jackniffed_size,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(pb);
      fftw_destroy_plan(pb);
      
      //copy back
      EXCLUDING_LOOP(i,j,ijack)
	{
	  double x=out[i][0]/out[0][0];
	  jck_corr[ijack*jackniffed_size+i]=x;
	  ave_corr[i]+=x;
	  err_corr[i]+=x*x;
	}
    }
  
  //convert to ave and err
  for(int i=0;i<jackniffed_size;i++)
    {
      ave_corr[i]/=njacks;
      err_corr[i]/=njacks;
      err_corr[i]-=ave_corr[i]*ave_corr[i];
      err_corr[i]=sqrt(err_corr[i]*(njacks-1));
    }
  
  fftw_free(in);
  fftw_free(out);
}

int main(int narg,char **arg)
{
  if(narg<2)
    {
      cerr<<"Use: "<<arg[0]<<" filein"<<endl;
      exit(1);
    }
  
  //load
  read(arg[1]);
  jackniffed_size=size-clust_size;
  
  //compute autocorr
  double *ave_corr=(double*)malloc(sizeof(double)*jackniffed_size);
  double *jck_corr=(double*)malloc(sizeof(double)*jackniffed_size*jackniffed_size);
  double *err_corr=(double*)malloc(sizeof(double)*jackniffed_size);
  autocorr(ave_corr,jck_corr,err_corr);
  
  //compute ave and non adjusted err
  double ave=0,err=0;
  for(int i=0;i<size;i++)
    {
      double x=data[i];
      ave+=x;
      err+=x*x;
    }
  ave/=size;
  err/=size;
  err-=ave*ave;
  err=sqrt(err/(size-1));
  
  //fix where to stop and plot autocorr
  int istop=0;
  ofstream autocorr_plot("/tmp/autocorr.xmg");
  autocorr_plot<<"@type xydy"<<endl;
  do
    {
      autocorr_plot<<istop<<" "<<ave_corr[istop]<<" "<<err_corr[istop]<<endl;
      istop++;
    }
  while(fabs(ave_corr[istop])>err_corr[istop] && istop<jackniffed_size-1);
  
  //compute tint across jacknives
  double stint=0,s2tint=0;
  for(int ijack=0;ijack<jackniffed_size;ijack++)
    {
      double tint=0;
      int i=1;
      while(fabs(jck_corr[ijack*jackniffed_size+i])>err_corr[i] && i<jackniffed_size-1)
	{
	  tint+=(jck_corr[ijack*jackniffed_size+i]+jck_corr[ijack*jackniffed_size+i-1])/2;
	  i++;
	}
      stint+=tint;
      s2tint+=tint*tint;
    }
  s2tint/=jackniffed_size;
  stint/=jackniffed_size;
  s2tint-=stint*stint;
  s2tint=sqrt(s2tint*(jackniffed_size-1));
  
  double tau=(2*stint-1)/2;
  cerr<<"tint: "<<stint<<" "<<s2tint<<endl;
  cout<<ave<<" +- "<<err*sqrt(2*tau+1)<<endl;
  
  free(data);
  free(ave_corr);
  free(jck_corr);
  free(err_corr);
  
  return 0;
}
