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

const int clust_size=5;
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
}

#define EXCLUDING_LOOP(i,j,ijack)					\
  for(int imin[2]={0,(ijack+1)*clust_size},imax[2]={ijack*clust_size,size},iter=0;iter<2;iter++) \
    for(int j=imin[iter],i=j-iter*clust_size;j<imax[iter];j++,i++)

void autocorr(double *ave_corr,double *err_corr)
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
  double *err_corr=(double*)malloc(sizeof(double)*jackniffed_size);
  autocorr(ave_corr,err_corr);
  
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
  
  //compute tint
  int i=1;
  double tint=0;
  ofstream autocorr_plot("/tmp/autocorr.xmg");
  autocorr_plot<<"@type xydy"<<endl;
  autocorr_plot<<0<<" "<<1<<" "<<0<<endl;
  while(fabs(ave_corr[i])>err_corr[i])
    {
      tint+=(ave_corr[i]+ave_corr[i-1])/2;
      autocorr_plot<<i<<" "<<ave_corr[i]<<" "<<err_corr[i]<<endl;
      i++;
    }
  tint=2*tint+1;
  cerr<<"tint: "<<tint<<endl;
  cout<<ave<<" +- "<<err*sqrt(tint)<<endl;
  
  free(data);
  free(ave_corr);
  free(err_corr);
  
  read_from_file("/tmp/in");
  
  return 0;
}
