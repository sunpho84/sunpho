#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdarg.h>
#include <fftw3.h>
#include <math.h>

int autocorr_debug=true;

const long unsigned int max_njacks=2048;

using namespace std;

class autocorr_data_t
{
public:
  double *data;
  void autocorr(double *ave_corr,double *jck_corr,double *err_corr);
  void compute_tint(double &med_tint,double &err_tint,const char *path);
  long unsigned int clust_size,size,int_size;
  bool allocated;
  int jackniffed_size;
  int njacks;
  autocorr_data_t(long unsigned int clust_size=1):clust_size(clust_size),size(0),allocated(false){}
  ~autocorr_data_t(){if(allocated) delete[] data;}
  double &operator[](int i){return data[i];}
  void ave_err(double &ave,double &err);
  void get_from(vector<double> &ext)
  {
    size=ext.size();
    data=new double[size];
    allocated=true;
    for(long unsigned int i=0;i<size;i++) data[i]=ext[i];
  }
  void point_to(vector<double> &ext)
  {
    data=&(ext[0]);
    size=ext.size();
  }
};

void autocorr_data_t::ave_err(double &ave,double &err)
{
  //compute ave and non adjusted err
  ave=err=0;
  for(long unsigned int i=0;i<size;i++)
    {
      ave+=data[i];
      err+=data[i]*data[i];
    }
  ave/=size;
  err/=size;
  err-=ave*ave;
  err=sqrt(err/(size-1));
}

#define EXCLUDING_LOOP(i,j,ijack)					\
  for(long unsigned int imin[2]={0,(ijack+1)*clust_size},imax[2]={ijack*clust_size,size},iter=0;iter<2;iter++) \
    for(long unsigned int j=imin[iter],i=j-iter*clust_size;j<imax[iter];j++,i++)

void autocorr_data_t::autocorr(double *ave_corr,double *jck_corr,double *err_corr)
{
  //allocate fftw
  fftw_complex *in=(fftw_complex*)fftw_malloc(jackniffed_size*sizeof(fftw_complex));
  fftw_complex *out=(fftw_complex*)fftw_malloc(jackniffed_size*sizeof(fftw_complex));
  fftw_plan pf=fftw_plan_dft_1d(jackniffed_size,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_plan pb=fftw_plan_dft_1d(jackniffed_size,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
  
  for(int i=0;i<jackniffed_size;i++) ave_corr[i]=err_corr[i]=0;
  
  for(int ijack=0;ijack<njacks;ijack++)
    {
      cerr<<ijack<<"/"<<njacks<<endl;
      
      //compute ave
      double ave=0;
      EXCLUDING_LOOP(i,j,ijack) ave+=(*this)[j];
      ave/=jackniffed_size;
      
      //copy in
      EXCLUDING_LOOP(i,j,ijack)
	{
	  in[i][0]=(*this)[j]-ave;
	  in[i][1]=0;
	}
      
      //take fftw
      fftw_execute(pf);
    
      //take module
      for(int i=0;i<jackniffed_size;i++)
	{
	  in[i][0]=out[i][0]*out[i][0]+out[i][1]*out[i][1];
	  in[i][1]=0;
	}

      //take fftw
      fftw_execute(pb);
      
      //copy back
      for(int i=0;i<jackniffed_size;i++)
	{
	  double x=out[i][0]/out[0][0];
	  jck_corr[ijack*jackniffed_size+i]=x;
	  ave_corr[i]+=x;
	  err_corr[i]+=x*x;
	}
      
      if(ijack==0)
	{
	  ofstream autocorr_plot("/tmp/autocorr.xmg");
	  autocorr_plot<<"@type xy"<<endl;
	  int i=0;
	  do
	    {
	      autocorr_plot<<i<<" "<<out[i][0]/out[0][0]<<endl;
	      i++;
	    }
	  while(out[i][0]>0);
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
  fftw_destroy_plan(pf);
  fftw_destroy_plan(pb);
}

//compute tint and err
void autocorr_data_t::compute_tint(double &med_tint,double &err_tint,const char *path=NULL)
{
  //adjust clust_size so to reach consistency
  bool loop;
  do
    {
      //enforce not to overshot
      long unsigned int max_clust_size=size/max_njacks;
      if(clust_size<max_clust_size)
	{
	  cout<<"Too many jacknives to make clust_size="<<clust_size<<", reducing to "<<max_njacks<<" jacknives"<<endl;
	  clust_size=max_clust_size;
	}
      
      //fix njacks according
      njacks=size/clust_size;
      int_size=njacks*clust_size;
      jackniffed_size=int_size-clust_size;
    
      cout<<" clust_size: "<<clust_size<<endl;
      cout<<" njacks: "<<njacks<<endl;
      cout<<" int_size: "<<int_size<<endl;
      cout<<" jackniffed_size: "<<jackniffed_size<<endl;
      
      //compute autocorr
      double *ave_corr=(double*)malloc(sizeof(double)*jackniffed_size);
      double *jck_corr=(double*)malloc(sizeof(double)*jackniffed_size*njacks);
      double *err_corr=(double*)malloc(sizeof(double)*jackniffed_size);
      autocorr(ave_corr,jck_corr,err_corr);
      
      //fix where to stop and plot autocorr
      int istop=0;
      ofstream autocorr_plot;
      if(path!=NULL) autocorr_plot.open(path);
      autocorr_plot<<"@type xydy"<<endl;
      do
	{
	  if(path!=NULL) autocorr_plot<<istop<<" "<<ave_corr[istop]<<" "<<err_corr[istop]<<endl;
	  istop++;
	}
      while(fabs(ave_corr[istop])>0.5*err_corr[istop] && istop<jackniffed_size-1);
      
      //compute tint across jacknives
      med_tint=0,err_tint=0;
      for(int ijack=0;ijack<njacks;ijack++)
	{
	  double tint=0;
	  int i=1;
	  do
	    {
	      tint+=(jck_corr[ijack*jackniffed_size+i]+jck_corr[ijack*jackniffed_size+i-1])/2;
	      i++;
	    }
	  while(fabs(jck_corr[ijack*jackniffed_size+i])>0.5*err_corr[i] && i<jackniffed_size-2);
	  med_tint+=tint;
	  err_tint+=tint*tint;
	}
      err_tint/=njacks;
      med_tint/=njacks;
      err_tint-=med_tint*med_tint;
      err_tint=sqrt(err_tint*(njacks-1)); 
      if(autocorr_debug) cout<<"tint: "<<med_tint<<" +- "<<err_tint<<endl;
      
      //free
      free(ave_corr);
      free(jck_corr);
      free(err_corr);
      
      loop=false;
      if(fabs(clust_size-2*med_tint)>=4*err_tint)
	{
	  long unsigned int new_clust_size=(long unsigned int)(2*med_tint+0.5);
	  if(new_clust_size!=clust_size && new_clust_size!=0 && new_clust_size>=max_clust_size && new_clust_size<size)
	    {
	      clust_size=new_clust_size;
	      loop=true;
	      if(autocorr_debug) cout<<" recomputing with cluster size: "<<clust_size<<endl;
	    }
	}
    }
  while(loop && med_tint>1);
}

void compute_tint(double &med_tint,double &err_tint,vector<double> &ext,int cluster_size=1,const char *path=NULL)
{
  autocorr_data_t aut(cluster_size);
  aut.point_to(ext);
  aut.compute_tint(med_tint,err_tint,path);
}
