#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "data.hpp"
#include "geometry.hpp"
#include "parameters.hpp"
#include "types.hpp"

#include <iostream>
#include <fftw3.h>
#include <omp.h>

using namespace std;

dcomplex *P_fftw,*C_fftw;
bool init_fftw_flag=false;
fftw_plan fw,bw;

void init_fftw()
{
  if(init_fftw_flag==false)
    {
      fftw_init_threads();
      fftw_plan_with_nthreads(omp_get_max_threads());
      
      //forward
      P_fftw=new dcomplex[N*N*V];
      int rank=2,sizes[2]={L,L};
      int howmany=N*N,ostride=howmany,istride=howmany;
      int idist=1,odist=1;
      fw=fftw_plan_many_dft(rank,sizes,howmany,
			    (fftw_complex*)P_fftw,NULL,istride,idist,
			    (fftw_complex*)P_fftw,NULL,ostride,odist,
			    FFTW_FORWARD,FFTW_MEASURE);
      
      //backward
      C_fftw=new dcomplex[V];
      bw=fftw_plan_dft(rank,sizes,(fftw_complex*)C_fftw,(fftw_complex*)C_fftw,FFTW_BACKWARD,FFTW_MEASURE);
      
      init_fftw_flag=true;
    }
}

void stop_fftw()
{
  if(init_fftw_flag==true)
    {
      delete[] P_fftw;
      delete[] C_fftw;
      
      fftw_destroy_plan(fw);
      fftw_destroy_plan(bw);
      init_fftw_flag=false;
    }
}

//compute correlation function
void compute_corr(double &mag0,double &mag1,double &mom2,double *out,double *outd,dcomplex *z)
{
  init_fftw();

  //fill P
  for(int s=0;s<V;s++)
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
	P_fftw[(s*N+i)*N+j]=conj(z[s*N+i])*z[s*N+j];
  
  //take FFT
  fftw_execute(fw);
  
  //take P(k)*P(-k)
  for(int s=0;s<V;s++)
    {
      //find -k
      coords c,c1;
      coords_of_site(c,s);
      c1[0]=(L-c[0])%L;
      c1[1]=(L-c[1])%L;
      int s1=site_of_coords(c1);
      
      //take the product
      C_fftw[s]=0;
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  C_fftw[s]+=P_fftw[(s*N+i)*N+j]*P_fftw[(s1*N+j)*N+i];
      C_fftw[s]/=V;
    }
  C_fftw[0]-=V/(double)N;
  
  //mark magnetic susceptibility ingredients
  mag0=C_fftw[0].real();
  mag1=(C_fftw[1].real()+C_fftw[L].real())/2;
  
  //take anti-fftw
  fftw_execute(bw);
  
  //compute second momentum
  mom2=0;
  for(int dx=-L/2;dx<L/2;dx++)
    for(int dy=-L/2;dy<L/2;dy++)
      {
	coords c={dx+L/2,dy+L/2};
	int x2=dx*dx+dy*dy;
	mom2+=x2*C_fftw[site_of_coords(c)].real();
      }
  
  //project to zero spatial momentum
  for(int x=0;x<L;x++)
    {
      double res=0,resd=0;
      for(int y=0;y<L;y++)
	{
	  coords c1={x,y},c2={y,x};
	  res+=C_fftw[site_of_coords(c1)].real()+C_fftw[site_of_coords(c2)].real();
	  coords c3={(x-y+L)%L,y%L},c4={(x+y)%L,y%L};
	  resd+=C_fftw[site_of_coords(c3)].real()+C_fftw[site_of_coords(c4)].real();
	}
      out[x]=res/V/2;
      outd[x]=sqrt(2)*resd/V/2;
    }
    
  //stop_fftw();
}

//compute correlation function
void compute_corre(double *out,double *outd,dcomplex *z)
{
  for(int x=0;x<L;x++)
    {
      double res=0,resd=0;
#pragma omp parallel for reduction(+:res,resd)
      for(int y=0;y<L;y++)
	for(int s0=0;s0<V;s0++)
	  {
	    coords c0;
	    coords_of_site(c0,s0);
	    
	    coords c={(x+c0[0])%L,(y+c0[1])%L};
	    int s=site_of_coords(c);
	    
	    dcomplex t=0.0;
	    for(int n=0;n<N;n++) t+=conj(z[s*N+n])*z[s0*N+n];
	    res+=norm(t);

	    coords c1={(x-y+L+c0[0])%L,(y+c0[1])%L};
	    int s1=site_of_coords(c1);
	    
	    dcomplex t1=0.0;
	    for(int n=0;n<N;n++) t1+=conj(z[s1*N+n])*z[s0*N+n];
	    resd+=norm(t1);
	  }
      out[x]=res/V-L/2;
      outd[x]=sqrt(2)*(resd/V-L/2);
    }
}
