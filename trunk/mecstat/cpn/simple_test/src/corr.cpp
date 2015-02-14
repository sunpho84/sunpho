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
      P_fftw=new dcomplex[2*N*N*L];
      {
	int rank=1,sizes[1]={L};
	int howmany=2*N*N,ostride=howmany,istride=howmany;
	int idist=1,odist=1;
	fw=fftw_plan_many_dft(rank,sizes,howmany,
			      (fftw_complex*)P_fftw,NULL,istride,idist,
			      (fftw_complex*)P_fftw,NULL,ostride,odist,
			      FFTW_FORWARD,FFTW_MEASURE);
      }
      
      //backward
      C_fftw=new dcomplex[2*L];
      {
	int rank=1,sizes[1]={L};
	int howmany=2,ostride=howmany,istride=howmany;
	int idist=1,odist=1;
	bw=fftw_plan_many_dft(rank,sizes,howmany,
			      (fftw_complex*)C_fftw,NULL,istride,idist,
			      (fftw_complex*)C_fftw,NULL,ostride,odist,
			      FFTW_BACKWARD,FFTW_MEASURE);
      }
      init_fftw_flag=true;
    }
}

//compute correlation function
void compute_corr(double &mag0,double &mag1,double &mom2,double *out,double *outd,dcomplex *z)
{
  init_fftw();
  
#pragma omp parallel for
  for(int i=0;i<2*N*N*L;i++) P_fftw[i]=0;

#pragma omp parallel for
  for(int t=0;t<L;t++)
    for(int x=0;x<L;x++)
      {
	int s=t*L+x;
	int tmx=(t-x+L)%L;
	
	for(int i=0;i<N;i++)
	  for(int j=0;j<N;j++)
	    {
	      dcomplex c=conj(z[s*N+i])*z[s*N+j];
	      P_fftw[((t*N+i)*N+j)*2+0]+=c;
	      P_fftw[((tmx*N+i)*N+j)*2+1]+=c;
	    }
      }
  
  fftw_execute(fw);
  
#pragma omp parallel for
  for(int t=0;t<L;t++)
    {
      C_fftw[0+2*t]=0;
      C_fftw[1+2*t]=0;
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  {
	    C_fftw[0+2*t]+=norm(P_fftw[((t*N+i)*N+j)*2+0]);
	    C_fftw[1+2*t]+=norm(P_fftw[((t*N+i)*N+j)*2+1]);
	  }
      C_fftw[0+2*t]/=V*L;
      C_fftw[1+2*t]/=V*L;
    }
  
  //subtract disconnected
  C_fftw[0]-=L/(double)N;
  C_fftw[1]-=L/(double)N;
  
  //note magnetization
  mag0=C_fftw[0].real()*L;
  mag1=C_fftw[2].real()*L;
  
  fftw_execute(bw);
  
  for(int t=0;t<L;t++)
    {
      out[t]=C_fftw[0+2*t].real();
      outd[t]=C_fftw[1+2*t].real()*sqrt(2);
    }
}

