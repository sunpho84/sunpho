#include "src/global.c"

void matrix_prod_complex(complex *a,complex *b,complex *c,int n)
{
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      {
	a[i*n+j][0]=a[i*n+j][1]=0;
	for(int k=0;k<n;k++)
	  complex_summassign_the_prod(a[i*n+j],b[i*n+k],c[k*n+j]);
      }  
}

int main()
{
  L=24;
  T=48;
  int njack=10;
  int nmass=5;
  init();
  
  ccss_propagator ***in=malloc(sizeof(ccss_propagator**)*(njack+1));
  for(int ijack=0;ijack<=njack;ijack++)
    {
      in[ijack]=malloc(sizeof(ccss_propagator*)*2);
      for(int r=0;r<2;r++) in[ijack][r]=malloc(sizeof(ccss_propagator)*nmass);
    }
  
  read_ccss_propagator_set(in,"/home/francesco/QCD/LAVORI/RI-MOM/test/%04d/fft",501,30,2,njack,nmass);
  
  cmom *a=malloc((njack+1)*sizeof(cmom));
  for(int ijack=0;ijack<=njack;ijack++) ccss_trace_ccss_propagator(a[ijack],in[ijack][0][0]);
  for(int imom=0;imom<nmom;imom++) printf("%lg %lg\n",a[0][imom][0],a[0][imom][1]);
  
  return 0;
}
