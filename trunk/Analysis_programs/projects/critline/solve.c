#include <stdio.h>
#include <math.h>

int main()
{
  const int d=2;
  double A[d*d]={0,1,1,0},b[d]={2,3},x[d];

  for(int i=0;i<d;i++)
    {
      double C=A[i*d+i];
      for(int j=i;j<d;j++) A[i*d+j]/=C;
      b[i]/=C;

      for(int k=i+1;k<d;k++)
	{
	  double C=A[k*d+i];
	  for(int j=i;j<d;j++) A[k*d+j]-=A[i*d+j]*C;
	  b[k]-=C*b[i];
	}
    }
  
  for(int k=d-1;k>=0;k--)
    {
      double S=0;
      for(int i=k+1;i<d;i++) S+=A[k*d+i]*x[i];
      x[k]=b[k]-S;
    } 

  for(int i=0;i<d;i++) printf("x[%d] = %lf\n",i,x[i]);
  
  return 0;
}
