#pragma once

void ludcmp_real(double *a,int n,int *indx)
{
  int imax=0;
  double vv[n];
  
  for(int i=0;i<n;i++)
    {
      double big=0;
      for(int j=0;j<n;j++) big=max(big,fabs(a[i*n+j]));
      vv[i]=1/big;
    }
  for(int j=0;j<n;j++)
    {
      for(int i=0;i<j;i++)
        {
          double sum=a[i*n+j];
          for(int k=0;k<i;k++) sum-=a[i*n+k]*a[k*n+j];
          a[i*n+j]=sum;
        }
      double big=0;
      for(int i=j;i<n;i++)
	{
	  double sum=a[i*n+j];
	  for(int k=0;k<j;k++) sum-=a[i*n+k]*a[k*n+j];
	  a[i*n+j]=sum;
	  double dum=vv[i]*fabs(sum);
	  if(dum>=big)
	    {
	      big=dum;
	      imax=i;
	    }
	}
      if(j!=imax+1)
	{
	  for(int k=0;k<n;k++)
	    {
	      double dum=a[imax*n+k];
	      a[imax*n+k]=a[j*n+k];
	      a[j*n+k]=dum;
	    }
        vv[imax]=vv[j];
      }
      indx[j]=imax;

      if(j!=(n-1))
	{
	  double dum=1/(a[j*n+j]);
	  for(int i=j+1;i<n;i++) a[i*n+j]*=dum;
	}
    }
}

void lubksb_real(double *a,int n,int *indx,double *b)
{
  int ii=-1;
  
  for(int i=0;i<n;i++)
    { 
      int ip=indx[i];
      double sum=b[ip];
      b[ip]=b[i];
      if(ii!=-1) for(int j=ii;j<i;j++) sum-=a[i*n+j]*b[j];
      else if(sum) ii=i;
      b[i]=sum; 
    }
  for(int i=n-1;i>=0;i--)
    {
      double sum=b[i];
      for(int j=i+1;j<n;j++) sum-=a[i*n+j]*b[j];
      b[i]=sum/a[i*n+i];
    }
}

void invert_matrix_real(double *out,double *in,int n)
{
  double col[n],*fuffa=malloc(sizeof(double)*n*n);
  int indx[n];
  memcpy(fuffa,in,sizeof(double)*n*n);

  ludcmp_real(fuffa,n,indx);
  for(int j=0;j<n;j++)
    {
      for(int i=0;i<n;i++) col[i]=0;
      col[j]=1;
      lubksb_real(fuffa,n,indx,col);
      for(int i=0;i<n;i++) out[i*n+j]=col[i];
    }
  
  free(fuffa);
}

/////////////////// complex case ////////////

void ludcmp_complex(complex *a,int n,int *indx)
{
  int imax=0;
  double vv[n];
  
  for(int i=0;i<n;i++)
    {
      double big=0;
      for(int j=0;j<n;j++) big=max(big,complex_fabs(a[i*n+j]));
      vv[i]=1/big;
    }
  for(int j=0;j<n;j++)
    {
      for(int i=0;i<j;i++)
        {
          complex sum;
	  complex_copy(sum,a[i*n+j]);
          for(int k=0;k<i;k++) complex_subtassign_the_prod(sum,a[i*n+k],a[k*n+j]);
          complex_copy(a[i*n+j],sum);
        }
      double big=0;
      for(int i=j;i<n;i++)
	{
	  complex sum;
	  complex_copy(sum,a[i*n+j]);
	  for(int k=0;k<j;k++) complex_subtassign_the_prod(sum,a[i*n+k],a[k*n+j]);
	  complex_copy(a[i*n+j],sum);
	  double dum=vv[i]*complex_fabs(sum);
	  if(dum>=big)
	    {
	      big=dum;
	      imax=i;
	    }
	}
      if(j!=imax+1)
	{
	  for(int k=0;k<n;k++)
	    {
	      complex dum;
	      complex_copy(dum,a[imax*n+k]);
	      complex_copy(a[imax*n+k],a[j*n+k]);
	      complex_copy(a[j*n+k],dum);
	    }
        vv[imax]=vv[j];
      }
      indx[j]=imax;

      if(j!=(n-1))
	{
	  complex dum;
	  complex_reciprocal(dum,(a[j*n+j]));
	  for(int i=j+1;i<n;i++) complex_prodassign(a[i*n+j],dum);
	}
    }
}

void lubksb_complex(complex *a,int n,int *indx,complex *b)
{
  int ii=-1;
  
  for(int i=0;i<n;i++)
    { 
      int ip=indx[i];
      complex sum;
      complex_copy(sum,b[ip]);
      complex_copy(b[ip],b[i]);
      if(ii!=-1) for(int j=ii;j<i;j++) complex_subtassign_the_prod(sum,a[i*n+j],b[j]);
      else if(sum[0]!=0||sum[1]!=0) ii=i;
      complex_copy(b[i],sum); 
    }
  for(int i=n-1;i>=0;i--)
    {
      complex sum;
      complex_copy(sum,b[i]);
      for(int j=i+1;j<n;j++) complex_subtassign_the_prod(sum,a[i*n+j],b[j]);
      complex rec;
      complex_reciprocal(rec,a[i*n+i]);
      complex_prod(b[i],sum,rec);
    }
}

void invert_matrix_complex(complex *out,complex *in,int n)
{
  complex col[n],*fuffa=malloc(sizeof(complex)*n*n);
  int indx[n];
  memcpy(fuffa,in,sizeof(complex)*n*n);

  ludcmp_complex(fuffa,n,indx);
  for(int j=0;j<n;j++)
    {
      memset(col,0,sizeof(complex)*n);
      col[j][0]=1;
      lubksb_complex(fuffa,n,indx,col);
      for(int i=0;i<n;i++) complex_copy(out[i*n+j],col[i]);
    }
  
  free(fuffa);
}

//////////////////////////////////////////////////////////////

void invert_ccss_propagator(ccss_propagator out,ccss_propagator in)
{
  for(int imom=0;imom<nmom;imom++)
    invert_matrix_complex((complex*)out[imom],(complex*)in,12);
}
