#include <algorithm>
#include <cmath>

//check if we can use c++11 features
#if __cplusplus > 199711L
 #define RANDOM_CPP11
#endif

#ifdef RANDOM_CPP11
 #include <random>
#else
 #include <stdlib.h>
#endif

#ifdef RANDOM_CPP11
std::mt19937_64 gen;
#endif
void random_init(int seed)
{
#ifdef RANDOM_CPP11
  gen.seed(seed);
#else
  srand(seed);
#endif
}

//generate a number distributed uniformly over the range [min,max)
double get_unif_double(double min,double max)
{
#ifdef RANDOM_CPP11
  std::uniform_real_distribution<double> distribution(min,max);
  return distribution(gen);
#else
  int i;
  do i=rand();while(i==RAND_MAX);
  return min+(double)rand()*(max-min)/RAND_MAX;
#endif
}

//like previously but on integer
int get_unif_int(int min,int max)
{
  if(min==max) return min;
#ifdef RANDOM_CPP11
  else
    {
      std::uniform_int_distribution<int> distribution(min,max-1);
      return distribution(gen);
    }
#else
  return (int)get_unif_double(min,max);
#endif
}

//extract a random number distributed according to a gaussian
double get_gauss_double(double ave,double sig)
{
#ifdef RANDOM_CPP11
  std::normal_distribution<double> distribution(ave,sig);
  return distribution(gen);
#else
  double q,r,x;
  static double y;
  static bool flag=true;
  
  if(flag)
    {
      r=sqrt(-2*log(1-get_unif_double(0,1)));
      q=2*M_PI*get_unif_double(0,1);
      
      x=r*cos(q);
      y=r*sin(q);
    }
  else x=y;
  
  flag=!flag;

  return x;
#endif
}
//fill a source with random gauss
void fill_source(double *b,double ave,double sigma,int n)
{for(int i=0;i<n;i++) b[i]=get_gauss_double(ave,sigma);}

//fill a unitary matrix to be used as base for rotation
void fill_unitary_matrix(double *u,int n)
{
  for(int i=0;i<n;i++)
    {
      //fill current line
      for(int j=0;j<n;j++) u[i*n+j]=get_gauss_double(0,1);
      //othonormalize w.r.t. past lines
      for(int j=0;j<i;j++) 
        {
          double norm=0;
          for(int k=0;k<n;k++) norm+=u[i*n+k]*u[j*n+k];
          for(int k=0;k<n;k++) u[i*n+k]-=u[j*n+k]*norm;
        }
      //normalize current line
      double norm=0;
      for(int k=0;k<n;k++) norm+=u[i*n+k]*u[i*n+k];
      norm=1/sqrt(norm);
      for(int k=0;k<n;k++) u[i*n+k]*=norm;
    }
}

//fill the matrix
void fill_defpos_symm_matrix(double *m,double cond_numb,int n)
{
  //fill a unitary matrix
  double *u=new double[n*n];
  fill_unitary_matrix(u,n);
  
  //prepare eigenvalues
  double *e=new double[n];
  e[0]=1;
  for(int i=1;i<n-1;i++) e[i]=exp(-get_unif_double(0,log(cond_numb)));
  e[n-1]=1/cond_numb;
  //make some swaps
  for(int iswap=0;iswap<3*n;iswap++) std::swap(e[get_unif_int(0,n)],e[get_unif_int(0,n)]);
  
  //rotate
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++) 
      {
        m[i*n+j]=0;
        for(int k=0;k<n;k++) m[i*n+j]+=u[k*n+i]*e[k]*u[k*n+j];
      }
  
  delete[] e;
  delete[] u;
}

//fill a completely random matrix with gaussian entries
void fill_gauss_matrix(double *m,double ave,double sigma,int n)
{
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      m[i*n+j]=get_gauss_double(ave,sigma);
}
