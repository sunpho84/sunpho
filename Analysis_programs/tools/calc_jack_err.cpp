#include "include.h"

const int njacks=3;
const int nboots=2000;
const int nsamples=1000;

mt19937_64 gen(5543234);

void boot(double &m,double &e,double *data)
{
  uniform_int_distribution<int> clust_id_dis(0,njacks-1);
  
  double sx=0,s2x=0;
  for(int iboot=0;iboot<nboots;iboot++)
    {
      double x=data[clust_id_dis(gen)];

      sx+=x;
      s2x+=x*x;
    }
  sx/=nboots;
  s2x/=nboots;
  
  m=sx;
  e=sqrt((s2x-sx*sx)/(njacks-1));
}

void naive(double &m,double &e,double *data)
{
  double sx=0,s2x=0;
  for(int ijack=0;ijack<njacks;ijack++)
    {
      double x=data[ijack];
      sx+=x;
      s2x+=x*x;
    }
  sx/=njacks;
  s2x/=njacks;

  m=sx;
  e=sqrt((s2x-sx*sx)/(njacks-1));
}

void jack(double &m,double &e,double *data)
{
  double sum=0;
  for(int ijack=0;ijack<njacks;ijack++) sum+=data[ijack];

  double sx=0,s2x=0;
  for(int ijack=0;ijack<njacks;ijack++)
    {
      double x=(sum-data[ijack])/(njacks-1);
      sx+=x;
      s2x+=x*x;
    }
  sx/=njacks;
  s2x/=njacks;
  
  m=sx;
  e=sqrt((s2x-sx*sx)*(njacks-1));
}

void print_ave_err(void (*fun)(double&,double&,double*data),double *data,const char *comm)
{
  double m,e;
  fun(m,e,data);
  cout<<comm<<" estimators:\t"<<m<<" "<<e<<endl;
}

void compare_err(void (*fun1)(double&,double&,double*data),void (*fun2)(double&,double&,double*data),double *data)
{
    double m1,e1;
    fun1(m1,e1,data);
    double m2,e2;
    fun2(m2,e2,data);
    cout<<e1-e2<<endl;
}

int main()
{
  srand(101);

    for(int isa=0;isa<nsamples;isa++)
      {
	
	//prepare the sample 
	double data[njacks];
	for(int i=0;i<njacks;i++) data[i]=(double)rand()/RAND_MAX;
	
	compare_err(naive,boot,data);
	
	//print_ave_err(naive,data,"naive");
	//print_ave_err(jack,data,"jack");
	//print_ave_err(boot,data,"boot");
      }
  return 0;
}
