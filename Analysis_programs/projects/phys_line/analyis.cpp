#include <math.h>
#include <iostream>

using namespace std;

double db0=11.972809241766459,f0=0.26772510275074468;
double r0=2.1975;
double l3=1.0423672088359199,l4=2.1490196003108331;
double xi=db0/(16*M_PI*M_PI*f0*f0);
double ml_phys=3.6e-3;

template <class T> T m_fun(T m)
{
  m*=r0;
  
  T m2=db0*m*(1+  xi*m*(log(db0*m)+l3));
  
  return sqrt(m2)/r0;
}

template <class T> T f_fun(T m)
{
  m*=r0;
  
  T f=   f0* (1-2*xi*m*(log(db0*m)-l4));
  
  return f/r0;
}

template <class T> T ratio(T m)
{return m_fun(m)/f_fun(m);}


double find(double x)
{
  double m=1e-5;
  double s=1e-5;
  
  do
    if(ratio(m+s)<=x)
      {
	m+=s;
	s*=2;
      }
    else s/=2;
  while(s>1e-14);
  
  return m;
}

int main()
{
  for(double i=0;i<5;i+=0.5)
    {
      double m=ml_phys*i;
      cout<<m<<" "<<ratio(m)<<" "<<m_fun(m)<<endl;
    }

  double ampi,afpi;
  cout<<"ampi? ";
  cin>>ampi;
  cout<<"afpi? ";
  cin>>afpi;
  double ml=find(ampi/afpi);
  double a=ampi/m_fun(ml);
  cout<<"ml: "<<ml*1000<<" MeV"<<endl;
  //cout<<"a: "<<a<<" GeV^-1"<<endl;
  cout<<"a^-1: "<<1/a<<" GeV"<<endl;
  cout<<"MPi: "<<ampi/a<<" GeV"<<endl;
  
  return 0;
}
