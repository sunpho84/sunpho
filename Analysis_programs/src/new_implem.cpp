#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

using namespace std;

const int njacks=16;
const int nmass=10;
const int T=48;

class jack : public vector<double>
{
public:
  void med_err(double &ym,double &y2m);
};

class jack_corr : public vector<jack>
{
public:
  jack_corr symmetrized(int pm);
};

jack operator+(const jack &a,const jack &b){jack c;for(int i=0;i<a.size();i++) c.push_back(a[i]+b[i]);return c;}
jack operator-(const jack &a,const jack &b){jack c;for(int i=0;i<a.size();i++) c.push_back(a[i]-b[i]);return c;}
jack operator*(const jack &a,const jack &b){jack c;for(int i=0;i<a.size();i++) c.push_back(a[i]*b[i]);return c;}
jack operator/(const jack &a,const jack &b){jack c;for(int i=0;i<a.size();i++) c.push_back(a[i]/b[i]);return c;}
jack_corr operator+(const jack_corr &a,const jack_corr &b)
{jack_corr c;for(int i=0;i<a.size();i++) c.push_back(a[i]+b[i]);return c;}
jack_corr operator-(const jack_corr &a,const jack_corr &b)
{jack_corr c;for(int i=0;i<a.size();i++) c.push_back(a[i]-b[i]);return c;}
jack_corr operator*(const jack_corr &a,const jack_corr &b)
{jack_corr c;for(int i=0;i<a.size();i++) c.push_back(a[i]*b[i]);return c;}
jack_corr operator/(const jack_corr &a,const jack_corr &b)
{jack_corr c;for(int i=0;i<a.size();i++) c.push_back(a[i]/b[i]);return c;}

jack operator+(const jack &a,double b){jack c;for(int i=0;i<a.size();i++) c.push_back(a[i]+b);return c;}
jack operator-(const jack &a,double b){jack c;for(int i=0;i<a.size();i++) c.push_back(a[i]-b);return c;}
jack operator*(const jack &a,double b){jack c;for(int i=0;i<a.size();i++) c.push_back(a[i]*b);return c;}
jack operator/(const jack &a,double b){jack c;for(int i=0;i<a.size();i++) c.push_back(a[i]/b);return c;}
jack_corr operator+(const jack_corr &a,double b)
{jack_corr c;for(int i=0;i<a.size();i++) c.push_back(a[i]+b);return c;}
jack_corr operator-(const jack_corr &a,double b)
{jack_corr c;for(int i=0;i<a.size();i++) c.push_back(a[i]-b);return c;}
jack_corr operator*(const jack_corr &a,double b)
{jack_corr c;for(int i=0;i<a.size();i++) c.push_back(a[i]*b);return c;}
jack_corr operator/(const jack_corr &a,double b)
{jack_corr c;for(int i=0;i<a.size();i++) c.push_back(a[i]/b);return c;}

jack operator+(double a,const jack &b){jack c;for(int i=0;i<b.size();i++) c.push_back(a+b[i]);return c;}
jack operator-(double a,const jack &b){jack c;for(int i=0;i<b.size();i++) c.push_back(a-b[i]);return c;}
jack operator*(double a,const jack &b){jack c;for(int i=0;i<b.size();i++) c.push_back(a*b[i]);return c;}
jack operator/(double a,const jack &b){jack c;for(int i=0;i<b.size();i++) c.push_back(a/b[i]);return c;}
jack_corr operator+(double a,const jack_corr &b)
{jack_corr c;for(int i=0;i<b.size();i++) c.push_back(a+b[i]);return c;}
jack_corr operator-(double a,const jack_corr &b)
{jack_corr c;for(int i=0;i<b.size();i++) c.push_back(a-b[i]);return c;}
jack_corr operator*(double a,const jack_corr &b)
{jack_corr c;for(int i=0;i<b.size();i++) c.push_back(a*b[i]);return c;}
jack_corr operator/(double a,const jack_corr &b)
{jack_corr c;for(int i=0;i<b.size();i++) c.push_back(a/b[i]);return c;}

template <class T> T operator+=(T &a,const T &b){return a=a+b;}
template <class T> T operator-=(T &a,const T &b){return a=a-b;}
template <class T> T operator*=(T &a,const T &b){return a=a*b;}
template <class T> T operator/=(T &a,const T &b){return a=a/b;}

template <class T> T sqr(const T &a){return a*a;}

template <class T> T operator+=( T &a,const double b){return a=a+b;}
template <class T> T operator-=( T &a,const double b){return a=a-b;}
template <class T> T operator*=( T &a,const double b){return a=a*b;}
template <class T> T operator/=( T &a,const double b){return a=a/b;}

void crash(const char *temp,...)
{
  char buffer[1024];
  va_list args;

  va_start(args,temp);
  vsprintf(buffer,temp,args);
  va_end(args);

  cerr<<"ERROR: "<<buffer<<endl;
  exit(1);
}

void jack::med_err(double &ym,double &y2m)
{
  ym=y2m=0;
  int n=this->size()-1;
  for(int i=0;i<n;i++)
    {
      double y=(*this)[i];
      ym+=y;
      y2m+=y*y;
    }
  y2m/=n;
  ym/=n;
  y2m=sqrt((y2m-ym*ym)*n/(n-1));
  ym=(*this)[n];
}
  
jack logj(const jack &a)
{jack c;c.resize(a.size());transform(a.begin(),a.end(),c.begin(),logl);return c;}

ifstream& operator>>(ifstream &in,jack &out)
{
  out.clear();
  for(int ijack=0;ijack<=njacks;ijack++)
    {
      double temp;
      in>>temp;
      out.push_back(temp);
    }
  return in;
}

ostream& operator<<(ostream &out,jack in)
{
  double m,e;
  in.med_err(m,e);
  
  out<<m<<" "<<e;
  
  return out;
}

/////

jack_corr jack_corr::symmetrized(int pm=1)
{
  jack_corr out;
  int T=this->size();
  for(int t1=0;t1<=T/2;t1++)
    {
      int t2=(T-t1)%T;
      out.push_back(((*this)[t1]+pm*(*this)[t2])/2);
    }
  return out;
}

ifstream& operator>>(ifstream &in,jack_corr &out)
{
  out.clear();
  for(int t=0;t<T;t++)
    {
      jack temp;
      in>>temp;
      out.push_back(temp);
    }
  return in;
}

ostream& operator<<(ostream &out,jack_corr in)
{
  int n=in.size();
  for(int t=0;t<n;t++)
    {
      out<<t<<" "<<in[t];
      if(t!=(n-1)) out<<endl;
    }
  
  return out;
}

jack_corr effective_mass(jack_corr in)
{
  jack_corr out;
  
  if(in.size()!=T/2+1) crash("corr must be half of the lattice");

  for(int t=0;t<T/2;t++)
    out.push_back(logj(in[t]/in[t+1]));
  
  return out;
}

int icombo(int r,int im1,int im2,int ith)
{return im1+nmass*(r+2*(im2+nmass*ith));}

jack_corr load_corr(const char *path,int im1,int im2,int ith)
{
  ifstream in(path);
  if(!in.good()) crash("File %s not found",path);
  
  int ic1=icombo(0,im1,im2,ith);
  int ic2=icombo(1,im1,im2,ith);
  
  int icorr=0;
  
  jack_corr c1,c2;
  while(icorr<=max(ic1,ic2) && in.good())
    {
      jack_corr ct;
      in>>ct;
      if(icorr==ic1) c1=ct;
      if(icorr==ic2) c2=ct;
      
      icorr++;
    };
  
  if(icorr!=max(ic1,ic2)+1) crash("corr not loaded");
  
  return (c1+c2)/2;
}

int main(int narg,char **arg)
{
  if(narg<5) crash("Use: %s path im1 im2 ith",arg[0]);
  int im1=atoi(arg[2]);
  int im2=atoi(arg[3]);
  int ith=atoi(arg[4]);
  
  cout.precision(16);
  
  cout<</*effective_mass*/(load_corr(arg[1],im1,im2,ith).symmetrized())<<endl;
  
  return 0;
}
