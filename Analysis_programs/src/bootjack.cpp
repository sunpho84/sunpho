#include <math.h>
#include <string.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>

#include "rand.cpp"

using namespace std;

class TYPE
{
public:
  int njack;
#ifdef BOOT
  int nboot;
#endif
  double *data;

#ifdef BOOT
  void create(int a,int b){nboot=a;njack=b;data=new double[a+1];for(int i=0;i<=a;i++) data[i]=0;}
  TYPE(const TYPE& in){create(in.nboot,in.njack);put(in.data);}
  explicit TYPE(){data=NULL;nboot=njack=0;}
  explicit TYPE(int a,int b){create(a,b);}
  explicit TYPE(int nb,int nj,int *in){double temp[nb+1];for(int i=0;i<nb+1;i++) temp[i]=in[i];TYPE(nb,nj,temp);}
  explicit TYPE(int nb,int nj,double *in){create(nb,nj);put(in);}
  void reallocate_if_necessary(int nb,int nj){if(nboot!=nb||njack!=nj){if(data!=NULL) delete[] data;create(nb,nj);}}
  TYPE operator=(const TYPE& in){reallocate_if_necessary(in.nboot,in.njack);put(in.data);return *this;}
#else
  void clusterize(int clust_size=1)
  {
    data[njack]=0;
    for(int ijack=0;ijack<njack;ijack++) data[njack]+=data[ijack];
    for(int ijack=0;ijack<njack;ijack++) data[ijack]=(data[njack]-data[ijack])/((njack-1)*clust_size);
    data[njack]/=njack*clust_size;
  }
  void create(int nj){njack=nj;data=new double[nj+1];for(int i=0;i<=nj;i++) data[i]=0;}
  TYPE(const TYPE& in){create(in.njack);put(in.data);}
  explicit TYPE(){data=NULL;njack=0;}
  explicit TYPE(int nj){create(nj);}
  explicit TYPE(int nj,int *in){double temp[nj+1];for(int i=0;i<nj+1;i++) temp[i]=in[i];TYPE(nj,temp);}
  explicit TYPE(int nj,double *in){create(nj);put(in);}
  void reallocate_if_necessary(int nj){if(njack!=nj){if(data!=NULL) delete[] data;create(nj);}}
  TYPE operator=(const TYPE& in){reallocate_if_necessary(in.njack);put(in.data);return *this;}
#endif
  ~TYPE(){if(data!=NULL) delete[]data;data=NULL;}
  
  double &operator[](int i){return data[i];}
  TYPE operator=(double a){for(int i=0;i<N+1;i++) data[i]=a;return *this;}
  
  double med(){return data[N];}
  double true_med(){double sx=0;for(int ij=0;ij<N;ij++) sx+=data[ij];return sx/N;}
  double err();
  
  void fill_gauss(double med,double sigma,int seed);
  void put(double* in){memcpy(data,in,sizeof(double)*(N+1));}
  void get(double* out){memcpy(out,data,sizeof(double)*(N+1));}
  
  TYPE write_to_file(string path);
  TYPE append_to_binfile(const char*,...);
  TYPE write_to_binfile(const char*,...);
  TYPE load(const char*,int);
};

void TYPE::fill_gauss(double med,double sigma,int seed)
{
  ran_gen extr(seed);
  for(int itype=0;itype<N;itype++) data[itype]=extr.get_gauss(med,sigma/sqrt(njack-1));
  data[N]=med;
}

#ifdef BOOT
boot fill_gauss(double med,double sigma,int seed,int nboots,int njacks)
{boot e(nboots,njacks);e.fill_gauss(med,sigma,seed);return e;}
#else
jack fill_gauss(double med,double sigma,int seed,int njacks)
{jack e(njacks);e.fill_gauss(med,sigma,seed);return e;}
#endif

double TYPE::err()
{
  if(njack<=1) return 0;
  
  double sx=0,s2x=0;
  
  for(int ij=0;ij<N;ij++)
    {
      sx+=data[ij];
      s2x+=data[ij]*data[ij];
    }
  sx/=N;
  s2x/=N;
  s2x-=sx*sx;
  
  return sqrt(fabs(s2x)*(njack-1));
}

TYPE TYPE::append_to_binfile(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  FILE *fout=open_file(buffer,"aw");
  int nw=fwrite(data,sizeof(double)*(N+1),1,fout);
  if(nw!=1)
    {
      cerr<<"Error appending to file "<<buffer<<endl;
      exit(1);
    }
  fclose(fout);
  
  return *this;
}

TYPE TYPE::write_to_binfile(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  FILE *fout=open_file(buffer,"w");
  int nw=fwrite(data,sizeof(double)*(N+1),1,fout);
  if(nw!=1)
    {
      cerr<<"Error appending to file "<<buffer<<endl;
      exit(1);
    }
  fclose(fout);
  
  return *this;
}

TYPE TYPE::load(const char *path,int i=0)
{
  FILE *fin=open_file(path,"r");
  if(fseeko(fin,(off_t)i*sizeof(double)*(N+1),SEEK_SET))
    {
      cerr<<"Error while searching for boot or jack "<<i<<" in file "<<path<<"!"<<endl;
      exit(1);
    }
  int nw=fread(data,sizeof(double),N+1,fin);
  if(nw!=N+1)
    {
      cerr<<"Error reading from file "<<path<<endl;
      exit(1);
    }
  fclose(fin);
  
  return *this;
}

ostream& operator<<(ostream &out,const TYPE obj)
{
  double med=TYPE(obj).med();
  double err=TYPE(obj).err();

  if(!std::isnan(med) && !std::isnan(err)) out<<med<<" "<<err;
  
  return out;
}

//functions

TYPE single_operator(const TYPE &a,double (*fun)(const double))
{
  int N=a.N;
  TYPE c(a);
  transform(a.data,a.data+N+1,c.data,ptr_fun(fun));

  return c;
}

TYPE pair_operator(const TYPE &a,const TYPE &b,double (*fun)(const double,const double))
{
  int N=a.N;
  if(b.N!=N)
    {
      cerr<<"Error, unmatched N: "<<N<<" "<<b.N<<"!"<<endl;
      exit(1);
    }
  TYPE c(a);
  transform(a.data,a.data+N+1,b.data,c.data,ptr_fun(fun));

  return c;
}

TYPE pair_operator(const TYPE &a,const double b,double (*fun)(const double,const double),int du)
{
  int N=a.N;
  TYPE c(a);
  
  for(int i=0;i<N+1;i++)
    if(du==2) c.data[i]=fun(a.data[i],b);
    else      c.data[i]=fun(b,a.data[i]);
  
  return c;
}

TYPE pair_operator(const TYPE &a,const double b,double (*fun)(const double,const double))
{return pair_operator(a,b,fun,2);}
TYPE pair_operator(const double a,const TYPE &b,double (*fun)(const double,const double))
{return pair_operator(b,a,fun,1);}

TYPE operator+(const TYPE &a,const TYPE &b){return pair_operator(a,b,double_summ);}
TYPE operator-(const TYPE &a,const TYPE &b){return pair_operator(a,b,double_subt);}
TYPE operator*(const TYPE &a,const TYPE &b){return pair_operator(a,b,double_prod);}
TYPE operator/(const TYPE &a,const TYPE &b){return pair_operator(a,b,double_frac);}

TYPE operator+=(TYPE &a,const TYPE &b){return a=a+b;}
TYPE operator-=(TYPE &a,const TYPE &b){return a=a-b;}
TYPE operator*=(TYPE &a,const TYPE &b){return a=a*b;}
TYPE operator/=(TYPE &a,const TYPE &b){return a=a/b;}

/////////////

TYPE operator+(const TYPE &a,const double b){return pair_operator(a,b,double_summ);}
TYPE operator-(const TYPE &a,const double b){return pair_operator(a,b,double_subt);}
TYPE operator*(const TYPE &a,const double b){return pair_operator(a,b,double_prod);}
TYPE operator/(const TYPE &a,const double b){return pair_operator(a,b,double_frac);}

TYPE operator+=( TYPE &a,const double b){return a=a+b;}
TYPE operator-=( TYPE &a,const double b){return a=a-b;}
TYPE operator*=( TYPE &a,const double b){return a=a*b;}
TYPE operator/=( TYPE &a,const double b){return a=a/b;}

TYPE operator+(const double a,const TYPE &b){return pair_operator(a,b,double_summ);}
TYPE operator-(const double a,const TYPE &b){return pair_operator(a,b,double_subt);}
TYPE operator*(const double a,const TYPE &b){return pair_operator(a,b,double_prod);}
TYPE operator/(const double a,const TYPE &b){return pair_operator(a,b,double_frac);}

TYPE operator+(const TYPE &a){return single_operator(a,unary_plus);}
TYPE operator-(const TYPE &a){return single_operator(a,unary_minus);}

TYPE sin(const TYPE &a){return single_operator(a,sin);}
TYPE cos(const TYPE &a){return single_operator(a,cos);}
TYPE tan(const TYPE &a){return single_operator(a,tan);}
TYPE asin(const TYPE &a){return single_operator(a,asin);}
TYPE acos(const TYPE &a){return single_operator(a,acos);}
TYPE atan(const TYPE &a){return single_operator(a,atan);}

TYPE sinh(const TYPE &a){return single_operator(a,sinh);}
TYPE cosh(const TYPE &a){return single_operator(a,cosh);}
TYPE tanh(const TYPE &a){return single_operator(a,tanh);}
TYPE asinh(const TYPE &a){return single_operator(a,asinh);}
TYPE acosh(const TYPE &a){return single_operator(a,acosh);}
TYPE atanh(const TYPE &a){return single_operator(a,atanh);}

TYPE fabs(const TYPE &a){return single_operator(a,fabs);}
TYPE exp(const TYPE &a){return single_operator(a,exp);}
TYPE log(const TYPE &a){return single_operator(a,log);}
TYPE sqr(const TYPE &a){return single_operator(a,sqr);}
TYPE sqrt(const TYPE &a){return single_operator(a,sqrt);}
TYPE pow(const TYPE &a,double b){return pair_operator(a,b,pow);}

#ifdef BOOT
TYPE boot_weighted_average(TYPE &a,TYPE &b)
#else
TYPE jack_weighted_average(TYPE &a,TYPE &b)
#endif
{
  double ea=a.err();
  double eb=b.err();

  double wa=1/(ea*ea);
  double wb=1/(eb*eb);

  return (a*wa+b*wb)/(wa+wb);
}

string smart_print(TYPE a)
{
  double m=a.med();
  double e=a.err();
  
  return smart_print(m,e);
}

string raw_print(TYPE a)
{
  ostringstream os;
  os.precision(16);
  for(int i=0;i<=a.N;i++) os<<a[i]<<" ";
  
  return os.str();
}


double med(TYPE x)
{
  double s=0;
  for(int i=0;i<x.N;i++) s+=x[i];
  return s/x.N;
}

double cov(TYPE x,TYPE y)
{return (med(x*y)-med(x)*med(y))*(x.njack-1);}

double var(TYPE x)
{return cov(x,x);}


TYPE TYPE::write_to_file(string path)
{
  ofstream out(path);
  out<<(*this)<<endl;
  
  return *this;
}
