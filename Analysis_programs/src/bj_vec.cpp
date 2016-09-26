#include <math.h>
#include <string.h>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#ifndef NOROOT
#include <TMatrixD.h>
#endif

using namespace std;

class VTYPE
{
public:
  int nel;
#ifdef BVEC
  int nboot;
#endif
  int njack;
  TYPE *data;
  
#ifdef BVEC
  void create(int ne,int nb,int nj){nel=ne;nboot=nb;njack=nj;data=new boot[nel];for(int iel=0;iel<nel;iel++) data[iel].create(nb,nj);}
  VTYPE(){data=NULL;nel=nboot=njack=0;}
  VTYPE(const VTYPE& a) : nel(a.nel),nboot(a.nboot),njack(a.njack){create(nel,nboot,njack);for(int iel=0;iel<nel;iel++) data[iel]=a.data[iel];}
  VTYPE(int ne,int nb,int nj){create(ne,nb,nj);}
  void reallocate_if_necessary(int ne,int nb,int nj){if(nel!=ne||njack!=nj||nb!=nboot){delete[] data;create(ne,nb,nj);}}
  VTYPE operator=(const VTYPE& a){
    reallocate_if_necessary(a.nel,a.nboot,a.njack);for(int iel=0;iel<nel;iel++) data[iel]=a.data[iel];return *this;}
#else
  void clusterize(int clust_size=1){for(int iel=0;iel<nel;iel++) data[iel].clusterize(clust_size);}
  void create(int ne,int nj){nel=ne;njack=nj;data=new jack[nel];for(int iel=0;iel<nel;iel++) data[iel].create(nj);}
  explicit VTYPE(){data=NULL;nel=njack=0;}
  VTYPE(const VTYPE& a) : nel(a.nel),njack(a.njack){create(nel,njack);for(int iel=0;iel<nel;iel++) data[iel]=a.data[iel];}
  VTYPE(int ne,int nj){create(ne,nj);}
  void reallocate_if_necessary(int ne,int nj){if(nel!=ne||njack!=nj){if(data!=NULL)delete[] data;create(ne,nj);}}
  VTYPE operator=(const VTYPE& a){
    reallocate_if_necessary(a.nel,a.njack);for(int iel=0;iel<nel;iel++) data[iel]=a.data[iel];return *this;}
#endif
  ~VTYPE(){if(data!=NULL) delete[] data;data=NULL;}
  
  void put(double **out){for(int iel=0;iel<nel;iel++) data[iel].put(out[iel]);}
  void put(double *out){for(int iel=0;iel<nel;iel++) data[iel].put(out+iel*(N+1));}
  void get(double *in){for(int iel=0;iel<nel;iel++) data[iel].get(in+iel*(N+1));}
  VTYPE load(string,int);
  VTYPE load(FILE *,int);
  VTYPE load(FILE *);
  VTYPE load_naz(string,int);
  void print_to_file(string);
  void print_to_file(ofstream &);
  void print_rel_err_to_file(string);
  
  TYPE& operator[](int i){return data[i];}
  VTYPE operator=(double in){for(int iel=0;iel<nel;iel++) data[iel]=in;return *this;}
  VTYPE operator=(const TYPE& a){
    if(data==NULL){cerr<<"Error, using unallocated vector!"<<endl;exit(1);}for(int iel=0;iel<nel;iel++) data[iel]=a;return *this;}

  
  VTYPE subset(int,int);
  VTYPE first_half();
  VTYPE simmetric();
  VTYPE inverted();
  VTYPE simmetrized(int parity);
  VTYPE shifted(int);
  VTYPE write_to_binfile(FILE *);
  VTYPE append_to_binfile(const char *format,...);
  VTYPE write_to_binfile(const char *format,...);
};

ostream& operator<<(ostream &out,const VTYPE &obj);

VTYPE VTYPE::write_to_binfile(FILE *fout)
{
  double out[nel*(N+1)];
  get(out);
  int nw=fwrite(out,sizeof(double),(N+1)*nel,fout);
  if(nw!=(N+1)*nel)
    {
      cerr<<"Error writing to file"<<endl;
      exit(1);
    }
  
  return *this;
}

VTYPE VTYPE::append_to_binfile(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  FILE *fout=open_file(buffer,"aw");
  write_to_binfile(fout);
  fclose(fout);
  
  return *this;
}

VTYPE VTYPE::write_to_binfile(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  FILE *fout=open_file(buffer,"w");
  write_to_binfile(fout);
  fclose(fout);
  
  return *this;
}

VTYPE VTYPE::load(FILE *fin)
{
  double *in=new double[nel*(N+1)];
  
  int stat=fread(in,sizeof(double),nel*(N+1),fin);
  if(stat!=nel*(N+1))
    {
      if(stat==EOF)
	{
	  fprintf(stderr,"Error, reached EOF while reading data!\n");
	  exit(1);
	}
      else
        {
	  fprintf(stderr,"Error while reading data from file, obtained: %d while reading %d elements\n",stat,nel*(N+1));
	  exit(1);
	}
    }
  
  put(in);
  
  delete[] in;
  
  return *this;
}  

VTYPE VTYPE::load(FILE *fin,int i)
{
  if(fseeko(fin,(off_t)i*sizeof(double)*nel*(N+1),SEEK_SET))
    {
      fprintf(stderr,"Error while searching for correlation %d!\n",i);
      exit(1);
    }
  
  load(fin);
  
  return *this;
}  

VTYPE VTYPE::load(string path,int i)
{
  if(debug_load) cout<<"Loading corr "<<i<<" from path "<<path<<endl;
  FILE *fin=open_file(path,"r");
  
  load(fin,i);
  fclose(fin);
  
  return (*this);
}

void VTYPE::print_to_file(ofstream &fout)
{
  int i=fout.precision();
  fout.precision(16);
  if(!(fout<<"@type xydy"<<endl)) crash("writing header");
  fout<<(*this);
  fout.precision(i);
}

void VTYPE::print_to_file(string path)
{
  ofstream fout(path);
  if(!fout.good()) crash("opening %s [print_to_file]",path.c_str());
  
  print_to_file(fout);
  
  fout.close();
}

void VTYPE::print_rel_err_to_file(string path)
{
  ofstream fout(path);
  if(!fout.good()) crash("opening %s [print_rel_err_to_file]",path.c_str());
  
  fout<<"@type xy"<<endl;
  for(int t=0;t<this->nel;t++) fout<<(*this)[t].med()/(*this)[t].err()<<endl;
  fout.close();
}

VTYPE VTYPE::load_naz(string path,int icorr)
{
  double in[nel*2*(N+1)];
  
  FILE *fin=open_file(path,"r");
  
  int off=(off_t)2*(icorr/2)*sizeof(double)*nel*(N+1);
  if(fseeko(fin,off,SEEK_SET))
    {
      fprintf(stderr,"Error while searching for correlation %d!\n",icorr);
      exit(1);
    }
  int stat=fread(in,sizeof(double)*2*nel*(N+1),1,fin);
  if(stat!=1)
    {
      if(stat==EOF)
	{
	  fprintf(stderr,"Error, reached EOF while reading data!\n");
	  exit(1);
	}
      else
        {
	  perror("Error while reading data");
	  exit(1);
	}
    }
  
  int ri=icorr%2;
  for(int iel=0;iel<nel;iel++)
    for(int ijack=0;ijack<N+1;ijack++)
      data[iel].data[ijack]=in[ijack+(N+1)*(ri+2*iel)];
  
  fclose(fin);
  
  return (*this);
}

#ifdef BVEC
bvec bvec_load(string path,int nel,int nboot,int njack,int i)
{
  bvec out(nel,nboot,njack);
  
  return out.load(path,i);
}

bvec bvec_load_naz(string path,int nel,int nboot,int njack,int i)
{
  bvec out(nel,nboot,njack);
  
  return out.load_naz(path,i);
}
#else
jvec jvec_load(string path,int nel,int njack,int i)
{
  jvec out(nel,njack);
  
  return out.load(path,i);
}

jvec jvec_load_naz(string path,int nel,int njack,int i)
{
  jvec out(nel,njack);
  
  return out.load_naz(path,i);
}
#endif

ostream& operator<<(ostream &out,const VTYPE &obj)
{
  for(int iel=0;iel<obj.nel;iel++)
    {
      double med=obj.data[iel].med();
      double err=obj.data[iel].err();
      if(!std::isnan(med) && !std::isnan(err))
	out<<iel<<" "<<obj.data[iel]<<endl;
    }
  
  return out;
}

VTYPE single_operator(const VTYPE &a,double (*fun)(const double))
{
  int nel=a.nel;
  int njack=a.njack;
#ifdef BVEC
  int nboot=a.nboot;
  VTYPE c(nel,nboot,njack);
#else
  VTYPE c(nel,njack);
#endif

  for(int iel=0;iel<nel;iel++) c.data[iel]=single_operator(a.data[iel],fun);

  return c;
}

VTYPE pair_operator(const VTYPE &a,const VTYPE &b,double (*fun)(const double,const double))
{
  int nel=a.nel;
  int njack=a.njack;
#ifdef BVEC
  int nboot=a.nboot;
  if(b.nboot!=nboot){cerr<<"Error, unmatched nboot!"<<endl;exit(1);}
  if(b.njack!=njack){cerr<<"Error, unmatched njack!"<<endl;exit(1);}
  if(b.nel!=nel){cerr<<"Error, unmatched nel: "<<nel<<" vs "<<b.nel<<"!"<<endl;exit(1);}
  bvec c(nel,nboot,njack);
#else
  if(b.njack!=njack){cerr<<"Error, unmatched njack!"<<endl;exit(1);}
  if(b.nel!=nel){cerr<<"Error, unmatched nel: "<<nel<<" vs "<<b.nel<<"!"<<endl;exit(1);}
  jvec c(nel,njack);
#endif
  
  for(int iel=0;iel<nel;iel++) c.data[iel]=pair_operator(a.data[iel],b.data[iel],fun);
  
  return c;
}

VTYPE pair_operator(const VTYPE &a,const double b,double (*fun)(const double,const double),int du)
{
  int nel=a.nel;
  int njack=a.njack;
#ifdef BVEC
  int nboot=a.nboot;
  bvec c(nel,nboot,njack);
#else
  jvec c(nel,njack);
#endif
  
  for(int iel=0;iel<nel;iel++) c.data[iel]=pair_operator(a.data[iel],b,fun,du);
  
  return c;
}

VTYPE pair_operator(const VTYPE &a,const double b,double (*fun)(const double,const double))
{return pair_operator(a,b,fun,2);}
VTYPE pair_operator(const double a,const VTYPE &b,double (*fun)(const double,const double))
{return pair_operator(b,a,fun,1);}

VTYPE pair_operator(const VTYPE &a,const TYPE b,double (*fun)(const double,const double),int du)
{
  int nel=a.nel;
#ifdef BVEC
  bvec c(nel,a.nboot,a.njack);
#else
  jvec c(nel,a.njack);
#endif
  
  for(int iel=0;iel<nel;iel++)
    if(du==2) c.data[iel]=pair_operator(a.data[iel],b,fun);
    else      c.data[iel]=pair_operator(b,a.data[iel],fun);
  
  return c;
}

VTYPE pair_operator(const VTYPE &a,const TYPE b,double (*fun)(const double,const double))
{return pair_operator(a,b,fun,2);}
VTYPE pair_operator(const TYPE a,const VTYPE &b,double (*fun)(const double,const double))
{return pair_operator(b,a,fun,1);}

VTYPE operator+(const VTYPE &a,const VTYPE &b){return pair_operator(a,b,double_summ);}
VTYPE operator-(const VTYPE &a,const VTYPE &b){return pair_operator(a,b,double_subt);}
VTYPE operator*(const VTYPE &a,const VTYPE &b){return pair_operator(a,b,double_prod);}
VTYPE operator/(const VTYPE &a,const VTYPE &b){return pair_operator(a,b,double_frac);}

VTYPE operator+=(VTYPE &a,const VTYPE &b){return a=a+b;}
VTYPE operator-=(VTYPE &a,const VTYPE &b){return a=a-b;}
VTYPE operator*=(VTYPE &a,const VTYPE &b){return a=a*b;}
VTYPE operator/=(VTYPE &a,const VTYPE &b){return a=a/b;}

///////////

VTYPE operator+(const VTYPE &a,const TYPE &b){return pair_operator(a,b,double_summ);}
VTYPE operator-(const VTYPE &a,const TYPE &b){return pair_operator(a,b,double_subt);}
VTYPE operator*(const VTYPE &a,const TYPE &b){return pair_operator(a,b,double_prod);}
VTYPE operator/(const VTYPE &a,const TYPE &b){return pair_operator(a,b,double_frac);}

VTYPE operator+(const TYPE &a,const VTYPE &b){return pair_operator(a,b,double_summ);}
VTYPE operator-(const TYPE &a,const VTYPE &b){return pair_operator(a,b,double_subt);}
VTYPE operator*(const TYPE &a,const VTYPE &b){return pair_operator(a,b,double_prod);}
VTYPE operator/(const TYPE &a,const VTYPE &b){return pair_operator(a,b,double_frac);}

VTYPE operator+=(VTYPE &a,const TYPE &b){return a=a+b;}
VTYPE operator-=(VTYPE &a,const TYPE &b){return a=a-b;}
VTYPE operator*=(VTYPE &a,const TYPE &b){return a=a*b;}
VTYPE operator/=(VTYPE &a,const TYPE &b){return a=a/b;}

////////////

VTYPE operator+(const VTYPE &a,const double b){return pair_operator(a,b,double_summ);}
VTYPE operator-(const VTYPE &a,const double b){return pair_operator(a,b,double_subt);}
VTYPE operator*(const VTYPE &a,const double b){return pair_operator(a,b,double_prod);}
VTYPE operator/(const VTYPE &a,const double b){return pair_operator(a,b,double_frac);}

VTYPE operator+=(VTYPE &a,const double b){return a=a+b;}
VTYPE operator-=(VTYPE &a,const double b){return a=a-b;}
VTYPE operator*=(VTYPE &a,const double b){return a=a*b;}
VTYPE operator/=(VTYPE &a,const double b){return a=a/b;}

VTYPE operator+(const double a,const VTYPE &b){return pair_operator(a,b,double_summ);}
VTYPE operator-(const double a,const VTYPE &b){return pair_operator(a,b,double_subt);}
VTYPE operator*(const double a,const VTYPE &b){return pair_operator(a,b,double_prod);}
VTYPE operator/(const double a,const VTYPE &b){return pair_operator(a,b,double_frac);}

////////////

VTYPE operator+(const VTYPE &a){return single_operator(a,unary_plus);}
VTYPE operator-(const VTYPE &a){return single_operator(a,unary_minus);}

VTYPE sin(const VTYPE &a){return single_operator(a,sin);}
VTYPE cos(const VTYPE &a){return single_operator(a,cos);}
VTYPE tan(const VTYPE &a){return single_operator(a,tan);}
VTYPE asin(const VTYPE &a){return single_operator(a,asin);}
VTYPE acos(const VTYPE &a){return single_operator(a,acos);}
VTYPE atan(const VTYPE &a){return single_operator(a,atan);}

VTYPE sinh(const VTYPE &a){return single_operator(a,sinh);}
VTYPE cosh(const VTYPE &a){return single_operator(a,cosh);}
VTYPE tanh(const VTYPE &a){return single_operator(a,tanh);}
VTYPE asinh(const VTYPE &a){return single_operator(a,asinh);}
VTYPE acosh(const VTYPE &a){return single_operator(a,acosh);}
VTYPE atanh(const VTYPE &a){return single_operator(a,atanh);}

VTYPE fabs(const VTYPE &a){return single_operator(a,fabs);}
VTYPE sqr(const VTYPE &a){return single_operator(a,sqr);}
VTYPE log(const VTYPE &a){return single_operator(a,log);}
VTYPE sqrt(const VTYPE &a){return single_operator(a,sqrt);}
VTYPE pow(const VTYPE &a,double b){return pair_operator(a,b,pow);}

//////////////

VTYPE VTYPE::subset(int tmin,int tmax)
{
  if(tmin<0||tmin>nel||tmax<0||tmax>nel)
    {
      cerr<<"Error, required impossible interval ["<<tmin<<","<<tmax<<"] of a "<<nel<<"long vector!"<<endl;
      exit(1);
    }
  int diff=tmax-tmin;
  
#ifdef BVEC
  bvec c(diff,nboot,njack);
#else
  jvec c(diff,njack);
#endif
  for(int iel=0;iel<diff;iel++)
    c.data[iel]=data[iel+tmin];
  
  return c;
}

VTYPE VTYPE::first_half()
{
  if(nel%2!=0)
    {
      cerr<<"Error, required the first half of an odd-length ("<<nel<<") VTYPEtor!"<<endl;
      exit(1);
    }
  
  return subset(0,nel/2+1);
}

VTYPE VTYPE::simmetric()
{
#ifdef BVEC
  bvec c(nel,nboot,njack);
#else
  jvec c(nel,njack);
#endif
  
  for(int iel=0;iel<nel;iel++)
    c.data[iel]=data[(nel-iel)%nel];

  return c;
}

VTYPE VTYPE::inverted()
{
#ifdef BVEC
  bvec c(nel,nboot,njack);
#else
  jvec c(nel,njack);
#endif
  
  for(int iel=0;iel<nel;iel++)
    c.data[iel]=data[nel-iel-1];

  return c;
}

VTYPE VTYPE::simmetrized(int parity)
{
  if(abs(parity)!=1&&parity)
    {
      cerr<<"Error, parity required for simmetrization: "<<parity<<endl;
      exit(1);
    }
  
  if(nel%2!=0)
    {
      cerr<<"Error, required to simmetrize an odd-length VTYPEtor!"<<endl;
      exit(1);
    }
  
#ifdef BVEC
  bvec c(nel/2+1,nboot,njack);
#else
  jvec c(nel/2+1,njack);
#endif
  
  for(int iel=0;iel<nel/2+1;iel++)
    c.data[iel]=(data[iel]+parity*data[(nel-iel)%nel])/(1+abs(parity));
  
  return c;
}

VTYPE VTYPE::shifted(int am)
{
  VTYPE a(*this);
  
  for(int iel=0;iel<nel;iel++)
    {
      int iso=(iel-am+nel)%nel;
      a[iel]=(*this)[iso];
    }
  
  return a;
}  

VTYPE rectangle_integrate(VTYPE in)
{
#ifdef BVEC
  bvec out(in.nel+1,in.nboot,in.njack);
#else
  jvec out(in.nel+1,in.njack);
#endif
  
  out[0]=0;
  for(int iel=1;iel<=in.nel;iel++)
    out[iel]=out[iel-1]+in[iel-1];
  
  return out;
}  

VTYPE paste(VTYPE a,VTYPE b)
{
#ifdef BVEC
  bvec c(a.nel+b.nel,a.nboot,a.njack);
#else
  jvec c(a.nel+b.nel,a.njack);
#endif
  
  for(int iel=0;iel<a.nel;iel++)
    c[iel]=a[iel];
  for(int iel=0;iel<b.nel;iel++)
    c[a.nel+iel]=b[iel];
  
  return c;
}  

VTYPE interleave(VTYPE a,VTYPE b)
{
  if(a.nel!=b.nel) crash("impossible to interleave %d and %d",a.nel,b.nel);
#ifdef BVEC
  bvec c(a.nel+b.nel,a.nboot,a.njack);
#else
  jvec c(a.nel+b.nel,a.njack);
#endif
  
  for(int iel=0;iel<a.nel;iel++)
    {
      c[2*iel+0]=a[iel];
      c[2*iel+1]=b[iel];
    }
  
  return c;
}  

string write_constant_fit_plot(VTYPE in,TYPE y,int tin,int tfin,int iset=0)
{
  ostringstream out;
  double ym=y.med(),dy=y.err();
  //error of the line
  out<<"@s"<<iset+0<<" line type 1"<<endl;
  out<<"@s"<<iset+0<<" line color "<<iset/3%7+8<<endl;
  out<<"@s"<<iset+0<<" fill color "<<iset/3%7+8<<endl;
  out<<"@s"<<iset+0<<" fill type 1"<<endl;
  out<<"@type xy"<<endl;
  out<<tin-0.4<<" "<<ym-dy<<endl<<tfin+0.4<<" "<<ym-dy<<endl;
  out<<tfin+0.4<<" "<<ym+dy<<endl<<tin-0.4<<" "<<ym+dy<<endl;
  out<<tin-0.4<<" "<<ym-dy<<endl;
  out<<"&"<<endl;
  //central line
  out<<"@s"<<iset+1<<" line color "<<iset/3%7+1<<endl;
  out<<"@type xy"<<endl;
  out<<tin-0.4<<" "<<ym<<endl<<tfin+0.4<<" "<<ym<<endl;
  //plot the original data with error
  out<<"&"<<endl;
  out<<"@type xydy"<<endl;
  out<<"@s"<<iset+2<<" line type 0"<<endl;
  out<<"@s"<<iset+2<<" symbol color "<<iset/3%7+1<<endl;
  out<<"@s"<<iset+2<<" errorbar color "<<iset/3%7+1<<endl;
  out<<"@s"<<iset+2<<" symbol "<<iset/3%7+1<<endl;
  out<<in;
  out<<"&"<<endl;
  
  return out.str();
}
void append_constant_fit_plot(string path,VTYPE in,TYPE y,int tin,int tfin,int iset)
{
  ofstream out(path,ios::app);
  if(!out.good()) crash("opening %s [append_constant_fit_plot]",path.c_str());
  out<<write_constant_fit_plot(in,y,tin,tfin,iset);
  out.close();
}
void write_constant_fit_plot(string path,VTYPE in,TYPE y,int tin,int tfin,int iset=0)
{
  ofstream out(path);
  if(!out.good()) crash("opening %s [write_constant_fit_plot]",path.c_str());
  out<<"@page size 800,600"<<endl;
  out.close();
  append_constant_fit_plot(path,in,y,tin,tfin,iset);
}

TYPE constant_fit(VTYPE in,int tin,int tfin,string path="",string path_ch2="")
{
  TYPE E(in.data[0]);

  E=0;
  double norm=0;
  
  //take weighted average
  for(int iel=max(tin,0);iel<=min(tfin,in.nel-1);iel++)
    {
      TYPE ele=in.data[iel];
      double err=in.data[iel].err();
      double weight=1/(err*err);
      if(!std::isnan(err)&&err!=0)
	{
	  E+=ele*weight;
	  norm+=weight;
	}
    }
  
  //take simply average in the other case
  if(norm==0)
    for(int iel=max(tin,0);iel<=min(tfin,in.nel-1);iel++)
      {
	norm=norm+1;
	E+=in.data[iel];
      }
  
  //normalize
  E/=norm;
  
  if(path!="") write_constant_fit_plot(path,in,E,tin,tfin);
  if(path_ch2!="")
    {
      ofstream out(path_ch2);
      out<<"t ((teo-data)/err)^2=chi2_contr"<<endl;
      out<<"================================="<<endl;
      double tot=0;
      int ndof=0;
      for(int iel=max(tin,0);iel<=min(tfin,in.nel-1);iel++)
	{
	  double contr=sqr((in.data[iel].med()-E.med())/in.data[iel].err());
	  tot+=contr;
	  ndof++;
	  out<<iel<<" (("<<in.data[iel].med()<<"-"<<E.med()<<")/"<<in.data[iel].err()<<")^2="<<contr<<endl;
	}
      out<<"================================="<<endl;
      out<<"Total chi2: "<<tot<<"/"<<ndof-1<<"="<<tot/(ndof-1)<<endl;
    }
  
  return E;
}

#ifndef NOROOT

TYPE correlated_constant_fit(VTYPE in,int tin,int tfin,string path="")
{
  tin=max(tin,0);
  tfin=min(tfin,in.nel);
  
  //build covariance matrix
  int nel=tfin-tin;
  TMatrixD cova(nel,nel);
  for(int i=0;i<nel;i++)
    for(int j=0;j<nel;j++)
      cova(i,j)=cov(in[i+tin],in[j+tin]);
  
  //take the inverse
  TMatrixD inv_cova=cova.Invert();
  
  TYPE E(in.data[0]),norm(in.data[0]);
  norm=E=0;
  
  //take weighted average
  for(int i=0;i<nel;i++)
    for(int j=0;j<nel;j++)
      {
	TYPE y=(in.data[i+tin]+in.data[j+tin])/2;
	
	E+=y*inv_cova(i,j);
	norm+=inv_cova(i,j);
      }
  
  //normalize
  E/=norm;
  
  if(path!="") write_constant_fit_plot(path,in,E,tin,tfin);
  
  return E;
}

#endif

string write_line_with_error(TYPE q,TYPE m,double xmin,double xmax,int npoints);
void linear_fit(TYPE &q,TYPE &m,double *x,VTYPE y,double xmin,double xmax,string plot_path="")
{
  double S,Sx,Sx2;
  TYPE Sxy(y[0]),Sy(y[0]);
  
  Sx2=S=Sx=0;
  Sxy=Sy=0;
  for(int iel=0;iel<y.nel;iel++)
    if(x[iel]>=xmin && x[iel]<=xmax)
      {
	double err=y[iel].err();
	double weight=1/(err*err);
	double xi=x[iel];
	TYPE yi=y[iel];
	
	S+=weight;
	Sx+=xi*weight;
	Sx2+=xi*xi*weight;
	Sxy+=xi*yi*weight;
	Sy+=yi*weight;
      }
  
  double delta=S*Sx2-Sx*Sx;
  m=(S*Sxy-Sx*Sy)/delta;
  q=(Sx2*Sy-Sxy*Sx)/delta;

  if(plot_path!="")
    {
      ofstream out(plot_path);
      out<<"@type xy"<<endl<<
	write_line_with_error(q,m,xmin-0.2,xmax+0.2,100)<<endl<<
	"&"<<endl<<
	"@type xydy"<<endl;
      
      for(int i=0;i<y.nel;i++)
	out<<x[i]<<" "<<y[i]<<endl;
    }
}

void linear_fit(TYPE &q,TYPE &m,VTYPE y,int tmin,int tmax,string plot_path="")
{
  double x[y.nel];
  for(int iel=0;iel<y.nel;iel++) x[iel]=iel;
  double xmin=max(0,tmin);
  double xmax=min(tmax+1,y.nel);
  
  linear_fit(q,m,x,y,xmin,xmax,plot_path);
}

VTYPE par_single_fun(double (*fun)(double,double*),VTYPE &x,VTYPE &par)
{
  int nx=x.nel;
  int npar=par.nel;
  
  VTYPE y(x);
  
  for(int ijack=0;ijack<x.njack+1;ijack++)
    {
      double dpar[npar];
      //create temporary double vector for pars
      for(int ipar=0;ipar<npar;ipar++) dpar[ipar]=par.data[ipar].data[ijack];
      for(int ix=0;ix<nx;ix++) y.data[ix].data[ijack]=fun(x.data[ix].data[ijack],dpar);
    }
  
  return y;
}
