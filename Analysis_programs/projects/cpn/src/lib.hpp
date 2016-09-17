#include <algorithm>
#include <cstdlib>
#include <cstdarg>
#include <cmath>
#include <fftw3.h>
#include <stdarg.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <vector>

using namespace std;
using namespace placeholders;

const int npart=4;

int part_size;
int nind_confs;
int *boot_id;

enum start_cond_t{COLD,HOT,LOAD};

int seed;
start_cond_t start_cond;
int nterm,nsweep,ncorr,nmicro,use_hmc;
int N;
double beta;
double g;
int L;
int nstout_lev;
double stout_rho;
int nhmc_steps;
int use_topo_pot;
int use_charge_pot;
double charge_pot,th_top;

#define crash(...) internal_crash(__LINE__,__VA_ARGS__)

vector<double> mag[2];
vector<double> ene;
vector<double> *corr;

int compute_corr_each;

//=========================== tools =================================

template <class T> T sqr(T in)
{return in*in;}

template<class T> T det3(T *a,T *b,T *c)
{
  T d;
  d= a[0]*(b[1]*c[2]-b[2]*c[1]);
  d+=b[0]*(c[1]*a[2]-c[2]*a[1]);
  d+=c[0]*(a[1]*b[2]-a[2]*b[1]);
  
  return d;
}

//crash promptin error message
void internal_crash(int line,const char *templ,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,templ);
  vsprintf(buffer,templ,args);
  va_end(args);
  
  fprintf(stderr,"ERROR on line %d: %s\n",line,buffer);
  exit(1);
}

//read an element from input file
template <class T> void read(T &out,ifstream &in,string is)
{
  string s;
  if(!(in>>s)) crash("impossible to read expect string \"%s\"",is.c_str());
  if(s!=is) crash("obtained %s while reading %s",s.c_str(),is.c_str());
  if(!(in>>out)) crash("reading data");
}

//open file
FILE *open_file(const char *path,const char *tag)
{
  FILE *file=fopen(path,tag);
  if(file==NULL) crash("problem opening file %s with %s",path,tag);
  return file;
}

//usual combination
string combine(const char *format,...)
{
  char buffer[1024];
  va_list args;
  
  va_start(args,format);
  vsprintf(buffer,format,args);
  va_end(args);
  
  return string(buffer);
}

//============================== types ============================

//returned as ave_err
typedef pair<double,double> ave_err_t;
ostream& operator<<(ostream &out,ave_err_t obj)
{return out<<obj.first<<" "<<obj.second;}
string smart_print(ave_err_t in)
{
  double m=in.first,e=in.second;
  double orim=m;
  int s=0;
  ostringstream o;
  
  if(e==0) o<<m<<"(0)";
  else
    if(e<1)
      {
        //print 1 digit of error or 2 if err starts with 1 or 2
        while(int(e)<=3)
          {
            e*=10;
            m*=10;
            s--;
          }
        
        char t[1000];
        sprintf(t,"%.*f",(unsigned int)abs(s),orim);
        
        o<<t<<"("<<int(e+0.5)<<")";
      }
    else
      {
        if(e>=3)
          {
            //count the numbr of digits to truncate
            int s=0;
	    o<<(int)(m+0.5)*pow(10,s)<<"("<<(int)(e+0.5)*pow(10,s)<<")";
	  }
	else
	  {
	    char t[1000];
	    sprintf(t,"%.1f(%.1f)",m,e);
	    o<<t;
	  }
      }
  
  return o.str();
}

void adjust_offset_slices(vector<double> *slices,int npart)
{
  int n=slices[0].size();
  
  //reset
  for(int islice=0;islice<npart;islice++)
    {
      double a=0;
      for(int i=0;i<n;i++) a+=slices[islice][i];
      a/=n;
      for(int i=0;i<n;i++) slices[islice][i]-=a;
    }
}

ave_err_t ave_err_part(function<double(int)> fun,int n=npart)
{
  double ave=0,err=0;
  for(int i=0;i<n;i++)
    {
      double a=fun(i);
      ave+=a;
      err+=a*a;
    }
  
  ave/=n;
  err/=n;
  err-=ave*ave;
  err=sqrt(err/(npart-1));
  
  return make_pair(ave,err);
}

//compute the average and the error
void reconstruct_band(vector<double> &ave,vector<double> &err,vector<double> *slices,int npart)
{
  int n=ave.size();
  
  for(int i=0;i<n;i++)
    {
      auto tmp=ave_err_part([slices,i,n](int islice)->double{return (slices[islice][i]+slices[islice][n-i-1])/2;});
      ave[i]=tmp.first;
      err[i]=tmp.second;
    }
}

int clust_size=100000;
int nboot=100;
struct boot_t : vector<double>
{
  void init(){resize(nboot);}
  boot_t(){init();}
  boot_t(double in){init();(*this)=in;}
  void populate(function<double(int)>);
  boot_t(function<double(int)> fun){init();populate(fun);}
  ave_err_t ave_err();
  double err() {return ave_err().second;}
  
  boot_t operator+(boot_t &oth)
  {boot_t out;transform(this->begin(),this->end(),oth.begin(),out.begin(),plus<double>());return out;}
  boot_t operator-(boot_t &oth)
  {boot_t out;transform(this->begin(),this->end(),oth.begin(),out.begin(),minus<double>());return out;}
  boot_t operator-()
  {boot_t out;transform(this->begin(),this->end(),out.begin(),bind(minus<double>(),0,_1));return out;}
  boot_t operator*(boot_t &oth)
  {boot_t out;transform(this->begin(),this->end(),oth.begin(),out.begin(),multiplies<double>());return out;}
  boot_t operator/(boot_t &oth)
  {boot_t out;transform(this->begin(),this->end(),oth.begin(),out.begin(),divides<double>());return out;}
  boot_t operator+(double oth)
  {boot_t out;transform(this->begin(),this->end(),out.begin(),bind(plus<double>(),_1,oth));return out;}
  boot_t operator-(double oth)
  {boot_t out;transform(this->begin(),this->end(),out.begin(),bind(minus<double>(),_1,oth));return out;}
  boot_t operator*(double oth)
  {boot_t out;transform(this->begin(),this->end(),out.begin(),bind(multiplies<double>(),_1,oth));return out;}
  boot_t operator/(double oth)
  {boot_t out;transform(this->begin(),this->end(),out.begin(),bind(divides<double>(),_1,oth));return out;}
  boot_t operator=(double in){for(auto &it : (*this)) it=in;return (*this);}
  boot_t operator+(boot_t &&oth){return (*this)+oth;}
  boot_t operator-(boot_t &&oth){return (*this)-oth;}
  boot_t operator*(boot_t &&oth){return (*this)*oth;}
  boot_t operator/(boot_t &&oth){return (*this)/oth;}
  boot_t operator+=(boot_t &oth)
  {transform(this->begin(),this->end(),oth.begin(),this->begin(),plus<double>());return (*this);}
  boot_t operator-=(boot_t &oth)
  {transform(this->begin(),this->end(),oth.begin(),this->begin(),minus<double>());return (*this);}
  boot_t operator*=(boot_t &oth)
  {transform(this->begin(),this->end(),oth.begin(),this->begin(),multiplies<double>());return (*this);}
  boot_t operator/=(boot_t &oth)
  {transform(this->begin(),this->end(),oth.begin(),this->begin(),divides<double>());return (*this);}
  boot_t operator+=(boot_t &&oth){return (*this)+=oth;}
  boot_t operator-=(boot_t &&oth){return (*this)-=oth;}
  boot_t operator*=(boot_t &&oth){return (*this)*=oth;}
  boot_t operator/=(boot_t &&oth){return (*this)/=oth;}
};

inline boot_t operator+(double a,boot_t &b) {return b+a;}
inline boot_t operator+(double a,boot_t &&b) {return a+b;}
inline boot_t operator-(double a,boot_t &b) {boot_t out;transform(b.begin(),b.end(),out.begin(),bind(minus<double>(),a,_1));return out;}
inline boot_t operator-(double a,boot_t &&b) {return a-b;}
inline boot_t operator*(double a,boot_t &b) {return b*a;}
inline boot_t operator*(double a,boot_t &&b) {return a*b;}
inline boot_t operator/(double a,boot_t &b) {boot_t out;transform(b.begin(),b.end(),out.begin(),bind(divides<double>(),a,_1));return out;}
inline boot_t operator/(double a,boot_t &&b) {return a/b;}
boot_t sqrt(boot_t &in){boot_t out;transform(in.begin(),in.end(),out.begin(),(double (*)(double))sqrt);return out;}
boot_t pow(boot_t &in,double a){boot_t out;transform(in.begin(),in.end(),out.begin(),bind((double (*)(double,double))pow,_1,a));return out;}
boot_t log(boot_t &in){boot_t out;transform(in.begin(),in.end(),out.begin(),(double (*)(double))log);return out;}
boot_t sqrt(boot_t &&in){return sqrt(in);}
boot_t pow(boot_t &&in,double b){return pow(in,b);}
boot_t log(boot_t &&in){return log(in);}
ostream& operator<<(ostream &out,boot_t &obj)
{for(auto &it : obj) out<<it<<endl;return out;}

class meta_pars_t
{
public:
  int after;
  int each;
  double coeff;
  double width;
  double barr;
  double force_out;
  double well_tempering;
  double bend;
  
  int ngrid;
  vector<double> weight,weight_slice[npart];
  vector<double> coll;
  vector<double> pote_slice[npart];
  vector<double> pote_ave,pote_err;
  
  boot_t conf_weight;
  boot_t conf_weight_slice[npart];
  boot_t corr_weight;
  
  meta_pars_t() : after(0),each(0),coeff(0.0),width(0.0),barr(0.0),force_out(0.0),well_tempering(0.0),bend(0.0) {}
  
  //read all the parameters
  void read_pars(ifstream &input,const char *tag)
  {
    read(after,input,combine("Chrono%sAfter",tag));
    read(each,input,combine("Chrono%sEach",tag));
    read(coeff,input,combine("Chrono%sCoeff",tag));
    read(width,input,combine("Chrono%sWidth",tag));
    read(barr,input,combine("Chrono%sBarr",tag));
    read(force_out,input,combine("Chrono%sForceOut",tag));
    read(bend,input,combine("Chrono%sBend",tag));
    read(well_tempering,input,combine("Chrono%sWellTempering",tag));
  };
  
  //setup
  void setup()
  {
    weight.resize(nsweep);
    for(int ipart=0;ipart<npart;ipart++) weight_slice[ipart].resize(nsweep);
    coll.resize(nsweep);
    
    //allocate
    if(width!=0) ngrid=(2*barr+width/2)/width;
    else ngrid=1;
    cout<<"Ngrid: "<<ngrid<<endl;
  }
  
  //do everything
  void mega_reweight(function<double(double)> fun=[](double x){return 1;})
  {
    pote_ave.resize(ngrid+1);
    pote_err.resize(ngrid+1);
    
    if(coeff!=0)
      {
	//reconstruct fractions of it
	for(int ipart=0;ipart<npart;ipart++)
	  {
	    pote_slice[ipart].resize(ngrid+1);
	    reconstruct_potential(pote_slice[ipart],coll,true,npart,ipart); //averaged
	  }
	plot_potential_slices("plots/reco_topo_potential_slice.xmg",pote_slice,npart,false);
	
	//reconstruct the average
	adjust_offset_slices(pote_slice,npart);
	reconstruct_band(pote_ave,pote_err,pote_slice,npart);
	plot_band("plots/reco_topo_potential.xmg",pote_ave,pote_err);
	ofstream umb("plots/umbrella.xmg");
	for(int igrid=0;igrid<=ngrid;igrid++) umb<<-barr+igrid*width<<" "<<pote_ave[igrid]<<endl;
      }
    reconstruct_potential(pote_ave,coll,true);
    
    //always compute weight
    cout<<"Computing reweighting factor"<<endl;
    compute_reweighting_factor(weight,pote_ave,coll,fun);
    double wmax=0;for(int iconf=nterm;iconf<nsweep;iconf++) wmax=max(weight[iconf],wmax); //here just to check
    double t=0;for(int iconf=nterm;iconf<nsweep;iconf++) t+=weight[iconf];
    conf_weight.populate([=](int iconf)->double{return weight[iconf];});
    cout<<"Tot weight: "<<smart_print(conf_weight.ave_err())<<" "<<t/wmax/(nsweep-nterm)<<endl;
    //do the same for slices
    if(coeff!=0)
      for(int ipart=0;ipart<npart;ipart++)
	{
	  compute_reweighting_factor(weight_slice[ipart],pote_slice[ipart],coll,fun);
	  conf_weight_slice[ipart].populate([=](int iconf)->double{return weight_slice[ipart][iconf];});
	  cout<<"Tot weight slice "<<ipart<<": "<<smart_print(conf_weight_slice[ipart].ave_err())<<endl;
	}
    corr_weight.populate([=](int iconf)->double{int icorr=iconf/compute_corr_each;return weight[icorr*compute_corr_each];});
    
    ofstream weight_out("plots/weight.xmg");
    weight_out.precision(12);
    wmax=0;
    for(int iconf=nterm;iconf<nsweep;iconf++) if(fabs(coll[iconf])<barr) wmax=max(wmax,weight[iconf]);
    for(int iconf=nterm;iconf<nsweep;iconf++) /*if(fabs(coll[iconf])<barr)*/ weight_out<<coll[iconf]<<" "<<weight[iconf]/wmax<<endl;
  }
  
  //interpolate using spline
  double interpolate_potential(double x,vector<double> &grid)
  {
    //take igrid
    int igrid=floor((x+barr)/width);
    
    //inside the barriers
    if(igrid>=0 && igrid<ngrid)
      {
	//interpolate
	double x0=igrid*width-barr;
	double m=(grid[igrid+1]-grid[igrid])/width;
	double q=grid[igrid]-m*x0;
	return q+m*x;
      }
    else
      if(igrid<0)
	return force_out*sqr(-x-barr)/2+grid[0];
      else
	return force_out*sqr(+x-barr)/2+grid[ngrid];
  }
  
  //write the potential to a file
  void plot_potential_slices(const char *path,vector<double> *grid,int nslices,bool subt=true,bool ave=false)
  {
    ofstream fout(path);
    for(int islice=0;islice<nslices;islice++)
      {
	double tmax=0;
	for(int igrid=0;igrid<=ngrid;igrid++) tmax=std::max(tmax,grid[islice][igrid]);
	for(double x=-barr-3*width;x<barr+width*3;x+=width/2)
	  {
	    if(ave) fout<<x<<" "<<tmax-(interpolate_potential(x,grid[islice])+interpolate_potential(-x,grid[islice]))/2<<endl;
	    else    fout<<x<<" "<<tmax-interpolate_potential(x,grid[islice])<<endl;
	  }
	fout<<"&"<<endl;
      }
  }
  
  //write band to a file
  void plot_band(const char *path,vector<double> &ave,vector<double> &err)
  {
    vector<double> low(ngrid+1),high(ngrid+1);
    for(size_t i=0;i<ave.size();i++)
      {
	low[i] =ave[i]-err[i];
	high[i]=ave[i]+err[i];
      }
    double m=*max_element(ave.begin(),ave.end());
    
    ofstream fout(path);
    for(double x=-barr-3*width;x<barr+width*3;x+=width/2) fout<<x<<" "<<m-interpolate_potential(x,low)<<endl;
    for(double x=-barr-3*width;x<barr+width*3;x+=width/2) fout<<-x<<" "<<m-interpolate_potential(x,high)<<endl;
  }
  
  //reconstruct the chronological potential
  void reconstruct_potential(vector<double> &out,vector<double> &in,bool ave=false,int npart=1,int ipart=0,int jpart=-1)
  {
    int init_reco_time=time(0);
    
    if(ngrid==1) out[0]=1;
    else
      {
	if(jpart==-1) jpart=ipart+1;
	
	for(auto &grid : out) grid=0;
	
	int part_size=(nsweep-nterm)/npart;
	int start=ipart*part_size+nterm;
	int end=jpart*part_size+nterm;
	int ori_weight=(end-start)/each;
	int weight=ori_weight;
	if((end-start)%each) crash("each %d does not divide (end-start)=%d for part %d/%d, nsweep might be %d",each,end-start,ipart,npart,part_size/each*each*npart+after);
	for(int isweep=after;isweep<end;isweep+=each)
	  {
	    double x=in[isweep];
	    double alpha=x/width;
	    int igrid=floor(alpha)+ngrid/2;
	    alpha=alpha-floor(alpha);
	    if(igrid>=0 && igrid<=ngrid) out[igrid]+=weight*(1-alpha)*coeff/ori_weight;
	    if(igrid+1>=0 && igrid+1<=ngrid) out[igrid+1]+=weight*alpha*coeff/ori_weight;
	    if(ave && isweep>=start) weight=weight-1;
	  }
	
	if(after<nsweep && ave && weight!=0) crash("weight=%d nsweep=%d after=%d each=%d ipart %d npart %d start=%d end=%d",weight,nsweep,after,each,ipart,npart,start,end);
      }
    
    //symmetrize
    for(int i=0;i<=ngrid;i++) out[i]=out[ngrid-i]=(out[i]+out[ngrid-i])/2;
    
    cout<<">>Potential reco time: "<<time(0)-init_reco_time<<" s"<<endl;
  }
  
  //compute the reweighting factor
  void compute_reweighting_factor(vector<double> &weight,vector<double> &grid,vector<double> &coll,function<double(double)> fun)
  {
    double tmax=*max_element(grid.begin(),grid.end());
    
    double tot_weight=0;
    for(int isweep=0;isweep<nsweep;isweep++)
      {
	if(coeff==0 || isweep<after) weight[isweep]=fun(coll[isweep]);
	else
	  {
	    double pot=(interpolate_potential(coll[isweep],grid)+interpolate_potential(-coll[isweep],grid))/2;
	    weight[isweep]=exp(pot-tmax)*fun(coll[isweep]);
	  }
	tot_weight+=weight[isweep];
      }
  }
};

meta_pars_t chrono_topo;
meta_pars_t chrono_charge;

#include "../../src/effmass.cpp"

//compute the effective mass of a function
vector<boot_t> effective_mass(vector<boot_t> in)
{
  vector<boot_t> out(L/2);
  
  for(int t=0;t<L/2;t++)
    {
      auto temp=log(in[t]/in[t+1]);
      auto tempme=temp.ave_err();
      
      for(int iboot=0;iboot<nboot;iboot++) out[t][iboot]=effective_mass(in[t][iboot],in[t+1][iboot],t,L/2,tempme.first,tempme.second);
    }
  
  return out;
}

template<class T> T compute_xi_g(T &n,T &d)
{return sqrt(n/d-1)/(2*sin(M_PI/L));}

//prepare the boot sample
void prepare_boot_sample()
{
  cout<<"Preparing boot samples"<<endl;
  
  //count the number of indipendent configurations
  nind_confs=(nsweep-nterm)/clust_size;
  cout<<"Number of independent configurations: "<<nind_confs<<endl;
  boot_id=new int[nboot*nind_confs];
  
  //draw
  uniform_int_distribution<int> clust_id_dis(0,nind_confs-1);
  uniform_int_distribution<int> clust_entry_dis(0,clust_size-1);
  uniform_int_distribution<int> conf_dis(nterm,nsweep-1);
  mt19937_64 gen(55432334);
  for(int i=0;i<nboot*nind_confs;i++)
    {
      int clust_id=clust_id_dis(gen);
      //int clust_entry=clust_entry_dis(gen);
      //int conf_id=nterm+clust_entry+clust_size*clust_id;
      ////int conf_id=conf_dis(gen);
      //if(conf_id>=nsweep) crash("conf_id %d >= nsweep %d -  clust_size %d - clust_entry %d - clust_id %d - nind_confs %d",
      //conf_id,nsweep,clust_size,clust_entry,clust_id,nind_confs);
      //boot_id[i]=conf_id;
      boot_id[i]=clust_id;
    }
}

inline int get_boot_id(int iboot,int ientry)
{return boot_id[ientry+nind_confs*iboot];}

//compute average and error according to boostrap
void boot_t::populate(function<double(int)> fun)
{
  //compute per clust
  double *per_clust=new double[nind_confs];
  for(int clust_id=0;clust_id<nind_confs;clust_id++)
    {
      per_clust[clust_id]=0;
      for(int clust_entry=0;clust_entry<clust_size;clust_entry++)
	{
	  int iconf_id=nterm+clust_entry+clust_size*clust_id;
	  per_clust[clust_id]+=fun(iconf_id);
	}
      per_clust[clust_id]/=clust_size;
    }
  
  //compute per boot
  for(int iboot=0;iboot<nboot;iboot++)
    {
      double per_boot_sum=0;
      for(int ientry=0;ientry<nind_confs;ientry++)
	{
	  int clust_id=get_boot_id(iboot,ientry);
	  per_boot_sum+=per_clust[clust_id];
	}
      per_boot_sum/=nind_confs;
      (*this)[iboot]=per_boot_sum;
    }
  
  delete[] per_clust;
}

//compute average and error
ave_err_t boot_t::ave_err()
{
  double ave=0,err=0;
  for(auto &it : *this)
    {
      ave+=it;
      err+=it*it;
    }
  
  ave/=nboot;
  err/=nboot;
  err-=ave*ave;
  err=sqrt(err);
  
  return make_pair(ave,err);
}

inline int ipart_fun(int iconf)
{return (iconf-nterm)/part_size;}

ave_err_t slice_ave_err(boot_t *in)
{
  double ave=0,err=0,ave_err=0;
  for(int ipart=0;ipart<npart;ipart++)
    {
      ave_err_t temp=in[ipart].ave_err();
      ave+=temp.first;
      ave_err+=temp.second*temp.second;
      err+=temp.first*temp.first;
    }
  
  ave/=npart;
  ave_err/=npart;
  err/=npart;
  
  err-=ave*ave;
  err/=npart-1;
  
  ave_err_t out;
  out.first=ave;
  out.second=sqrt(ave_err+err);
  
  return out;
}

//bayesian average and erorr of binomial
double bayes_bin_ave(int k,int n)
{return (1.0+k)/(2.0+n);}

//bayesian average and erorr of binomial
double bayes_bin_err(int k,int n)
{return sqrt((1.0+k)*(1.0-k+n)/((2.0+n)*(2.0+n)*(3.0+n)));}

template <class T> vector<T> lin_solve(vector<double> &A,vector<T> &b)
{
  int d=b.size();
  vector<T> x(d);
  
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
      T S=0;
      for(int i=k+1;i<d;i++) S+=A[k*d+i]*x[i];
      x[k]=b[k]-S;
    }
  
  return x;
}

template <class T> class poly_fit_res_t
{
public:
  int ndof;
  vector<T> coeffs;
  T chi2;
};

template <class T> T poly(double x,vector<T> coeffs)
{
  T out=0;
  for(size_t iord=0;iord<coeffs.size();iord++) out+=pow(x,iord)*coeffs[iord];
  return out;
}

template<class T> T poly_fit_chi2(vector<double> &x,vector<T> &y,vector<double> &e,vector<T> &coeffs,double xmin=-1e+300,double xmax=+1e+300)
{
  T chi2=0;
  for(size_t p=0;p<y.size();p++)
    if((x[p]<=xmax && x[p]>=xmin)||(x[p]>=-xmax && x[p]<=-xmin)) chi2+=sqr((poly(x[p],coeffs)-y[p])/e[p]);
  return chi2;
}

template <class T> poly_fit_res_t<T> poly_fit(vector<double> &x,vector<T> &y,vector<double> &e,int d,double xmin=-1e+300,double xmax=+1e+300)
{
  int np=y.size();
  
  vector<double> Al(2*d+1);
  for(auto &it :Al) it=0.0;
  vector<T> c(d+1);
  for(auto &it : c) it=0.0;
  
  int npoints_used=0;
  for(int p=0;p<np;p++)
    if((x[p]<=xmax && x[p]>=xmin)||(x[p]>=-xmax && x[p]<=-xmin))
      {
	npoints_used++;
	
        //calculate the weight
        double w=pow(e[p],-2);
        //compute Al and c
        for(int f=0;f<=2*d;f++)
          {
            Al[f]+=w;
            if(f<=d) c[f]+=y[p]*w;
            w*=x[p];
          }
      }
  
  vector<double> A((d+1)*(d+1));
  for(int i=0;i<=d;i++)
    for(int j=0;j<=d;j++)
      A[i*(d+1)+j]=Al[i+j];
  
  //build output
  poly_fit_res_t<T> out;
  out.coeffs=lin_solve(A,c);
  out.chi2=poly_fit_chi2(x,y,e,out.coeffs,xmin,xmax);
  out.ndof=npoints_used-(d+1);
  
  return out;
}
template <class T> poly_fit_res_t<T> poly_fit(vector<double> &x,vector<T> &y,int d,double xmin=-1e+300,double xmax=+1e+300)
{
  vector<double> e;
  for(auto &it : y) e.push_back(it.err());
  
  return poly_fit(x,y,e,d,xmin,xmax);
}

/////////////////////////////////////////////////////// read //////////////////////////////////////////////////

//scan the whole input file
void read_input(bool read_charge_flag)
{
  ifstream input("input");
  if(!input.good()) crash("error opening input");
  read(N,input,"N");
  read(L,input,"L");
  read(beta,input,"Beta");
  g=1/(N*beta);
  cout<<"g: "<<g<<endl;
  read(seed,input,"Seed");
  read(nsweep,input,"NSweep");
  string start_cond_str;
  read(start_cond_str,input,"StartCond");
  read(nterm,input,"NTerm");
  part_size=(nsweep-nterm)/npart;
  if(part_size*npart+nterm!=nsweep) crash("npart=%d does not divide (nsweeps-nterm)=%d",npart,nsweep-nterm);
  read(compute_corr_each,input,"ComputeCorrEach");
  ncorr=nsweep/compute_corr_each;
  if(ncorr*compute_corr_each!=nsweep) crash("nsweep=%d not multiple of compute_corr_each=%d",nsweep,compute_corr_each);
  read(use_hmc,input,"UseHMC");
  if(!use_hmc) crash("Have to use HMC");
  read(nhmc_steps,input,"NhmcSteps");
  read(use_topo_pot,input,"UseTopoPot");
  switch(use_topo_pot)
    {
    case 2:
      chrono_topo.read_pars(input,"Topo");
      break;
    default:
      chrono_topo.coeff=0;
      break;
    }
  read(nstout_lev,input,"NStoutLev");
  
  read(stout_rho,input,"StoutRho");
  if(read_charge_flag)
    {
      read(use_charge_pot,input,"UseChargePot");
      read(charge_pot,input,"ChPot");
      cout<<"1/g: "<<1/g<<", S: "<<sinh(charge_pot/L)/g<<endl;
      switch(use_charge_pot)
	{
	case 2:
	  chrono_charge.read_pars(input,"Charge");
	  break;
	default:
	  chrono_charge.coeff=0;
	  break;
	}
    }
  
  input.close();
}

//read energy
void read_energy()
{
  int a=0,n=0;
  cout<<"Reading energy"<<endl;
  ene.resize(nsweep);
  FILE *finene=open_file("energy","r");
  for(int isweep=0;isweep<nsweep;isweep++)
    {
      int iconf;
      double e;
      if(fscanf(finene,"%d %lg",&iconf,&e)!=2) crash("reading energy for iconf %d",isweep);
      
      if(iconf!=isweep) crash("reading energy for iconf %d, expected to be %d",iconf,isweep);
      
      ene[isweep]=e;
      
      if(isweep && isweep>=nterm)
	{
	  n++;
	  if(e!=ene[isweep-1]) a++;
	}
    }
  cout<<"Acceptance: "<<(double)a/n<<endl;
  
  fclose(finene);
}

//read magnetization
void read_magnetization()
{
  cout<<"Reading magnetization"<<endl;
  mag[0].resize(ncorr);
  mag[1].resize(ncorr);
  FILE *finmag=open_file("mag","r");
  for(int icorr=0;icorr<ncorr;icorr++)
    {
      int iconf;
      double mr,mi;
      if(fscanf(finmag,"%d %lg %lg",&iconf,&mr,&mi)!=3) crash("reading magnetization for icorr %d",icorr);
      
      if(iconf!=icorr*compute_corr_each) crash("reading magnetization for iconf %d, expected to be %d",iconf,icorr*compute_corr_each);
      
      mag[0][icorr]=mr;
      mag[1][icorr]=mi;
    }
  fclose(finmag);
}

//read correlation
void read_correlation()
{
  corr=new vector<double>[L/2+1];
  for(int i=0;i<L/2+1;i++) corr[i].resize(ncorr);
  
  cout<<"Reading correlation"<<endl;
  FILE *fincorr=open_file("corr","r");
  for(int icorr=0;icorr<ncorr;icorr++)
    for(int t=0;t<L/2+1;t++)
      {
	int isweep=icorr*compute_corr_each;
	int iconf,u;
	double c;
	if(fscanf(fincorr,"%d %d %lg",&iconf,&u,&c)!=3) crash("reading correlation for icorr %d (conf %d), t=%d",icorr,isweep,u);
	
	if(iconf!=isweep) crash("with iconf reading correlation for iconf %d, t %d, expected to be %d",iconf,u,isweep);
	if(u!=t) crash("with t reading correlation for iconf %d, t %d, expected to be %d",iconf,u,t);
      
      corr[t][icorr]=c;
      }
  fclose(fincorr);
}

