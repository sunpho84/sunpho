#include <algorithm>
#include <cstdlib>
#include <cstdarg>
#include <cmath>
#include <fftw3.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <sstream>

#define crash(...) internal_crash(__LINE__,__VA_ARGS__)

using namespace std;
using namespace placeholders;

//============================== types ============================

const int predef_nstout_lev=2;

enum start_cond_t{COLD,HOT,LOAD};

//returned as ave_err
typedef pair<double,double> boot_res_t;
ostream& operator<<(ostream &out,boot_res_t obj)
{return out<<obj.first<<" "<<obj.second;}
string smart_print(boot_res_t in)
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

int clust_size=1000;
int nboot=1000;
struct boot_t : vector<double>
{
  void init(){resize(nboot);}
  boot_t(){init();}
  boot_t(double in){init();(*this)=in;}
  void populate(function<double(int)>);
  void populate(function<double(int,std::vector<double>&)>);
  boot_t(function<double(int)> fun){init();populate(fun);}
  boot_t(function<double(int,std::vector<double>&)> fun){init();populate(fun);}
  boot_res_t ave_err();
  boot_t operator+(boot_t &oth)
  {boot_t out;transform(this->begin(),this->end(),oth.begin(),out.begin(),plus<double>());return out;}
  boot_t operator-(boot_t &oth)
  {boot_t out;transform(this->begin(),this->end(),oth.begin(),out.begin(),minus<double>());return out;}
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

//============================== var ==============================

int seed;
start_cond_t start_cond;
int nterm,nsweep,ncorr,nmicro,use_hmc;
int N;
double beta;
int L;
int nstout_lev;
double stout_rho;
int nhmc_steps;
int use_topo_pot;
double th_top;
int chrono_topo_after;
int chrono_topo_each;
double chrono_topo_coeff;
double chrono_topo_width;
double chrono_topo_barr;
double chrono_topo_force_out;
double chrono_topo_well_tempering;
double chrono_topo_bend;
int compute_corr_each;

vector<double> mag[2];
vector<double> ene;
vector<double> *corr;
vector<double> non_geo_Q[predef_nstout_lev+1];
vector<int> geo_Q;

const int npart=2;
vector<double> topo_pote_slice[npart];
vector<double> weight[npart];

int max_geo_Q;
int min_geo_Q;

int nind_confs;
int *boot_id;

int ngrid;

//=========================== tools =================================

boot_t get_reweighted(function<double(int)> fun,int frac=1)
{
  boot_t out;
  for(int ipart=1;ipart<npart;ipart++)
    {
      boot_t temp([frac,ipart,fun](int iconf)->double{return fun(iconf)*weight[ipart][iconf/frac];});
      for(int iboot=ipart-1;iboot<nboot;iboot+=npart-1) out[iboot]=temp[iboot];
    }
  
  return out;
}

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

//scan the whole input file
void read_input()
{
  ifstream input("input");
  if(!input.good()) crash("error opening input");
  read(N,input,"N");
  read(L,input,"L");
  read(beta,input,"Beta");
  read(seed,input,"Seed");
  read(nsweep,input,"NSweep");
  string start_cond_str;
  read(start_cond_str,input,"StartCond");
  read(nterm,input,"NTerm");
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
      read(chrono_topo_after,input,"ChronoTopoAfter");
      read(chrono_topo_each,input,"ChronoTopoEach");
      read(chrono_topo_coeff,input,"ChronoTopoCoeff");
      read(chrono_topo_width,input,"ChronoTopoWidth");
      read(chrono_topo_barr,input,"ChronoTopoBarr");
      read(chrono_topo_force_out,input,"ChronoTopoForceOut");
      read(chrono_topo_bend,input,"ChronoTopoBend");
      read(chrono_topo_well_tempering,input,"ChronoTopoWellTempering");
      break;
    default:
      chrono_topo_coeff=0;
      break;
    }
  read(nstout_lev,input,"NStoutLev");
  if(nstout_lev!=predef_nstout_lev) crash("Compiled for nstout_lev=%d",predef_nstout_lev);
  read(stout_rho,input,"StoutRho");

  input.close();
}

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

//read all the data files
void read_data()
{
  int init_read_time=time(0);

  //allocate data
  for(int ipart=0;ipart<npart;ipart++) weight[ipart].resize(nsweep);
  for(int ilev=0;ilev<=predef_nstout_lev;ilev++) non_geo_Q[ilev].resize(nsweep);
  geo_Q.resize(nsweep);
  ene.resize(nsweep);
  mag[0].resize(ncorr);
  mag[1].resize(ncorr);
  corr=new vector<double>[L/2+1];
  for(int i=0;i<L/2+1;i++) corr[i].resize(ncorr);
  
  ///////////////////////////////////////////////////////////
  
  //read energy
  cout<<"Reading energy"<<endl;
  FILE *finene=open_file("energy","r");
  for(int isweep=0;isweep<nsweep;isweep++)
    {
      int iconf;
      double e;
      if(fscanf(finene,"%d %lg",&iconf,&e)!=2) crash("reading energy for iconf %d",isweep);
      
      if(iconf!=isweep) crash("reading energy for iconf %d, expected to be %d",iconf,isweep);
      
      ene[isweep]=e;
    }
  fclose(finene);
  
  ///////////////////////////////////////////////////////////
  
  //read magnetization
  cout<<"Reading magnetization"<<endl;
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
  
  ///////////////////////////////////////////////////////////
  
  //read correlation
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
  
  ///////////////////////////////////////////////////////////
  
  //read topology
  cout<<"Reading topology"<<endl;
  FILE *fintop=open_file("topology","r");
  
  min_geo_Q=+100000;
  max_geo_Q=-100000;
  for(int isweep=0;isweep<nsweep;isweep++)
    for(int jlev=0;jlev<=nstout_lev;jlev++)
      {
	int iconf,ilev;
	double geo,non_geo;
	//if(fscanf(fintop,"%d %d %lg %lg",&iconf,&ilev,&geo,&non_geo)!=4) crash("reading topo at isweep %d ilev %d",isweep,ilev);
	char siconf[16],silev[3],sgeo[20],snon_geo[20];
	if(fscanf(fintop,"%s %s %s %s",siconf,silev,sgeo,snon_geo)==0) crash("reading topo at isweep %d ilev %",isweep,jlev);
	iconf=atoi(siconf);
	ilev=atoi(silev);
	geo=strtod(sgeo,NULL);
	non_geo=strtod(snon_geo,NULL);
	
	if(isweep!=iconf) crash("isweep=%d!=iconf=%d",isweep,iconf);
	if(ilev!=jlev) crash("jlev=%d!=ilev=%d",jlev,ilev);
	
	//approx to closest integer the geometrical
	geo_Q[isweep]=round(geo);
	non_geo_Q[ilev][isweep]=non_geo;
	
	//find min/max
	min_geo_Q=min(min_geo_Q,geo_Q[isweep]);
	max_geo_Q=max(max_geo_Q,geo_Q[isweep]);
      }
  fclose(fintop);
  
  cout<<">>Reading time: "<<time(0)-init_read_time<<" s"<<endl;
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
boot_res_t boot_t::ave_err()
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

//bayesian average and erorr of binomial
double bayes_bin_ave(int k,int n)
{return (1.0+k)/(2.0+n);}

//bayesian average and erorr of binomial
double bayes_bin_err(int k,int n)
{return sqrt((1.0+k)*(1.0-k+n)/((2.0+n)*(2.0+n)*(3.0+n)));}

//================================= analysis ================================

//compute the autocorrelation time of the topological charge
double compute_topo_tint(const char *path,auto &Q)
{
  fftw_complex *in=(fftw_complex*)fftw_malloc((nsweep-nterm)*sizeof(fftw_complex));
  fftw_complex *out=(fftw_complex*)fftw_malloc((nsweep-nterm)*sizeof(fftw_complex));
  
  //fill and take fft
  for(int i=0;i<nsweep-nterm;i++)
    {
      in[i][0]=Q[i+nterm];
      in[i][1]=0;
    }
  fftw_plan pf=fftw_plan_dft_1d(nsweep-nterm,in,out,FFTW_FORWARD,FFTW_ESTIMATE);  
  fftw_execute(pf);
  fftw_destroy_plan(pf);
  
  //take the square, remove zero mode and put normalization
  double norm=1.0/sqr(nsweep-nterm);
  for(int i=0;i<nsweep-nterm;i++)
    {
      in[i][0]=(out[i][0]*out[i][0]+out[i][1]*out[i][1])*norm;
      in[i][1]=0;
    }
  in[0][0]=in[0][1]=0;
    
  //compute back
  fftw_plan pb=fftw_plan_dft_1d(nsweep-nterm,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(pb);
  fftw_destroy_plan(pb);
  
  //print the corr
  ofstream out_corr(path);
  for(int i=0;i<(nsweep-nterm)/2;i++) out_corr<<out[i][0]/out[0][0]<<endl;
  
  free(in);
  free(out);
  
  return 0;
}

//compute the autocorrelation time of the mag
double compute_mag_tint(const char *path,auto &M)
{
  int n=nsweep/compute_corr_each;
  n-=n%3;
  fftw_complex *in=(fftw_complex*)fftw_malloc(2*n/3*sizeof(fftw_complex));
  fftw_complex *out=(fftw_complex*)fftw_malloc(2*n/3*sizeof(fftw_complex));
					       
  //fill and take fft
  for(int i=n/3;i<n;i++)
    {
      in[i-n/3][0]=M[i];
      in[i-n/3][1]=0;
    }
  fftw_plan pf=fftw_plan_dft_1d(2*n/3,in,out,FFTW_FORWARD,FFTW_ESTIMATE);  
  fftw_execute(pf);
  fftw_destroy_plan(pf);
  
  //take the square, remove zero mode and put normalization
  double norm=1.0/sqr(2*n/3);
  for(int i=0;i<2*n/3;i++)
    {
      in[i][0]=(out[i][0]*out[i][0]+out[i][1]*out[i][1])*norm;
      in[i][1]=0;
    }
  in[0][0]=in[0][1]=0;
    
  //compute back
  fftw_plan pb=fftw_plan_dft_1d(2*n/3,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(pb);
  fftw_destroy_plan(pb);
  
  //print the corr
  ofstream out_corr(path);
  for(int i=0;i<2*n/3/2;i++) out_corr<<i*compute_corr_each<<" "<<out[i][0]/out[0][0]<<endl;
  
  free(in);
  free(out);
  
  return 0;
}

//allocate and set x
void setup_chrono_topo_potential()
{
  //allocate
  ngrid=(2*chrono_topo_barr+chrono_topo_width/2)/chrono_topo_width;
  cout<<"Ngrid: "<<ngrid<<endl;
}

//reconstruct the chronological topological potential
void reconstruct_chrono_topo_potential(vector<double> &out,bool ave=false,int npart=1,int ipart=0,int jpart=-1)
{
  if(jpart==-1) jpart=ipart+1;
  
  for(auto &grid : out) grid=0;
  
  int init_reco_time=time(0);
  
  int part_size=(nsweep-chrono_topo_after)/npart;
  int start=ipart*part_size+chrono_topo_after;
  int end=jpart*part_size+chrono_topo_after;
  int ori_weight=(end-start)/chrono_topo_each;
  int weight=ori_weight;
  if((end-start)%chrono_topo_each) crash("chrono_topo_each %d does not divide (end-start)=%d for part %d/%d, nsweep might be %d",
			 chrono_topo_each,end-start,ipart,npart,part_size/chrono_topo_each*chrono_topo_each*npart+chrono_topo_after);
  for(int isweep=chrono_topo_after;isweep<end;isweep+=chrono_topo_each)
    {
      double Q=non_geo_Q[nstout_lev][isweep];
      
      int igrid=floor(Q/chrono_topo_width)+ngrid/2;
      double alpha=Q/chrono_topo_width;
      alpha=alpha-floor(alpha);
      if(igrid>=0 && igrid<=ngrid) out[igrid]+=weight*(1-alpha)*chrono_topo_coeff/ori_weight;
      if(igrid+1>=0 && igrid+1<=ngrid) out[igrid+1]+=weight*alpha*chrono_topo_coeff/ori_weight;
      if(ave && isweep>=start) weight=weight-1;
    }
  
  if(ave && weight!=0) crash("weight=%d nsweep=%d chrono_topo_after=%d each=%d ipart %d npart %d start=%d end=%d",weight,nsweep,chrono_topo_after,chrono_topo_each,ipart,npart,start,end);

  cout<<">>Potential reco time: "<<time(0)-init_reco_time<<" s"<<endl;
}

//compute the average and the error
void reconstruct_band_topo_potential(vector<double> &ave,vector<double> &err,vector<double> *slices,int npart)
{
  //reset
  for(int i=0;i<=ngrid;i++) ave[i]=err[i]=0;
  for(int islice=1;islice<npart;islice++)
    {
      double a=0;
      for(int i=0;i<=ngrid;i++) a+=slices[islice][i];
      a/=ngrid+1;
      for(int i=0;i<=ngrid;i++) slices[islice][i]-=a;
    }

  for(int islice=1;islice<npart;islice++)
    for(int i=0;i<=ngrid;i++)
      {
	double a=(slices[islice][i]+slices[islice][ngrid-i])/2;
	ave[i]+=a;
	err[i]+=a*a;
      }
  for(int i=0;i<=ngrid;i++)
    {
      ave[i]/=npart-1;
      err[i]/=npart-1;
      err[i]-=ave[i]*ave[i];
      err[i]=sqrt(err[i]);
    }
}

//interpolate using spline
double interpolate_topopotential(double Q,vector<double> &topo_grid)
{
  //take igrid
  int igrid=floor((Q+chrono_topo_barr)/chrono_topo_width);
  
  //inside the barriers
  if(igrid>=0 && igrid<ngrid)
    {
      //interpolate
      double x0=igrid*chrono_topo_width-chrono_topo_barr;
      double m=(topo_grid[igrid+1]-topo_grid[igrid])/chrono_topo_width;
      double q=topo_grid[igrid]-m*x0;
      return q+m*Q;
    }
  else
    if(igrid<0)
      return chrono_topo_force_out*sqr(-Q-chrono_topo_barr)/2+topo_grid[0];
    else
      return chrono_topo_force_out*sqr(+Q-chrono_topo_barr)/2+topo_grid[ngrid];
}

//write the potential to a file
void plot_chrono_topo_potential_slices(const char *path,vector<double> *topo_grid,int nslices,bool subt=true,bool ave=false)
{
  ofstream fout(path);
  for(int islice=0;islice<nslices;islice++)
    {
      double tmax=0;
      for(int igrid=0;igrid<=ngrid;igrid++) tmax=std::max(tmax,topo_grid[islice][igrid]);
      for(double Q=-chrono_topo_barr-3*chrono_topo_width;Q<chrono_topo_barr+chrono_topo_width*3;Q+=chrono_topo_width/2) 
	{
	  if(ave) fout<<Q<<" "<<tmax-(interpolate_topopotential(Q,topo_grid[islice])+interpolate_topopotential(-Q,topo_grid[islice]))/2<<endl;
	  else    fout<<Q<<" "<<tmax-interpolate_topopotential(Q,topo_grid[islice])<<endl;
	}
      fout<<"&"<<endl;
    }
}

//write band to a file
void plot_chrono_topo_potential_band(const char *path,vector<double> &ave,vector<double> &err)
{
  vector<double> low(ngrid+1),high(ngrid+1);
  transform(ave.begin(),ave.end(),err.begin(),low.begin(),minus<double>());
  transform(ave.begin(),ave.end(),err.begin(),high.begin(),plus<double>());
  double m=*max_element(ave.begin(),ave.end());
  
  ofstream fout(path);
  for(double Q=-chrono_topo_barr-3*chrono_topo_width;Q<chrono_topo_barr+chrono_topo_width*3;Q+=chrono_topo_width/2) 
    fout<<Q<<" "<<m-interpolate_topopotential(Q,low)<<endl;
  for(double Q=-chrono_topo_barr-3*chrono_topo_width;Q<chrono_topo_barr+chrono_topo_width*3;Q+=chrono_topo_width/2) 
    fout<<-Q<<" "<<m-interpolate_topopotential(Q,high)<<endl;
}

//compute the reweighting factor
void compute_reweigting_factor()
{
  cout<<"Computing reweighting factor"<<endl;

  for(int ipart=1;ipart<npart;ipart++)
    {
      double tmax=*max_element(topo_pote_slice[ipart].begin(),topo_pote_slice[ipart].end());
      
      for(int isweep=0;isweep<nsweep;isweep++)
	if(chrono_topo_coeff==0 || isweep<chrono_topo_after) weight[ipart][isweep]=1;
	else
	  {
	    double pot=(interpolate_topopotential(non_geo_Q[nstout_lev][isweep],topo_pote_slice[ipart])+
			interpolate_topopotential(-non_geo_Q[nstout_lev][isweep],topo_pote_slice[ipart]))/2;
	    weight[ipart][isweep]=exp(pot-tmax);
	  }
    }
}

//compute bare histogram of Q and reweighted ones
void draw_Q_histograms()
{
  for(int rew=0;rew<2;rew++)
    {
      map<int,boot_t> histo_Q;
      for(int Q=min_geo_Q;Q<=max_geo_Q;Q++)
	if(rew) histo_Q[Q]=get_reweighted([Q](int isweep)->double{return (geo_Q[isweep]==Q);});
	else    histo_Q[Q].populate([Q](int isweep)->double{return (geo_Q[isweep]==Q);});
      
      //normalize factor
      boot_t norm=0;
      for(auto &it : histo_Q) norm+=it.second;
      for(auto &it : histo_Q) it.second/=norm;

      //write the histo
      ofstream histo_geo_Q(rew?"plots/histo_geo_Q_rew.xmg":"plots/histo_geo_Q.xmg");
      histo_geo_Q<<"@type xydy"<<endl;
      for(auto it : histo_Q) histo_geo_Q<<it.first<<" "<<it.second.ave_err()<<endl;
    }
}

//compute the non-integer part of the noise
void plot_noise(double *Q_reno)
{
  cout<<">>Writing noise"<<endl;

  ofstream out_sigma("plots/sigma_noise.xmg");
  
  //allocate noise and plan
  fftw_complex *noise=(fftw_complex*)fftw_malloc((nsweep-nterm)*sizeof(fftw_complex));
  fftw_plan pf=fftw_plan_dft_1d(nsweep-nterm,noise,noise,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_plan pb=fftw_plan_dft_1d(nsweep-nterm,noise,noise,FFTW_BACKWARD,FFTW_ESTIMATE);

  //extract the noise from each corr
  FILE *fout_noise=open_file("plots/noise.xmg","w");
  FILE *fout_power=open_file("plots/noise_power_spectrum.xmg","w");
  FILE *fout_corre=open_file("plots/noise_correlation.xmg","w");
  for(int ilev=0;ilev<=nstout_lev;ilev++)
    {
      //compute the noise and print it
      for(int isweep=nterm;isweep<nsweep;isweep++)
	{
	  double n=geo_Q[isweep]-non_geo_Q[ilev][isweep]/Q_reno[ilev];
	  while(n<-0.5) n+=1;
	  while(n>+0.5) n-=1;
	  noise[isweep-nterm][0]=n;
	  noise[isweep-nterm][1]=0;
	  fprintf(fout_noise,"%lg\n",n);
	}
      fprintf(fout_noise,"&\n");
      
      //take fft and print
      fftw_execute(pf);
      for(int isweep=0;isweep<nsweep-nterm;isweep++)
	{
	  double n2=sqr(noise[isweep][0])+sqr(noise[isweep][1]);
	  noise[isweep][0]=n2/(nsweep-nterm);
	  noise[isweep][1]=0;
	  
	  if(isweep<2000) fprintf(fout_power,"%lg\n",n2);
	}
      fprintf(fout_power,"&\n");
            
      //take back fft and print
      fftw_execute(pb);
      for(int isweep=0;isweep<2000;isweep++) fprintf(fout_corre,"%lg\n",noise[isweep][0]/noise[0][0]);
      fprintf(fout_corre,"&\n");
      
      out_sigma<<ilev<<" "<<sqrt(noise[0][0]/(nsweep-nterm))<<endl;
    }

  out_sigma.close();
  
  fclose(fout_noise);
  fclose(fout_power);
  fclose(fout_corre);
  fftw_destroy_plan(pf);
  fftw_destroy_plan(pb);
  free(noise);
}

int main()
{
  int init_time=time(0);
  read_input();
  
  //read the data and reconstruct reweigthing factor
  read_data();
  cout<<"NSweeps: "<<nsweep<<endl;
  
  //setup x y and minmax for chrono topo potential
  setup_chrono_topo_potential();
      
  vector<double> topo_pote_ave(ngrid+1),topo_pote_err(ngrid+1);
  if(chrono_topo_coeff!=0)
    {
      //reconstruct fractions of it
      for(int ipart=0;ipart<npart;ipart++)
	{
	  topo_pote_slice[ipart].resize(ngrid+1);
	  reconstruct_chrono_topo_potential(topo_pote_slice[ipart],true,npart,ipart); //averaged
	}
      plot_chrono_topo_potential_slices("plots/reco_topo_potential_slice.xmg",topo_pote_slice,npart,false);
      
      //reconstruct the average
      reconstruct_band_topo_potential(topo_pote_ave,topo_pote_err,topo_pote_slice,npart);
      plot_chrono_topo_potential_band("plots/reco_topo_potential.xmg",topo_pote_ave,topo_pote_err);
    }
  
  //always compute weight
  compute_reweigting_factor();
  
  //compute autocorrelation time
  compute_topo_tint("plots/topo_correlation.xmg",geo_Q);
  compute_topo_tint("plots/nongeo_topo_correlation.xmg",non_geo_Q[nstout_lev]);
  
  //draw random indices
  prepare_boot_sample();

  boot_t tot_weight=get_reweighted([](int iconf)->double{return 1;});
  boot_t tot_weight_corr=get_reweighted([](int iconf)->double{return 1;},compute_corr_each);
  cout<<"Tot weight: "<<smart_print(tot_weight.ave_err())<<endl;
  cout<<"Tot weight_corr: "<<smart_print(tot_weight_corr.ave_err())<<endl;

  //===================== now analysis can start ==================
  
  draw_Q_histograms();

  //write the topology
  if(nsweep<1000000)
    {
      cout<<"Writing geo topology"<<endl;
      ofstream out_topo("plots/geo_topology.xmg");
      for(int i=nterm;i<nsweep;i++) out_topo<<geo_Q[i]<<endl;
      ofstream out_topo_non_geo("plots/non_geo_topology.xmg");
      for(int i=nterm;i<nsweep;i++) out_topo_non_geo<<non_geo_Q[nstout_lev][i]<<endl;
    }
  
  //compute average renormalization factors and Q2
  double Q_reno[nstout_lev+1];
  boot_t geo_Q1([](int iconf)->double{return geo_Q[iconf];});
  boot_t geo_Q2([](int iconf)->double{return sqr(geo_Q[iconf]);});
  boot_t geo_Q1_rew=get_reweighted([](int iconf)->double{return geo_Q[iconf];})/tot_weight;
  boot_t geo_Q2_rew=get_reweighted([](int iconf)->double{return sqr(geo_Q[iconf]);})/tot_weight;
  boot_t geo_Q4([](int iconf)->double{return sqr(sqr(geo_Q[iconf]));});
  boot_t geo_Q6([](int iconf)->double{return pow(geo_Q[iconf],6);});
  cout<<"Geo topo susceptibility: "<<smart_print((geo_Q2-geo_Q1*geo_Q1).ave_err())<<endl;
  cout<<"Geo topo susceptibility rew: "<<smart_print((geo_Q2_rew-geo_Q1_rew*geo_Q1_rew).ave_err())<<endl;
  cout<<"Geo topo chi4: "<<((geo_Q4-3*sqr(geo_Q2))/geo_Q2).ave_err()<<endl;
  cout<<"Geo topo chi6: "<<((geo_Q6-15*geo_Q2*geo_Q4+30*pow(geo_Q2,3))/geo_Q2).ave_err()<<endl;
  for(int ilev=0;ilev<=nstout_lev;ilev++)
    {
      boot_t geo_Q_non_geo_Q([&](int iconf)->double{return geo_Q[iconf]*non_geo_Q[ilev][iconf];});
      boot_res_t reno=(geo_Q_non_geo_Q/geo_Q2).ave_err();
      Q_reno[ilev]=reno.first;
      cout<<"Renormalization constant["<<ilev<<"]: "<<smart_print(reno)<<"   "<<endl;
    }
  
  //compute noise variance
  for(int ilev=0;ilev<=nstout_lev;ilev++)
    {
      double z=(std::isnan(Q_reno[ilev])?1:Q_reno[ilev]);
      boot_t nonint([z,ilev](int iconf)->double{return non_geo_Q[ilev][iconf]/z-geo_Q[iconf];});
      boot_t nonint2([z,ilev](int iconf)->double{return sqr(non_geo_Q[ilev][iconf]/z-geo_Q[iconf]);});
      cout<<"Noise variance["<<ilev<<"]: "<<smart_print(sqrt(nonint2-nonint*nonint).ave_err())<<"   "<<endl;
    }
  
  //compute occupancy per topological sector
  vector<boot_t> occupancy(max_geo_Q-min_geo_Q+1);
  boot_t total_occupancy=0;
  ofstream occupancy_out("plots/occupancy_ps.xmg");
  occupancy_out<<"@type xydy"<<endl;
  for(int iQ=min_geo_Q;iQ<=max_geo_Q;iQ++)
    {
      int jQ=iQ-min_geo_Q;
      occupancy[jQ].populate([iQ](int iconf)->int{return (geo_Q[iconf]==iQ);});
      total_occupancy+=occupancy[jQ];
      auto occupancy_ave_err=occupancy[jQ].ave_err();
      if(!std::isnan(occupancy_ave_err.second)) occupancy_out<<iQ<<" "<<occupancy_ave_err<<endl;
    }
  cout<<"Check, total occupancy: "<<smart_print(total_occupancy.ave_err())<<endl;

  ///compute the observables sector per sector
  boot_t ene_ps;
  vector<boot_t> mag_ps(max_geo_Q-min_geo_Q+1);
  vector<boot_t> xi_g_ps(max_geo_Q-min_geo_Q+1);
  ofstream ene_ps_out("plots/ene_ps.xmg");
  ofstream mag_ps_out("plots/mag_ps.xmg");
  ofstream mag_ps_rescaled_out("plots/mag_ps_rescaled.xmg");
  ofstream xi_g_ps_out("plots/xi_g_ps.xmg");
  ofstream out_eff_ps("plots/effective_mass_ps.xmg");
  ene_ps_out<<"@type xydy"<<endl;
  mag_ps_out<<"@type xydy"<<endl;
  mag_ps_rescaled_out<<"@type xydy"<<endl;
  xi_g_ps_out<<"@type xydy"<<endl;
  out_eff_ps<<"@type xydy"<<endl;
  for(int iQ=min_geo_Q;iQ<=max_geo_Q;iQ++)
    {
      int jQ=iQ-min_geo_Q;
      int kQ=-iQ-min_geo_Q;
      
      //compute occupancy
      boot_t occ([iQ](int iconf)->int{return (geo_Q[iconf/compute_corr_each*compute_corr_each]==iQ);});
      
      //compute energy and normalize
      boot_t ene_ps([iQ](int iconf)->double{return ene[iconf]*(geo_Q[iconf]==iQ);});
      ene_ps/=occ;
      auto ene_ps_ave_err=ene_ps.ave_err();
      if(!std::isnan(ene_ps_ave_err.second)) ene_ps_out<<iQ<<" "<<ene_ps_ave_err<<endl;
      
      //compute correlation function
      vector<boot_t> correlation_ps(L/2+1);
      for(int t=0;t<L/2+1;t++)
	{
	  correlation_ps[t].populate([iQ,t](int iconf)->double{return corr[t][iconf/compute_corr_each]*
		(geo_Q[iconf/compute_corr_each*compute_corr_each]==iQ);});
	  correlation_ps[t]/=occ;
	}
      auto effective_mass_ps=effective_mass(correlation_ps);
      for(int t=0;t<L/2;t++)
	{
	  auto effective_mass_ave_err=effective_mass_ps[t].ave_err();
	  if(!std::isnan(effective_mass_ave_err.second)) out_eff_ps<<t<<" "<<effective_mass_ave_err<<endl;
	}
      out_eff_ps<<"&"<<endl;
      
      //compute magnetization and normalize
      boot_t xi_g_num([iQ](int iconf)->double{return mag[0][iconf/compute_corr_each]*
	    (geo_Q[iconf/compute_corr_each*compute_corr_each]==iQ);});
      mag_ps[jQ]=xi_g_num/occ;
      auto mag_ps_ave_err=mag_ps[jQ].ave_err();
      if(!std::isnan(mag_ps_ave_err.second)) mag_ps_out<<iQ<<" "<<mag_ps_ave_err<<endl;
      if(iQ>=0 && kQ>=0)
	{
	  auto mag_ps_symm_ave_err=(0.5*(mag_ps[jQ]+mag_ps[kQ])).ave_err();
	  if(!std::isnan(mag_ps_symm_ave_err.second)) mag_ps_rescaled_out<<pow((double)iQ/L,2)<<" "<<mag_ps_symm_ave_err<<endl;
	}
      
      //compute xi_g
      boot_t xi_g_den([iQ](int iconf)->double{return mag[1][iconf/compute_corr_each]*
	    (geo_Q[iconf/compute_corr_each*compute_corr_each]==iQ);});
      xi_g_ps[jQ]=compute_xi_g(xi_g_num,xi_g_den);
      auto xi_g_ps_ave_err=xi_g_ps[jQ].ave_err();
      if(!std::isnan(xi_g_ps_ave_err.second)) xi_g_ps_out<<iQ<<" "<<xi_g_ps_ave_err<<endl;
  }
  
  //plot the noise
  if(nsweep<1000000) plot_noise(Q_reno);

  //compute energy
  boot_t energy([](int iconf)->double{return ene[iconf];});
  boot_t energy_rew=get_reweighted([](int iconf)->double{return ene[iconf];})/tot_weight;
  cout<<"Energy: "<<smart_print(energy.ave_err())<<endl;
  cout<<"Energy_rew: "<<smart_print((energy_rew).ave_err())<<endl;

  //reweight in beta
  if(0)
    for(double db=0.001;db<0.01;db+=0.001)
      {
	double e=energy.ave_err().first;
	boot_t energy_rew([e,db](int iconf)->double{return ene[iconf]*exp(-2*db*L*L*N*(ene[iconf]-e));});
	boot_t rew_den([e,db](int iconf)->double{return exp(-2*db*L*L*N*(ene[iconf]-e));});
	cout<<"Energy reweighted to "<<beta+db<<": "<<smart_print((energy_rew/rew_den).ave_err())<<" "<<smart_print(rew_den.ave_err())<<endl;
      }
  
  //compute magnetization
  boot_t magnetization([](int iconf)->double{return mag[0][iconf/compute_corr_each];});
  boot_t magnetization_rew=get_reweighted([](int iconf)->double{return mag[0][iconf/compute_corr_each];},compute_corr_each)/tot_weight_corr;
  
  cout<<"Magnetization: "<<smart_print(magnetization.ave_err())<<endl;
  cout<<"Magnetization_rew: "<<smart_print(magnetization_rew.ave_err())<<endl;
  
  //compute autocorrelation time
  compute_mag_tint("plots/mag_correlation.xmg",mag[0]);  

  //compute xi_g
  boot_t xi_g_den([](int iconf)->double{int icorr=iconf/compute_corr_each;return mag[1][icorr];});
  boot_t xi_g_den_rew=get_reweighted([](int iconf)->double{int icorr=iconf/compute_corr_each;return mag[1][icorr];},compute_corr_each)/tot_weight_corr;
  boot_t xi_g=compute_xi_g(magnetization,xi_g_den);
  boot_t xi_g_rew=compute_xi_g(magnetization_rew,xi_g_den_rew);
  cout<<"Xi_g: "<<smart_print(xi_g.ave_err())<<", L/Xi_g: "<<smart_print((L/xi_g).ave_err())<<endl;
  cout<<"Xi_g_rew: "<<smart_print(xi_g_rew.ave_err())<<", L/Xi_g_rew: "<<smart_print((L/xi_g_rew).ave_err())<<endl;
  
  //compute correlation function
  vector<boot_t> correlation(L/2+1);
  for(int t=0;t<L/2+1;t++) correlation[t].populate([t](int iconf)->double{return corr[t][iconf/compute_corr_each];});
  
  //print the correlation
  ofstream out_corr("plots/correlation.xmg");
  out_corr<<"@type xydy"<<endl;
  for(int t=0;t<L/2+1;t++) out_corr<<t<<" "<<correlation[t].ave_err()<<endl;
  out_corr.close();
  
  //compute and print the effective mass
  auto eff=effective_mass(correlation);
  ofstream out_eff("plots/effective_mass.xmg");
  out_eff<<"@type xydy"<<endl;
  for(int t=0;t<L/2;t++) out_eff<<t<<" "<<eff[t].ave_err()<<endl;
  out_eff.close();  
  
  //compute the number of transitions
  int base_trans=0;
  int ntrans=0;
  int ntrue_trans=0;
  int prev_Q=geo_Q[nterm];
  int prev_true_Q=geo_Q[nterm];
  int same_Q=0;
  for(int iconf=nterm;iconf<nsweep;iconf++)
    {
      int cur_Q=geo_Q[iconf];
      
      //if we changed mark down
      if(prev_Q!=cur_Q)
	{
	  ntrans++;
	  same_Q=0;
	}
      
      //increase the time since which we are in the same Q
      same_Q++;
      
      //if we are in the same Q since 10, this is a state
      if(same_Q==10)
	{
	  if(prev_true_Q!=cur_Q) ntrue_trans++;
	  prev_true_Q=cur_Q;
	}

      //mark previous step
      prev_Q=cur_Q;
      base_trans++;
    }
  
  double freq_trans=bayes_bin_ave(ntrans,base_trans),err_trans=bayes_bin_err(ntrans,base_trans);
  double freq_true_trans=bayes_bin_ave(ntrue_trans,base_trans),err_true_trans=bayes_bin_err(ntrue_trans,base_trans);
  cout<<"Frequency of transitions: "<<freq_trans<<" "<<err_trans<<" = exp("<<log(freq_trans)<<"["<<err_trans/freq_trans<<"])"<<endl;
  cout<<"Frequency of true transitions: "<<freq_true_trans<<" "<<err_true_trans<<" = exp("<<log(freq_true_trans)<<"["<<err_true_trans/freq_true_trans<<"])"<<endl;
  
  cout<<">>Time: "<<time(0)-init_time<<endl;  
  
  return 0;
}
