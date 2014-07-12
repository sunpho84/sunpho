#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>
#include <fstream>

#include "action.hpp"
#include "close.hpp"
#include "corr.hpp"
#include "data.hpp"
#include "debug.hpp"
#include "geometry.hpp"
#include "hmc.hpp"
#include "init.hpp"
#include "lambda.hpp"
#include "metro.hpp"
#include "micro.hpp"
#include "overheat.hpp"
#include "random.hpp"
#include "routines.hpp"
#include "staples.hpp"
#include "stout.hpp"
#include "topology.hpp"
#include "types.hpp"
#include "zeta.hpp"

#include <omp.h>

using namespace std;

int mu_stout=1;
int site_stout=4;

void gauge_transform(int s)
{
  double t=rand()/((double)RAND_MAX)*2*M_PI;
  
  dcomplex et(cos(t),sin(t));
  
  for(int n=0;n<N;n++) zeta[s*N+n]*=conj(et);
  lambda[s*NDIMS+0]*=et;
  lambda[s*NDIMS+1]*=et;
  lambda[neighdw(s,0)*NDIMS+0]*=conj(et);
  lambda[neighdw(s,1)*NDIMS+1]*=conj(et);
}
void gauge_transform()
{for(int s=0;s<V;s++) gauge_transform(s);}

//compute functional with stout
double compute_functional()
{
  stout_lambda_whole_stack(lambda_stout,stout_rho,nstout_lev,lambda);
  return topology(lambda_stout[nstout_lev]);
}

//compute the force
void check_stout_force()
{
  int site_check=site_stout;
  int mu_check=mu_stout;
  
  cout<<lambda[site_check*NDIMS+mu_check]<<endl;  
  //gauge_transform();
  cout<<lambda[site_check*NDIMS+mu_check]<<endl;
  compute_topological_force(fomega,stout_rho,nstout_lev,lambda);
  
  double pre_act=compute_functional();
  
  //change
  dcomplex pre=lambda[site_check*NDIMS+mu_check];
  double eps=1.e-6;
  lambda[site_check*NDIMS+mu_check]*=dcomplex(cos(eps),sin(eps));
  
  double post_act=compute_functional();
  lambda[site_check*NDIMS+mu_check]=pre;
      
  double f=-(post_act-pre_act)/eps;
  std::cout<<"mu: "<<mu_check<<" fnu: -("<<post_act<<"-"<<pre_act<<")/"<<eps<<"="<<-(post_act-pre_act)<<"/"<<eps<<"="<<
    f<<" fan: "<<fomega[site_check*NDIMS+mu_check]<<", ratio-1: "<<fomega[site_check*NDIMS+mu_check]/f-1<<std::endl;
}

//crash promptin error message
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

//read an element from input file
template <class T> void read(T &out,ifstream &in,string is)
{
  string s;
  if(!(in>>s)) crash("impossible to read expect string \"%s\"",is.c_str());
  if(s!=is) crash("obtained %s while reading %s",s.c_str(),is.c_str());
  if(!(in>>out)) crash("reading data");
}

int main()
{
  //read parameters
  ifstream input("input");
  if(!input.good()) crash("opening input");
  read(N,input,"N");
  read(L,input,"L");
  read(beta,input,"Beta");
  g=1/(N*beta);
  int seed;
  read(seed,input,"Seed");
  int nsweep,nterm,nmicro;
  read(nsweep,input,"NSweep");
  read(nterm,input,"NTerm");
  read(nmicro,input,"NMicro");
  read(nstout_lev,input,"NStoutLev");
  read(stout_rho,input,"StoutRho");
  
  //initialize
  init(HOT,seed);
  
  //thermalization
  for(int isweep=0;isweep<nterm;isweep++)
  {
    for(int imicro=0;imicro<nmicro;imicro++) micro_sweep();
    overheat_sweep();
  }
  
  //print number of threads
#pragma omp parallel
  {
#pragma omp single
    cout<<omp_get_num_threads()<<" threads"<<endl;
  }
  
  //open output files
  ofstream energy_file("energy");
  ofstream topology_file("topology");
  ofstream corr_file("corr");
  ofstream corrd_file("corrd");
  ofstream mag_file("mag");
  energy_file.precision(12);
  topology_file.precision(12);
  corr_file.precision(12);
  corrd_file.precision(12);
  mag_file.precision(12);
  
  //sweep with overheat-micro
  int init_time=time(0);
  for(int isweep=1;isweep<=nsweep;isweep++)
    {
      //metro_sweep();
      for(int imicro=0;imicro<nmicro;imicro++) micro_sweep();
      overheat_sweep();
      
      //hmc_update();
      
      //compute energy and topological charge
      energy_file<<isweep<<" "<<energy(zeta,lambda)/V/NDIMS<<endl;
      double topo_sim=geometric_topology_simplified(zeta);
      
      //compute topologycal charge
      stout_lambda_whole_stack(lambda_stout,stout_rho,nstout_lev,lambda);
      for(int ilev=0;ilev<=nstout_lev;ilev++)
	{
	  double topo_num=topology(lambda_stout[ilev]);
	  if(use_topo_pot==2 && ilev==nstout_lev)
	  {
	    chrono_topo_past_values.push_back(+topo_num);
	    chrono_topo_past_values.push_back(-topo_num);
	    if(isweep%500==0) draw_chrono_topo_potential();
	  }
	  
	  topology_file<<isweep<<" "<<ilev<<" "<<
	    //topo<<" "<<
	    topo_sim<<" "<<
	    topo_num<<endl;
	}
      
      //compute the correlation function
      double mag0,mag1,corr[L],corrd[L];
      compute_corr(mag0,mag1,corr,corrd,zeta);
      for(int i=0;i<=L/2;i++) corr_file<<isweep<<" "<<i<<" "<<corr[i]<<endl;
      for(int i=0;i<=L/2;i++) corrd_file<<isweep<<" "<<i/sqrt(2)<<" "<<corrd[i]<<endl;
      mag_file<<isweep<<" "<<mag0<<" "<<mag1<<endl;
    }
  
  //write lasted time
  cout<<time(0)-init_time<<" s"<<endl;
  
  //finalize
  close();
  
  return 0;
}
