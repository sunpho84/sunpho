#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <algorithm>
#include <iostream>
#include <fstream>

#include "data.hpp"
#include "debug.hpp"
#include "geometry.hpp"
#include "random.hpp"
#include "routines.hpp"
#include "types.hpp"
#include "zeta.hpp"

using namespace std;

//system
int init_time;

///////////////////////////////////////////////////////////////////////////////////////////////

//return the result of the scalar product of two lambda
double get_lambda_scalprod(dcomplex a,dcomplex b)
{return (conj(a)*b).real();}

//return the norm of a lambda
double get_lambda_norm(dcomplex &l)
{return sqrt(norm(l));}

//reunitarize a lambda
double lambda_unitarize(dcomplex &l)
{
  double n=get_lambda_norm(l);
  l*=1/n;
  return n;
}

//return the deviation from unitarity of a lambda
inline double check_lambda_unitarity(dcomplex &l)
{return fabs(get_lambda_norm(l)-1);}

//control that each element is unitary
void check_system_unitarity(double res=TINY)
{
  //check all sites
  for(int site=0;site<V;site++)
    {
      double dev_zeta=check_zeta_unitarity(zeta(site));
      if(dev_zeta>res) CRASH("zeta norm for site %d deviates from 1 by %lg",site,dev_zeta);

      //check all links
      for(int mu=0;mu<NDIMS;mu++)
	{
	  double dev_lambda=check_lambda_unitarity(lambda(site)[mu]);
	  if(dev_lambda>res) CRASH("lambda norm for site %d mu %d deviates from 1 by %lg",site,mu,dev_lambda);
	}
    }
}

//compute the total energy or action
double energy()
{
  double res=0;
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      {
	int site_up=neighup(site,mu);
	for(int n=0;n<N;n++) res+=(conj(zeta(site_up)[n])*zeta(site)[n]*lambda(site)[mu]).real();
	}
  
  return -(2*res-2*V*NDIMS);
}
double action()
{return energy()/g;}

//return the energy/action of a single site
//NB: the total action will be half the sum of the energy of all sites!
double site_energy(int site)
{
  double res=0;
  for(int mu=0;mu<NDIMS;mu++)
    {
      int site_up=neighup(site,mu);
      for(int n=0;n<N;n++) res+=(conj(zeta(site_up)[n])*zeta(site)[n]*lambda(site)[mu]).real();
      int site_dw=neighdw(site,mu);
      for(int n=0;n<N;n++) res+=(conj(zeta(site_dw)[n])*zeta(site)[n]*conj(lambda(site_dw)[mu])).real();
    }
  
  return -(2*res-4*NDIMS);
}
inline double site_action(int site)
{return site_energy(site)/g;}

//return the energy/action of a single link
double link_energy(int site,int mu)
{
  double res=0;
  
  int site_up=neighup(site,mu);
  for(int n=0;n<N;n++) res+=(conj(zeta(site_up)[n])*zeta(site)[n]*lambda(site)[mu]).real();
  
  return -(2*res-2);
}
inline double link_action(int site,int mu)
{return link_energy(site,mu)/g;}

//compute the staple of zeta
void site_staple(dcomplex *staple,int site)
{
  for(int n=0;n<N;n++) staple[n]=0;

  for(int mu=0;mu<NDIMS;mu++)
    {
      int site_up=neighup(site,mu);
      for(int n=0;n<N;n++) staple[n]+=(double)2.0*zeta(site_up)[n]*conj(lambda(site)[mu]);
      int site_dw=neighdw(site,mu);
      for(int n=0;n<N;n++) staple[n]+=(double)2.0*zeta(site_dw)[n]*lambda(site_dw)[mu];
    }
}
inline double site_staple_energy(int site,dcomplex *staple)
{
  double res=0;
  for(int n=0;n<N;n++) res+=(staple[n]*conj(zeta(site)[n])).real();
  return -(res-4*NDIMS);
}
inline double site_staple_action(int site,dcomplex *staple)
{return site_staple_energy(site,staple)/g;}

//compute the staple of lambda
void link_staple(dcomplex &staple,int site,int mu)
{
  staple=0;
  int site_up=neighup(site,mu);
  for(int n=0;n<N;n++) staple+=(double)2.0*conj(zeta(site)[n])*zeta(site_up)[n];
}
inline double link_staple_energy(int site,int mu)
{
  dcomplex staple;
  link_staple(staple,site,mu);
  return -((staple*conj(lambda(site)[mu])).real()-2);
}
inline double link_staple_action(int site,int mu)
{return link_staple_energy(site,mu)/g;}

//return inner product of zeta
void get_P(dcomplex *P,dcomplex *z)
{
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      P[i*N+j]=conj(z[i])*z[j];
}

//return the angle of scalar product between two zetas
double arg_an(dcomplex *a,dcomplex *b)
{
  dcomplex res=0;
  for(int n=0;n<N;n++) res+=conj(a[n])*b[n];
  //cout<<log(res)<<" "<<log(conj(res))<<endl;
  return arg(res);
}
dcomplex sc(dcomplex *a,dcomplex *b)
{
  dcomplex res=0;
  for(int n=0;n<N;n++) res+=conj(a[n])*b[n];
  return res;
}

//return the geometric definition of topology
double geometric_topology_simplified()
{
  double topo=0;
  for(int n=0;n<V;n++)
    {
      int mu=0,nu=1;
      int nmu=neighup(n,mu);
      int nnu=neighup(n,nu);
      int nmu_nu=neighup(nmu,nu);
      
      topo+=arg(sc(zeta(nmu_nu),zeta(n))*sc(zeta(nmu),zeta(nmu_nu))*sc(zeta(n),zeta(nmu)))+
	arg(sc(zeta(nnu),zeta(n))*sc(zeta(nmu_nu),zeta(nnu))*sc(zeta(n),zeta(nmu_nu)));
    }
  
  return topo/(2*M_PI);
}
double geometric_topology()
{
  int mu=0,nu=1;
  double topo=0;
  for(int n=0;n<V;n++)
    {
      dcomplex P1[N*N],P2[N*N],P3[N*N];
      get_P(P1,zeta(n));
      get_P(P3,zeta(neighup(neighup(n,mu),nu)));

      dcomplex c;
      
      c=0;
      get_P(P2,zeta(neighup(n,mu)));
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  for(int k=0;k<N;k++)
	    c+=P3[i*N+j]*P2[j*N+k]*P1[k*N+i];
      topo+=log(c).imag();
      
      c=0;
      get_P(P2,zeta(neighup(n,nu)));
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  for(int k=0;k<N;k++)
	    c+=P2[i*N+j]*P3[j*N+k]*P1[k*N+i];
      topo+=log(c).imag();
    }
  
  return topo/(2*M_PI);
}

double topology()
{
  int mu=0,nu=1;
  double topo=0;
  for(int n=0;n<V;n++)
    topo+=(lambda(n)[mu]*lambda(neighup(n,mu))[nu]*conj(lambda(neighup(n,nu))[mu]*lambda(n)[nu])).imag();
  
  return topo;
}

//initialize the code
void init(int cond,int seed)
{
  init_time=time(0);
  
#ifdef GOOD_GENERATOR
  //init the random generators
  rd=new random_device();
  gen=new mt19937_64((*rd)());
  gen->seed(seed);
  dis=new uniform_real_distribution<double>;
#else
  const int im1=2147483563,ia1=40014;
  const int iq1=53668,ir1=12211;
  int j,k;
  
  //initialization
  gen.idum=seed;
  gen.idum=std::max(gen.idum+1,1);
  gen.idum2=gen.idum;
  for(j=RAN2_NTAB+7;j>=0;j--)
    {
      k=gen.idum/iq1;
      gen.idum=ia1*(gen.idum-k*iq1)-k*ir1;
      if(gen.idum<0) gen.idum+=im1;
      if(j<RAN2_NTAB) gen.iv[j]=gen.idum;
    }
  gen.iy=gen.iv[0];
#endif
    
  //geometry
  V=1;
  for(int mu=0;mu<NDIMS;mu++) V*=L;
  cout<<"Volume: "<<V<<endl;
  neigh_data=new int[V*NDIMS*2];
  
  //loop over sites
  for(int site=0;site<V;site++)
    {
      //get the original coordinates
      coords c;
      coords_of_site(c,site);
      
      //loop over directions
      for(int mu=0;mu<NDIMS;mu++)
	{
	  //save original
	  int o=c[mu];
	  
	  //backward
	  c[mu]=(o+L-1)%L;
	  neighdw(site,mu)=site_of_coords(c);
	  
	  //forward
	  c[mu]=(o+L+1)%L;
	  neighup(site,mu)=site_of_coords(c);
	  
	  //restore original
	  c[mu]=o;
	}
    }
  
  //Zeta and Lambda
  zeta_data=new dcomplex[N*V];
  lambda_data=new dcomplex[V*NDIMS];

  //set the system to hot state
  init_system_to(cond);
}

//update zeta with metropolis
void metro_update_site(int site)
{
  //copy zeta
  dcomplex ori[N];
  for(int n=0;n<N;n++) ori[n]=zeta(site)[n];
  
  //change of action
  double ori_ac=site_action(site);
  set_ON_to_rnd(zeta(site));
  double fin_ac=site_action(site);
  
  //accept or not?
  double diff_ac=fin_ac-ori_ac;
  double t=exp(-diff_ac);
  double p=get_unif_double(1);
  if(p>t) for(int n=0;n<N;n++) zeta(site)[n]=ori[n];
}

//update lambda with metropolis
void metro_update_link(int site,int mu)
{
  //copy lambdaa
  dcomplex ori=lambda(site)[mu];
  
  //change of action
  double ori_ac=link_action(site,mu);
  set_U1_to_rnd(lambda(site)[mu]);
  double fin_ac=link_action(site,mu);
  
  //accept or not?
  double diff_ac=fin_ac-ori_ac;
  double t=exp(-diff_ac);
  double p=get_unif_double(1);
  if(p>t) lambda(site)[mu]=ori;
}

//update a site using microcanonical
void micro_update_site(int site)
{
  //get the staple, its norm and energy
  dcomplex staple[N];
  site_staple(staple,site);
  double staple_norm=get_zeta_norm(staple);
  double staple_energy=get_zeta_scalprod(zeta(site),staple);
  
  //extract site
  for(int n=0;n<N;n++) zeta(site)[n]=2*staple_energy/sqr(staple_norm)*staple[n]-zeta(site)[n];

  //reunitarize
  zeta_unitarize(zeta(site));
}

//update a site using overrelaxion/heatbath
void overheat_update_site(int site)
{
  //get the staple, its norm and energy
  dcomplex staple[N];
  site_staple(staple,site);
  double staple_norm=get_zeta_norm(staple);
  double staple_energy=get_zeta_scalprod(zeta(site),staple);
  
  //compute theta in the simple way
  double ctheta_old=staple_energy/staple_norm;
  double theta_old=acos(ctheta_old);
  
  //if theta is too small we switch to alternative method
  if(fabs(theta_old)<1.e-4) 
    {
      dcomplex diff[N];
      for(int n=0;n<N;n++) diff[n]=zeta(site)[n]-staple[n]/staple_norm;
      theta_old=asin(get_zeta_scalprod(diff,diff));
    }
  
  //it theta is too small the algorithm is undefined
  if(fabs(theta_old)>5.e-8)
    {
      //extract theta and compute its cos
      double theta_new=get_theta(beta*N*staple_norm,N);
      double ctheta_new=cos(theta_new);
      
      //extract remaining components
      double a=ctheta_new/staple_norm,b=ctheta_old/staple_norm,c=sin(theta_new)/sin(theta_old);
      for(int n=0;n<N;n++) zeta(site)[n]=a*staple[n]-(zeta(site)[n]-b*staple[n])*c;
      
      //reunitarize
      zeta_unitarize(zeta(site));
    }
  else cout<<"skipping site "<<site<<": "<<theta_old<<endl;
}

//update a link using microcanonical
void micro_update_link(int site,int mu)
{
  //get the staple
  dcomplex staple;
  link_staple(staple,site,mu);

  //compute the staple norm and energy
  double staple_norm=sqrt(norm(staple));
  double staple_energy=get_lambda_scalprod(lambda(site)[mu],staple);
  
  //extract link
  lambda(site)[mu]=2*staple_energy/sqr(staple_norm)*staple-lambda(site)[mu];

  //reunitarize
  lambda_unitarize(lambda(site)[mu]);
}

//update a link using overrelaxion/heatbath
void overheat_update_link(int site,int mu)
{
  //get the staple
  dcomplex staple;
  link_staple(staple,site,mu);

  //compute the staple norm and energy
  double staple_norm=sqrt(norm(staple));
  double staple_energy=get_lambda_scalprod(lambda(site)[mu],staple);
  
  //compute theta in the simple way
  double ctheta_old=staple_energy/staple_norm;
  double theta_old=acos(ctheta_old);
  
  //if theta is too small we switch to alternative method
  if(fabs(theta_old)<1.e-4) 
    {
      dcomplex diff=lambda(site)[mu]-staple/staple_norm;
      theta_old=asin(get_lambda_scalprod(diff,diff));
    }
  
  //it theta is too small the algorithm is undefined
  if(fabs(theta_old)>5.e-8)
    {
      //extract theta and compute its cos
      double theta_new=get_theta_1(beta*N*staple_norm);
      double ctheta_new=cos(theta_new);
      
      //extract remaining components  
      double a=ctheta_new/staple_norm,b=ctheta_old/staple_norm,c=sin(theta_new)/sin(theta_old);
      lambda(site)[mu]=a*staple-(lambda(site)[mu]-b*staple)*c;
      
      //reunitarize
      lambda_unitarize(lambda(site)[mu]);
    }
}

//sweep all the lattice
void metro_sweep()
{
  //loop over sites
  for(int site=0;site<V;site++)
    {
      metro_update_site(site);
      for(int mu=0;mu<NDIMS;mu++) metro_update_link(site,mu);
    }
}

//sweep all the lattice with microcanonical
void micro_sweep()
{
  //loop over sites
  for(int site=0;site<V;site++)
    {
      micro_update_site(site);
      for(int mu=0;mu<NDIMS;mu++) micro_update_link(site,mu);
    }
}

//sweep all the lattice with overrelaxation/heatbath
void overheat_sweep()
{
  //loop over sites
  for(int site=0;site<V;site++)
    {
      overheat_update_site(site);
      for(int mu=0;mu<NDIMS;mu++) overheat_update_link(site,mu);
    }
}

//perform a hybid monte carlo update
void hmc_update()
{
  //allocate momenta
  dcomplex *pi=new dcomplex[V*N];
  double *omega=new double[V*NDIMS];
  
  //draw momenta
  for(int site=0;site<V;site++)
    {
      
    }
  
  delete[] pi;
  delete[] omega;
}

//close the code
void close()
{
  delete[] lambda_data;
  delete[] zeta_data;
  
  delete[] neigh_data;
  
#ifdef GOOD_GENERATOR
  delete dis;
  delete gen;
  delete rd;
#endif
}

int main()
{
  //initialize
  init(HOT,100);

  ofstream energy_file("energy");
  energy_file.precision(16);
  ofstream topology_file("topology");
  
  //sweep with overheat-micro
  int nsweep=1000000;
  for(int isweep=1;isweep<=nsweep;isweep++)
    {
      //metro_sweep();
      for(int imicro=0;imicro<3;imicro++) micro_sweep();
      overheat_sweep();
      
      double topo_sim=geometric_topology_simplified();
      double topo_num=topology();
      
      energy_file<<energy()/V/NDIMS<<endl;
      topology_file<<topo_sim<<" "<<topo_num<<endl;
      
      //write time progress
      if(isweep%(nsweep/100)==0) cout<<isweep*100/nsweep<<"%, "<<time(0)-init_time<<" s"<<endl;
    }
  
  //finalize
  close();
  
  return 0;
}
