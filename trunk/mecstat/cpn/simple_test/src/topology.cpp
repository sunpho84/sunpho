#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "action.hpp"
#include "data.hpp"
#include "geometry.hpp"
#include "parameters.hpp"
#include "routines.hpp"
#include "staples.hpp"
#include "stout.hpp"
#include "tools.hpp"
#include "types.hpp"
#include "zeta.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>

int ngrid;
vector<double> topo_grid;

//return the geometric definition of topology
double geometric_topology_simplified(dcomplex *z)
{
  double topo=0;
  int mu=0,nu=1;
#pragma omp parallel for reduction(+:topo)
  for(int n=0;n<V;n++)
    {
      int nmu=neighup(n,mu);
      int nnu=neighup(n,nu);
      int nmu_nu=neighup(nmu,nu);
      
      topo+=
	arg(get_zeta_compl_scalprod(z+nmu_nu*N,z+n*N)*
	    get_zeta_compl_scalprod(z+nmu*N,z+nmu_nu*N)*
	    get_zeta_compl_scalprod(z+n*N,z+nmu*N))+
	arg(get_zeta_compl_scalprod(z+nnu*N,z+n*N)*
	    get_zeta_compl_scalprod(z+nmu_nu*N,z+nnu*N)*
	    get_zeta_compl_scalprod(z+n*N,z+nmu_nu*N));
    }
  
  return topo/(2*M_PI);
}

//alternative edition
double geometric_topology(dcomplex *z)
{
  int mu=0,nu=1;
  double topo=0;
#pragma omp parallel for reduction(+:topo)
  for(int n=0;n<V;n++)
    {
      dcomplex P1[N*N],P2[N*N],P3[N*N];
      get_zeta_P(P1,z+n*N);
      get_zeta_P(P3,z+neighup(neighup(n,mu),nu)*N);

      dcomplex c;
      
      c=0;
      get_zeta_P(P2,z+neighup(n,mu)*N);
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  for(int k=0;k<N;k++)
	    c+=P3[i*N+j]*P2[j*N+k]*P1[k*N+i];
      topo+=log(c).imag();
      
      c=0;
      get_zeta_P(P2,z+neighup(n,nu)*N);
      for(int i=0;i<N;i++)
	for(int j=0;j<N;j++)
	  for(int k=0;k<N;k++)
	    c+=P2[i*N+j]*P3[j*N+k]*P1[k*N+i];
      topo+=log(c).imag();
    }
  
  return topo/(2*M_PI);
}

//gauge version
double topology(dcomplex *l)
{
  int mu=0,nu=1;
  double topo=0;
#pragma omp parallel for reduction(+:topo)
  for(int n=0;n<V;n++)
    topo+=(l[n*NDIMS+mu]*l[neighup(n,mu)*NDIMS+nu]*conj(l[neighup(n,nu)*NDIMS+mu]*l[n*NDIMS+nu])).imag();
  
  return topo/(2*M_PI);
}

//update the history-dependent potential
void update_chrono_potential(double Q)
{
  int igrid=floor(Q/chrono_topo_width)+ngrid/2;
  double alpha=Q/chrono_topo_width;
  alpha=alpha-floor(alpha);
  if(igrid>=0 && igrid<=ngrid) topo_grid[igrid]+=(1-alpha)*chrono_topo_coeff;
  if(igrid+1>=0 && igrid+1<=ngrid) topo_grid[igrid+1]+=alpha*chrono_topo_coeff;
}

//compute the derivative of the topological potential
double compute_theta_pot_der(double Q)
{
  //take igrid
  int igrid=floor((Q+chrono_topo_barr)/chrono_topo_width);
  
  //inside the barriers
  if(igrid>=0 && igrid<ngrid)
    return (topo_grid[igrid+1]-topo_grid[igrid])/chrono_topo_width;
  else
    if(igrid<0)
      return -chrono_topo_force_out*(-Q-chrono_topo_barr);
    else
      return +chrono_topo_force_out*(+Q-chrono_topo_barr);
}
//wrapper
double compute_theta_pot_der(dcomplex *l)
{return compute_theta_pot_der(topology(l));}

//compute the topodynamical potential using past history
double compute_theta_pot(double Q)
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
 
//compute the topodynamical potential using past history
double compute_theta_pot(dcomplex *l)
{return compute_theta_pot(topology(l));}

//draw the chronological topological potential
void draw_chrono_topo_potential()
{
  //write
  ofstream fout("topo_potential");
  fout.precision(16);
  for(int i=0;i<=ngrid;i++) fout<<-chrono_topo_barr+i*chrono_topo_width<<" "<<topo_grid[i]<<endl;
  fout.close();
}
void load_chrono_topo_potential()
{
  //read
  ifstream fin("topo_potential");
  if(!fin.good()) crash("opening \"topo_potential\"");
  for(int igrid=0;igrid<=ngrid;igrid++)
    {
      double xread;
      fin>>xread>>topo_grid[igrid];
      if(!fin.good()) crash("reading line %d of \"topo_potential\"",igrid);
      int jgrid=floor((xread+chrono_topo_barr+chrono_topo_width/2)/chrono_topo_width);
      if(igrid!=jgrid) crash("found %d (%lg) when expecting %d",jgrid,xread,igrid);
    }
  fin.close();
}

//draw the chronological topological force
void draw_chrono_topo_force()
{
  double Q_min=-chrono_topo_barr*1.1;
  double Q_max=+chrono_topo_barr*1.1;
  double Q_diff=Q_max-Q_min;
  int n=ceil(Q_diff/chrono_topo_width*10);
  if(n==0) n=1;
  double dQ=Q_diff/n;
  
  //compute 
  double *Qy=new double[n+1];
  double *Qz=new double[n+1];
#pragma omp parallel for
  for(int i=0;i<=n;i++)
    {
      Qy[i]=compute_theta_pot_der(Q_min+i*dQ);
      Qz[i]=(compute_theta_pot(Q_min+i*dQ+dQ/10)-compute_theta_pot(Q_min+i*dQ-dQ/10))/(dQ/5);
    }
  
  //write
  ofstream fout("topo_force");
  for(int i=0;i<=n;i++) fout<<Q_min+i*dQ<<" "<<Qy[i]<<endl;
  fout<<"&"<<endl;
  for(int i=0;i<=n;i++) fout<<Q_min+i*dQ<<" "<<Qz[i]<<endl;
  fout.close();
  
  delete[] Qy;
  delete[] Qz;
}

//compute the force w.r.t topological term
void compute_unsmeared_topological_derivative(dcomplex *der,dcomplex *l)
{
  //compute potential
  double pot;
  if(use_topo_pot==2) pot=compute_theta_pot_der(l);
  else pot=th_top;

  int sign[2]={-1,+1};
#pragma omp parallel for
  for(int s=0;s<V;s++)
    for(int mu=0;mu<2;mu++)
      {
	topo_staple(der[s*NDIMS+mu],l,s,mu);
	der[s*NDIMS+mu]*=sign[mu]*pot/(2*M_PI);
      }
}

//finish computation of the force
void finish_topological_force(double *f,dcomplex *st,dcomplex *l)
{
#pragma omp parallel for
  for(int s=0;s<V;s++)
    for(int mu=0;mu<NDIMS;mu++)
      f[s*NDIMS+mu]=(l[s*NDIMS+mu]*st[s*NDIMS+mu]).real();
}

//compute and finish
void compute_unsmeared_topological_force(double *f,dcomplex *l)
{
  compute_unsmeared_topological_derivative(topo_staples_data,l);
  finish_topological_force(f,topo_staples_data,l);
  
#ifdef DEBUG_HMC
  int site=0;
  double eps=1.e-6;
  double pre_act=topo_action(lambda);
  for(int mu=0;mu<NDIMS;mu++)
    {
      dcomplex pre=lambda[site*NDIMS+mu];
      lambda[site*NDIMS+mu]*=dcomplex(cos(eps),sin(eps));
      
      double post_act=topo_action(lambda);
      lambda[site*NDIMS+mu]=pre;
      
      double f=-(post_act-pre_act)/eps;
      std::cout<<"mu: "<<mu<<" fnu: -("<<post_act<<"-"<<pre_act<<")/"<<eps<<"="<<-(post_act-pre_act)<<"/"<<eps<<"="<<
	f<<" fan: "<<fomega[site*NDIMS+mu]<<" "<<fomega[site*NDIMS+mu]/f<<std::endl;
    }
#endif
}

//compute the force for topological charge, smeared or not
void compute_topological_force(double *f,double rho,int nlev,dcomplex *l)
{
  if(nlev==0) compute_unsmeared_topological_force(f,l);
  else
    {
      stout_lambda_whole_stack(lambda_stout,rho,nlev,l);
      compute_unsmeared_topological_derivative(topo_staples_data,lambda_stout[nlev]);
      stout_remap_force(topo_staples_data,rho,nlev,lambda_stout);
      finish_topological_force(f,topo_staples_data,l);
    }
}
