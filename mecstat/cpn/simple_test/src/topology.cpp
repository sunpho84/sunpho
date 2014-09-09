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
#include "types.hpp"
#include "zeta.hpp"

#include <iostream>
#include <fstream>
#include <algorithm>

std::vector<double> chrono_topo_past_values;
std::vector<double> chrono_topo_past_weight;

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

//compute the topological potential according t
double compute_theta_pot_der(double Q)
{
  //define prefactors
  double pref_gauss=-chrono_topo_coeff/sqr(chrono_topo_width);
  double pref_bend=chrono_topo_coeff*chrono_topo_bend;
  
  //set to zero results and compute total number of contributions
  double topote_der=0;
  int nchrono=chrono_topo_past_values.size();

  //inside the barrier
  if(Q>-chrono_topo_barr && Q<+chrono_topo_barr)
    {
#pragma omp parallel for reduction(+:topote_der)
      for(int i=0;i<nchrono;i++)
	{
	  double q=chrono_topo_past_values[i];
	  double w=chrono_topo_past_weight[i];
	  double diff=Q-q,f=diff/chrono_topo_width;
	  double cont=w*(pref_gauss*diff+pref_bend*Q)*exp(-f*f/2+chrono_topo_bend*Q*Q/2);
	  topote_der+=cont;
	}
    }
  
  //compute parabolic barrier (overwrite gausian piece because it must be set to zero)
  if(Q<-chrono_topo_barr) topote_der=-chrono_topo_force_out*(-Q-chrono_topo_barr);
  if(Q>+chrono_topo_barr) topote_der=+chrono_topo_force_out*(+Q-chrono_topo_barr);
  
  return topote_der;
}
//wrapper
double compute_theta_pot_der(dcomplex *l)
{return compute_theta_pot_der(topology(l));}

//compute the topodynamical potential using past history
double compute_theta_pot(double Q)
{
  //compute parabolic barrier derivative
  double harm_potential=0;
  if(Q<-chrono_topo_barr) harm_potential=chrono_topo_force_out*sqr(-Q-chrono_topo_barr)/2;
  if(Q>+chrono_topo_barr) harm_potential=chrono_topo_force_out*sqr(+Q-chrono_topo_barr)/2;
  
  //put inside the barrier in any case
  if(Q<-chrono_topo_barr) Q=-chrono_topo_barr;
  if(Q>+chrono_topo_barr) Q=+chrono_topo_barr;
  
  //set to zero results and compute total number of contributions
  double gauss_topotential=0;
  int nchrono=chrono_topo_past_values.size();
  
  //inside the barrier
#pragma omp parallel for reduction(+:gauss_topotential)
  for(int i=0;i<nchrono;i++)
    {
      double q=chrono_topo_past_values[i];
      double w=chrono_topo_past_weight[i];
      double diff=Q-q,f=diff/chrono_topo_width;
      double cont=w*exp(-f*f/2+chrono_topo_bend*Q*Q/2);
      gauss_topotential+=cont;
    }
  
  //add correct normalization
  gauss_topotential*=chrono_topo_coeff;
  
  return gauss_topotential+harm_potential;
}

//compute the topodynamical potential using past history
double compute_theta_pot(dcomplex *l)
{return compute_theta_pot(topology(l));}

//draw the chronological topological potential
void draw_chrono_topo_potential()
{
  double Q_min=*std::min_element(chrono_topo_past_values.begin(),chrono_topo_past_values.end());
  double Q_max=*std::max_element(chrono_topo_past_values.begin(),chrono_topo_past_values.end());
  double Q_diff=Q_max-Q_min;
  int n=ceil(Q_diff/chrono_topo_width*20);
  if(n==0) n=1;
  double dQ=Q_diff/n;
  
  //compute 
  double *Qy=new double[n+1];
#pragma omp parallel for
  for(int i=0;i<=n;i++) Qy[i]=compute_theta_pot(Q_min+i*dQ);
  
  //write
  ofstream fout("topo_potential");
  for(int i=0;i<=n;i++) fout<<Q_min+i*dQ<<" "<<Qy[i]<<endl;
  fout.close();

  delete[] Qy;
}

//draw the chronological topological force
void draw_chrono_topo_force()
{
  double Q_min=*std::min_element(chrono_topo_past_values.begin(),chrono_topo_past_values.end());
  double Q_max=*std::max_element(chrono_topo_past_values.begin(),chrono_topo_past_values.end());
  double Q_diff=Q_max-Q_min;
  int n=ceil(Q_diff/chrono_topo_width*20);
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
