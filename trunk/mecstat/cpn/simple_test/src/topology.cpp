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

//return the geometric definition of topology
double geometric_topology_simplified(dcomplex *z)
{
  double topo=0;
  int mu=0,nu=1;
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
  for(int n=0;n<V;n++)
    topo+=(l[n*NDIMS+mu]*l[neighup(n,mu)*NDIMS+nu]*conj(l[neighup(n,nu)*NDIMS+mu]*l[n*NDIMS+nu])).imag();
  
  return topo/(2*M_PI);
}

//compute the topological potential according t
double compute_theta_pot_der(dcomplex *l)
{
  double topote=0,pref=-chrono_topo_coeff/sqr(chrono_topo_width);
  double Q=topology(l);
  
  for(std::vector<double>::iterator it=chrono_topo_past_values.begin();it!=chrono_topo_past_values.end();it++)
    {
      double q=*it;
      double diff=Q-q,f=diff/chrono_topo_width,cont=pref*diff*exp(-f*f/2);
      topote+=cont;
      //cout<<"contribution: Q="<<Q<<", q="<<q<<", "<<cont<<endl;
    }
  
  return topote;
}

//compute the topodynamical potential using past history
double compute_theta_pot(double Q)
{
  double topotential=0;
  for(std::vector<double>::iterator it=chrono_topo_past_values.begin();it!=chrono_topo_past_values.end();it++)
    {
      double q=*it;
      double diff=Q-q,f=diff/chrono_topo_width,cont=exp(-f*f/2);
      topotential+=cont;
    }
  topotential*=chrono_topo_coeff;
    
  return topotential;
}

//compute the topodynamical potential using past history
double compute_theta_pot(dcomplex *l)
{return compute_theta_pot(topology(l));}

//draw the chronological topological potential
void draw_chrono_topo_potential()
{
  ofstream fout("topo_potential");
    
  double Q_min=*std::min_element(chrono_topo_past_values.begin(),chrono_topo_past_values.end());
  double Q_max=*std::max_element(chrono_topo_past_values.begin(),chrono_topo_past_values.end());
  double Q_diff=Q_max-Q_min;
  int n=ceil(Q_diff/chrono_topo_width*20);
  double dQ=Q_diff/n;
    
  for(double Q=Q_min;Q<=Q_max;Q+=dQ) fout<<Q<<" "<<compute_theta_pot(Q)<<endl;
  fout.close();
}

//compute the force w.r.t topological term
void compute_unsmeared_topological_derivative(dcomplex *der,dcomplex *l)
{
  //compute potential
  double pot;
  if(use_topo_pot==2) pot=compute_theta_pot_der(l);
  else pot=th_top;

  int sign[2]={-1,+1};
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
