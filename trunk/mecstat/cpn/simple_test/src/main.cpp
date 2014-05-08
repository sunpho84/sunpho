#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>
#include <fstream>

#include "action.hpp"
#include "close.hpp"
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

using namespace std;

int mu_stout=1;
int site_stout=4;
double stout_coeff=0.05;

void gauge_transform()
{
  double t=0.4;
  dcomplex et(cos(t),sin(t));
  
  lambda[site_stout*NDIMS+0]*=et;
  lambda[site_stout*NDIMS+1]*=et;
  lambda[neighdw(site_stout,0)*NDIMS+0]*=conj(et);
  lambda[neighdw(site_stout,1)*NDIMS+1]*=conj(et);
}

void fake_stout(dcomplex *out,dcomplex *in)
{
  //produce the stouted conf
  stout_lambda(out,stout_coeff,in);
  
  //copy everyting but link 0 dir 0 back
  for(int n=0;n<V;n++)
    for(int mu=0;mu<NDIMS;mu++)
      if(0 && (n!=site_stout||mu!=mu_stout))
	out[n*NDIMS+mu]=in[n*NDIMS+mu];
}

//compute functional with stout
double compute_functional()
{
  fake_stout(lambda_old,lambda);
  return topology(lambda_old);
}

//compute derivative
void compute_link_derivative(dcomplex *staple,dcomplex *l)
{
  int sign[2]={-1,+1};
  for(int s=0;s<V;s++)
    for(int mu=0;mu<NDIMS;mu++)
      {
        topo_staple(staple[s*NDIMS+mu],l,s,mu);
        staple[s*NDIMS+mu]*=sign[mu]/(2*M_PI);
      }
}

void stout_remap_force(dcomplex *&staple,dcomplex *&staple_supp,dcomplex *l)
{
  for(int site=0;site<V;site++)
    for(int mu=0;mu<NDIMS;mu++)
      {
	int nu=!mu;
	
	int A=site,B=neighup(A,nu);//,C=neighup(B,mu);
	int D=neighup(A,mu),F=neighdw(A,nu),E=neighup(F,mu);
	
	double Q=stout_coeff*(conj(l[A*NDIMS+mu])*(l[A*NDIMS+nu]*l[B*NDIMS+mu]*conj(l[D*NDIMS+nu])+
						   conj(l[F*NDIMS+nu])*l[F*NDIMS+mu]*l[E*NDIMS+nu]
						   )).imag();
	
	dcomplex FW=lambda[D*NDIMS+nu]*conj(lambda[B*NDIMS+mu]*lambda[A*NDIMS+nu]);
	dcomplex BW=conj(lambda[E*NDIMS+nu]*lambda[F*NDIMS+mu])*lambda[F*NDIMS+nu];
	
	double LA_mu=(lambda_old[A*NDIMS+mu]*staple[A*NDIMS+mu]).real();
	double LA_nu=(lambda_old[A*NDIMS+nu]*staple[A*NDIMS+nu]).real();
	double LB_mu=(lambda_old[B*NDIMS+mu]*staple[B*NDIMS+mu]).real();
	double LD_nu=(lambda_old[D*NDIMS+nu]*staple[D*NDIMS+nu]).real();
	double LF_nu=(lambda_old[F*NDIMS+nu]*staple[F*NDIMS+nu]).real();
	double LF_mu=(lambda_old[F*NDIMS+mu]*staple[F*NDIMS+mu]).real();
	double LE_nu=(lambda_old[E*NDIMS+nu]*staple[E*NDIMS+nu]).real();
	
	dcomplex eiQ(cos(Q),sin(Q));
	
	staple_supp[A*NDIMS+mu]=
	  eiQ*staple[A*NDIMS+mu]+stout_coeff*
	       (-FW*(LA_mu-(+LA_nu+LB_mu-LD_nu))
		-BW*(LA_mu-(-LF_nu+LF_mu+LE_nu)));
      }
  
  swap(staple,staple_supp);
}

//finish computation of the force
void finish_force(double *f,dcomplex *staple)
{
  for(int s=0;s<V;s++)
    for(int mu=0;mu<NDIMS;mu++)
      f[s*NDIMS+mu]=(lambda[s*NDIMS+mu]*staple[s*NDIMS+mu]).real();
}

//compute the whole force
void compute_force()
{
  dcomplex *staple=new dcomplex[V*NDIMS];
  dcomplex *staple_supp=new dcomplex[V*NDIMS];
  
  fake_stout(lambda_old,lambda);
  compute_link_derivative(staple,lambda_old);
  stout_remap_force(staple,staple_supp,lambda);
  finish_force(fomega,staple);
  
  delete[] staple_supp;
  delete[] staple;
}

//compute the force
void check_stout_force()
{
  int site_check=site_stout;
  int mu_check=mu_stout;
  
  cout<<lambda[site_check*NDIMS+mu_check]<<endl;  
  gauge_transform();
  cout<<lambda[site_check*NDIMS+mu_check]<<endl;
  compute_force();
  
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

int main()
{
  //initialize
  init(HOT,100);
  
  for(int itraj=0;itraj<10;itraj++)
  {
    for(int imicro=0;imicro<3;imicro++) micro_sweep();
    overheat_sweep();
  }

  //CRASH("");
  ofstream energy_file("energy");
  energy_file.precision(16);
  ofstream topology_file("topology");
  
  //sweep with overheat-micro
  int nsweep=1;//0000;
  for(int isweep=1;isweep<=nsweep;isweep++)
    {
      //metro_sweep();
      //for(int imicro=0;imicro<3;imicro++) micro_sweep();
      //overheat_sweep();
      
      hmc_update();
      
      double topo_sim=geometric_topology_simplified(zeta);
      //double topo=geometric_topology(zeta);
      
      //compute topologycal charge and energy
      copy_lambda_conf(lambda_old,lambda);
      int nsto=10;
      for(int isto=0;isto<nsto;isto++)
	{
	  double topo_num=topology(lambda_old);
	  if(use_topo_pot==2 && isto==0)
	    {
	      chrono_topo_past_values.push_back(topo_num);
	      chrono_topo_past_values.push_back(-topo_num);
	      if(isweep%40==0) draw_chrono_topo_potential();
	    }
	  
	  energy_file<<isweep<<" "<<isto<<" "<<energy(lambda_old,zeta)/V/NDIMS<<endl;
	  topology_file<<isweep<<" "<<isto<<" "<<
	    //topo<<" "<<
	    topo_sim<<" "<<
	    topo_num<<endl;
	  
	  stout_lambda(lambda_old,0.05,lambda_old);
	}

      //write time progress
      //if(isweep%(nsweep/100)==0) cout<<isweep*100/nsweep<<"%, "<<time(0)-init_time<<" s"<<endl;
    }
  
  check_stout_force();
  
  //finalize
  close();
  
  return 0;
}
