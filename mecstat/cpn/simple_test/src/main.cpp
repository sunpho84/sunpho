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

void gauge_transform()
{
  double t=0.4;
  dcomplex et(cos(t),sin(t));
  
  lambda[site_stout*NDIMS+0]*=et;
  lambda[site_stout*NDIMS+1]*=et;
  lambda[neighdw(site_stout,0)*NDIMS+0]*=conj(et);
  lambda[neighdw(site_stout,1)*NDIMS+1]*=conj(et);
}

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
  int nsweep=100000;
  for(int isweep=1;isweep<=nsweep;isweep++)
    {
      //metro_sweep();
      //for(int imicro=0;imicro<3;imicro++) micro_sweep();
      //overheat_sweep();
      
      hmc_update();
      
      double topo_sim=geometric_topology_simplified(zeta);
      //double topo=geometric_topology(zeta);
      
      //compute topologycal charge and energy
      energy_file<<isweep<<" "<<energy(lambda,zeta)/V/NDIMS<<endl;
      stout_lambda_whole_stack(lambda_stout,stout_rho,nstout_lev,lambda);
      for(int ilev=0;ilev<=nstout_lev;ilev++)
	{
	  double topo_num=topology(lambda_stout[ilev]);
	  if(use_topo_pot==2 && ilev==nstout_lev)
	  {
	    chrono_topo_past_values.push_back(topo_num);
	    chrono_topo_past_values.push_back(-topo_num);
	    if(isweep%40==0) draw_chrono_topo_potential();
	  }
	  
	  topology_file<<isweep<<" "<<ilev<<" "<<
	    //topo<<" "<<
	  topo_sim<<" "<<
	  topo_num<<endl;
	}

      //write time progress
      //if(isweep%(nsweep/100)==0) cout<<isweep*100/nsweep<<"%, "<<time(0)-init_time<<" s"<<endl;
    }
  
  //check_stout_force();
  
  //finalize
  close();
  
  return 0;
}
