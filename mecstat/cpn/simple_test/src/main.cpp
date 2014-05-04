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
  int nsweep=1000000;
  for(int isweep=1;isweep<=nsweep;isweep++)
    {
      //metro_sweep();
      //for(int imicro=0;imicro<3;imicro++) micro_sweep();
      //overheat_sweep();
      hmc_update();
      double topo_sim=geometric_topology_simplified(zeta);
      //double topo=geometric_topology(zeta);
      
      copy_lambda_conf(lambda_old,lambda);
      int nsto=10;
      for(int isto=0;isto<nsto;isto++)
	{
	  double topo_num=topology(lambda_old);
	  
	  energy_file<<isweep<<" "<<isto<<" "<<energy(lambda_old,zeta)/V/NDIMS<<endl;
	  topology_file<<isweep<<" "<<isto<<" "<<
	    //topo<<" "<<
	    topo_sim<<" "<<
	    topo_num<<endl;
	  
	  stout_lambda(lambda_old,0.05,lambda_old);
	}

      //write time progress
      if(isweep%(nsweep/100)==0) cout<<isweep*100/nsweep<<"%, "<<time(0)-init_time<<" s"<<endl;
    }
  
  //finalize
  close();
  
  return 0;
}