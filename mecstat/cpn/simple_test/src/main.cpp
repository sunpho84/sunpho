#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <fstream>
#include <iostream>

#include "cpn.hpp"

using namespace std;

int main()
{
  //read input and initialize
  read_pars_t read_pars;
  read_input(read_pars,"input");
  int base_isweep;
  init(base_isweep,read_pars);
  
  //open output files
  ios::openmode mode=ios::out;
  if(read_pars.start_cond==LOAD) mode|=ios::app;
  ofstream energy_file("energy",mode);
  ofstream weight_file("topo_weight",mode);
  ofstream topology_file("topology",mode);
  ofstream corr_file("corr",mode);
  ofstream corrd_file("corrd",mode);
  ofstream mom2_file("mom2",mode);
  ofstream mag_file("mag",mode);
  energy_file.precision(PREC);
  weight_file.precision(PREC);
  topology_file.precision(PREC);
  corr_file.precision(PREC);
  corrd_file.precision(PREC);
  mag_file.precision(PREC);
  
  //sweep with overheat-micro
  timing_t tot_time,sweep_time,energy_time,geo_topo_time,topo_time,corr_time;
  
  tot_time.start();
  int isweep;
  for(isweep=base_isweep;isweep<base_isweep+read_pars.nsweep;isweep++)
    {
      //compute geometrical topological charge
      geo_topo_time.start();
      double topo_sim=geometric_topology_simplified(zeta);
      geo_topo_time.stop();
      
      //compute energy
      energy_time.start();
      energy_file<<isweep<<" "<<energy(zeta,lambda)/V/NDIMS<<endl;
      energy_time.stop();
      
      //compute topologycal charge
      topo_time.start();
      stout_lambda_whole_stack(lambda_stout,stout_rho,nstout_lev,lambda);
      for(int ilev=0;ilev<=nstout_lev;ilev++)
	{
	  double topo_num=topology(lambda_stout[ilev]);
	  if(use_topo_pot==2 && ilev==nstout_lev && isweep>=chrono_topo_after)
	    {
	      double w=exp(-chrono_topo_well_tempering*compute_theta_pot(+topo_num));
	      chrono_topo_past_values.push_back(+topo_num);
	      chrono_topo_past_weight.push_back(w);
	      weight_file<<w<<endl;
	      if(isweep%DRAW_EACH==0)
		{
		  draw_chrono_topo_potential(true);
		  draw_chrono_topo_force();
		}
	    }
	  
	  topology_file<<isweep<<" "<<ilev<<" "<<
	    topo_sim<<" "<<
	    topo_num<<endl;
	}
      topo_time.stop();
      
      //compute the correlation function
      if(isweep%compute_corr_each==0)
	{
	  corr_time.start();
	  double mag0,mag1,mom2,corr[L],corrd[L];
	  compute_corr(mag0,mag1,mom2,corr,corrd,zeta);
	  for(int i=0;i<=L/2;i++) corr_file<<isweep<<" "<<i<<" "<<corr[i]<<endl;
	  for(int i=0;i<=L/2;i++) corrd_file<<isweep<<" "<<i/sqrt(2)<<" "<<corrd[i]<<endl;
	  mom2_file<<isweep<<" "<<mom2<<endl;
	  mag_file<<isweep<<" "<<mag0<<" "<<mag1<<endl;
	  corr_time.stop();
	}
  
      //sweep
      sweep_time.start();
      if(read_pars.use_hmc) hmc_update(isweep<read_pars.nterm);
      else
	{
	  for(int imicro=0;imicro<read_pars.nmicro;imicro++) micro_sweep();
	  overheat_sweep();
	}
      sweep_time.stop();
    }
  
  tot_time.stop();
  
  //write lasted time
  cout<<"Tot time: "<<tot_time<<endl;
  cout<<"Sweep time: "<<sweep_time<<endl;
  cout<<"Geo topo time: "<<geo_topo_time<<endl;
  cout<<"Topo time: "<<topo_time<<endl;
  cout<<"Energy time: "<<energy_time<<endl;
  cout<<"Corr time: "<<corr_time<<endl;
  
  //write the conf
  write_conf("conf",isweep);
  
  print_rand_stat();
  
  //finalize
  close();
  
  return 0;
}
