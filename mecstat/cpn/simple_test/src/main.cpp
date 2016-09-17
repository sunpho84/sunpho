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
  ofstream charge_file("charge",mode);
  ofstream energy_file("energy",mode);
  ofstream polyakov_file("polyakov",mode);
  ofstream topology_file("topology",mode);
  ofstream corr_file("corr",mode);
  ofstream corrd_file("corrd",mode);
  ofstream mom2_file("mom2",mode);
  ofstream mag_file("mag",mode);
  charge_file.precision(PREC);
  energy_file.precision(PREC);
  polyakov_file.precision(PREC);
  topology_file.precision(PREC);
  corr_file.precision(PREC);
  corrd_file.precision(PREC);
  mag_file.precision(PREC);
  
  //sweep with overheat-micro
  timing_t tot_time,sweep_time,charge_time,energy_time,geo_topo_time,topo_time,corr_time;
  
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
      
      //compute charge
      charge_time.start();
      dcomplex C=charge(zeta,lambda);
      if(use_charge_pot==2)
	{
	  chrono_charge.update(isweep,C.imag());
	  if(isweep%DRAW_EACH==0)
	    {
	      chrono_charge.save();
	      chrono_charge.draw_force("charge_force");
	    }
	}
      charge_file<<isweep<<" "<<C.real()<<" "<<C.imag()<<endl;
      charge_time.stop();
      
      //compute polyakov loop
      dcomplex poly=polyakov(lambda);
      polyakov_file<<isweep<<" "<<poly.real()<<" "<<poly.imag()<<endl;
      
      //compute topological charge
      topo_time.start();
      stout_lambda_whole_stack(lambda_stout,stout_rho,nstout_lev,lambda);
      for(int ilev=0;ilev<=nstout_lev;ilev++)
	{
	  double topo_num=topology(lambda_stout[ilev]);
	  if(use_topo_pot==2 && ilev==nstout_lev)
	    {
	      chrono_topo.update(isweep,+topo_num);
	      if(isweep%DRAW_EACH==0)
		{
		  chrono_topo.save();
		  chrono_topo.draw_force("topo_force");
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
	  //compute_corr_alt(corr,zeta);
	  compute_corr(mag0,mag1,mom2,corr,corrd,zeta);
	  for(int i=0;i<=L/2;i++) corr_file<<isweep<<" "<<i<<" "<<corr[i]<<endl;
	  for(int i=0;i<=L/2;i++) corrd_file<<isweep<<" "<<i/sqrt(2)<<" "<<corrd[i]<<endl;
	  mom2_file<<isweep<<" "<<mom2<<endl;
	  mag_file<<isweep<<" "<<mag0<<" "<<mag1<<endl;
	  corr_time.stop();
	}
      
      //sweep
      sweep_time.start();
      switch(read_pars.use_hmc)
	{
	case 1:
	   hmc_update(isweep<read_pars.nterm);
	   break;
	case 0:
	  for(int imicro=0;imicro<read_pars.nmicro;imicro++) micro_sweep();
	  overheat_sweep();
	  break;
	default:
	  crash("unkwnown update %d",read_pars.use_hmc);
	}
      sweep_time.stop();
    }
  
  tot_time.stop();
  
  //write lasted time
  cout<<"Acc: "<<nacc/(double)read_pars.nsweep<<endl;
  cout<<"Tot time: "<<tot_time<<endl;
  cout<<"Sweep time: "<<sweep_time<<endl;
  cout<<"Geo topo time: "<<geo_topo_time<<endl;
  cout<<"Topo time: "<<topo_time<<endl;
  cout<<"Energy time: "<<energy_time<<endl;
  cout<<"charge time: "<<charge_time<<endl;
  cout<<"Corr time: "<<corr_time<<endl;
  
  //write the conf
  write_conf("conf",isweep);
  if(use_topo_pot==2) chrono_topo.save();
  if(use_charge_pot==2) chrono_charge.save();
  
  print_rand_stat();
  
  //finalize
  close();
  
  return 0;
}
