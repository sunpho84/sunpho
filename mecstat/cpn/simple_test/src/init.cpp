#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <fstream>
#include <iostream>
#include <omp.h>

#include "data.hpp"
#include "geometry.hpp"
#include "random.hpp"
#include "staples.hpp"
#include "tools.hpp"
#include "topology.hpp"

#define EXTERN_INIT

#include "init.hpp"

using namespace std;

//read the input file
void read_input(read_pars_t &read_pars,const char *path)
{
  //read parameters
  ifstream input(path);
  if(!input.good()) crash("opening input");
  read(N,input,"N");
  read(L,input,"L");
  read(beta,input,"Beta");
  g=1/(N*beta);
  read(read_pars.seed,input,"Seed");
  read(read_pars.nsweep,input,"NSweep");
  string start_cond_str;
  read(start_cond_str,input,"StartCond");
  if(file_exists("conf")) read_pars.start_cond=LOAD;
  else
    if(start_cond_str=="COLD") read_pars.start_cond=COLD;
    else
      if(start_cond_str=="HOT") read_pars.start_cond=HOT;
      else
	if(start_cond_str=="LOAD") read_pars.start_cond=LOAD;
	else crash("Unkwnown start cond %s, use: COLD, HOT, LOAD",start_cond_str.c_str());
  read(read_pars.nterm,input,"NTerm");
  read(compute_corr_each,input,"ComputeCorrEach");
  read(read_pars.use_hmc,input,"UseHMC");
  switch(read_pars.use_hmc)
    {
    case -1:
      break;
    case 0:
       read(read_pars.nmicro,input,"NMicro");
       break;
    case 1:
      read_pars.nmicro=3;
      read(nhmc_steps,input,"NhmcSteps");
      break;
    default:
      crash("unknwon update method %d",read_pars.use_hmc);
    }
  read(use_topo_pot,input,"UseTopoPot");
  switch(use_topo_pot)
    {
    case 0:
      break;
    case 1:
      read(th_top,input,"ThTop");
      break;
    case 2:
      if(!read_pars.use_hmc) crash("must use hmc");
      read(chrono_topo_after,input,"ChronoTopoAfter");
      read(chrono_topo_each,input,"ChronoTopoEach");
      read(chrono_topo_coeff,input,"ChronoTopoCoeff");
      read(chrono_topo_width,input,"ChronoTopoWidth");
      read(chrono_topo_barr,input,"ChronoTopoBarr");
      read(chrono_topo_force_out,input,"ChronoTopoForceOut");
      read(chrono_topo_bend,input,"ChronoTopoBend");
      read(chrono_topo_well_tempering,input,"ChronoTopoWellTempering");
      break;
    }
  read(nstout_lev,input,"NStoutLev");
  read(stout_rho,input,"StoutRho");
}

//initialize the system to hot
void init_system_to_hot()
{
  for(int site=0;site<V;site++)
    {
      //fill the Lambda
      for(int mu=0;mu<NDIMS;mu++)
        set_U1_to_rnd(lambda[site*NDIMS+mu],site);
      
      //fill the Zeta
      set_ON_to_rnd(zeta+site*N,site);
    }
}

//initialize to cold
void init_system_to_cold()
{
#pragma omp parallel for
  for(int site=0;site<V;site++)
    {
      //fill the Lambda
      for(int mu=0;mu<NDIMS;mu++) lambda[site*NDIMS+mu]=1;
      
      //fill the Zeta
      for(int n=0;n<N;n++) zeta[site*N+n]=(n==0);
    }
}

//initialize the code
void init(int &base_isweep,read_pars_t &read_pars)
{
  //take init time
  init_time=time(0);

  //write nthreads
#pragma omp parallel
  {
#pragma omp single
    cout<<omp_get_num_threads()<<" threads"<<endl;
  }
  
  //find parity
  npar=1;
  do npar++;
  while(L%npar!=0);
  cout<<"Parity: "<<npar<<endl;
  
  //geometry
  V=1;
  for(int mu=0;mu<NDIMS;mu++) V*=L;
  cout<<"Volume: "<<V<<endl;
  neigh_data=new int[V*NDIMS*2];
  
  //parity subdivision
  V_per_par=V/npar;
  lx_of_par=new int[V];
  int nof_par[npar];
  for(int ipar=0;ipar<npar;ipar++) nof_par[ipar]=ipar*V_per_par;
  
  //loop over sites
  for(int site=0;site<V;site++)
    {
      //get the original coordinates
      coords c;
      coords_of_site(c,site);
      
      //get parity
      int ipar=0;
      for(int mu=0;mu<NDIMS;mu++) ipar+=c[mu];
      ipar%=npar;
      
      //mark it in the list of parity
      lx_of_par[nof_par[ipar]++]=site;
      
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
  
  //check parity
  for(int ipar=0;ipar<npar;ipar++)
    if(nof_par[ipar]!=(ipar+1)*V_per_par)
      crash("expected %d obtained %d",(ipar+1)*V_per_par,nof_par[ipar]);
  
  //Zeta and Lambda
  zeta=new dcomplex[N*V];
  lambda=new dcomplex[V*NDIMS];

  //Zeta and Lambda for hmc copy
  zeta_old=new dcomplex[N*V];
  lambda_old=new dcomplex[V*NDIMS];
  
  //allocate momenta
  pi=new dcomplex[V*N];
  omega=new double[V*NDIMS];
  
  //allocate force
  fpi=new dcomplex[V*N];
  fomega=new double[V*NDIMS];

  //allocate topo staples
  topo_staples_data=new dcomplex[V*NDIMS];
  topo_staples_supp_data=new dcomplex[V*NDIMS];
  
  //allocate stout lambda
  lambda_stout=new dcomplex*[nstout_lev+1];
  for(int istout_lev=1;istout_lev<=nstout_lev;istout_lev++) lambda_stout[istout_lev]=new dcomplex[V*NDIMS];
  
  //seed the generator
  gen=new mt19937_64[V+1];
  for(int site=0;site<=V;site++)
    gen[site].seed(read_pars.seed+site);
  
  //init according to start condition
  bool conf_loaded=false;
  if(!file_exists("conf"))
    {
      base_isweep=0;
      switch(read_pars.start_cond)
	{
	case HOT:  init_system_to_hot();break;
	case COLD: init_system_to_cold();break;
	case LOAD: crash("conf not exists!");break;
	}
    }
  else
    {
      read_conf(base_isweep,"conf"); //remember rnd gen reinit
      conf_loaded=true;
    }

  if(use_topo_pot==2)
    {
      ngrid=(2*chrono_topo_barr+chrono_topo_width/2)/chrono_topo_width;
      topo_grid.resize(ngrid+1);
      for(auto &grid : topo_grid) grid=0;
      
      //load the potential
      if(conf_loaded) load_chrono_topo_potential();
    }
  
  nacc=0;
}
