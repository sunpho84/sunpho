#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <fstream>
#include <iostream>
#include <omp.h>

#include "geometry.hpp"
#include "random.hpp"
#include "tools.hpp"

#define EXTERN_INIT

#include "init.hpp"

using namespace std;

//read the input file
void read_input(const char *path)
{
  //read parameters
  ifstream input(path);
  if(!input.good()) crash("opening input");
  read(N,input,"N");
  read(L,input,"L");
  read(seed,input,"Seed");
  read(nsweep,input,"NSweep");
  string start_cond_str;
  read(start_cond_str,input,"StartCond");
  if(file_exists("conf")) start_cond=LOAD;
  else
    if(start_cond_str=="COLD") start_cond=COLD;
    else
      if(start_cond_str=="HOT") start_cond=HOT;
      else
	if(start_cond_str=="LOAD") start_cond=LOAD;
	else crash("Unkwnown start cond %s, use: COLD, HOT, LOAD",start_cond_str.c_str());
  read(nterm,input,"NTerm");
  read(nhmc_steps,input,"NhmcSteps");
}

//initialize the code
void init()
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
  
  //seed the generator
  gen=new mt19937_64[V+1];
  for(int site=0;site<=V;site++)
    gen[site].seed(seed+site);
}
