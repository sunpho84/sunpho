#include "global_variables.hpp"
#include "system.hpp"

#include <algorithm>
#include <iostream>
#include <stdlib.h>

using namespace std;

//initialization
system_t::system_t(int seed)
{
  //start the global random generator
  glb_rnd_gen.init(seed);
  
  //allocate the grid of random generator, one for site, and start them
  loc_rnd_gen=new rnd_gen_t[nsites];
  int internal_seed=(int)glb_rnd_gen.rnd_get_unif(0,RAND_MAX);
  for(int site=0;site<nsites;site++) loc_rnd_gen[site].init(internal_seed+site);
  
  //allocate configuration and cluster
  spins=new spin_t[nsites];
  next_cluster=new int[nsites];
  curr_cluster=new int[nsites];
  
  //generate the configuration
  for(int site=0;site<nsites;site++) spins[site]=(bool)(loc_rnd_gen[site].rnd_get_unif(0,1)>=0.5);
  
  //count the "up" sites and number of parallel sites
  glb_par_link=glb_up_spins=0;
  for(int site=0;site<nsites;site++)
    {
      glb_up_spins+=spins[site];
      for(int mu=0;mu<nmu;mu++) if(spins[neighs[site][mu]]==spins[site]) glb_par_link++;
    }
}

//change a single cluster using Wolf algorithm
int system_t::change_single_cluster()
{
  //cluster size, for future reference
  int cluster_size=1;
  
  //take the site and flip it
  int site=(int)glb_rnd_gen.rnd_get_unif(0,nsites);
  bool nse=!spins[site];
  spins[site]=nse;
  
  //update parallel sites surrounding
  for(int mu=0;mu<2*nmu;mu++)
    if(spins[neighs[site][mu]]==nse) glb_par_link++;
    else glb_par_link--;
  
  //set in the cluster
  next_cluster[0]=site;
  int pc=1;
  
  do
    {
      //swap current cluster with next
      swap(next_cluster,curr_cluster);
      int gc=pc;
      
      pc=0;
      for(int g=0;g<gc;g++)
	{
	  site=curr_cluster[g];
	  
	  //add neighbors
	  for(int mu=0;mu<2*nmu;mu++)
	    {
	      int p=neighs[site][mu];
	      
	      //if p is antiparallel, consider it
	      if(spins[p]!=nse)
		{
		  //try to add it
		  if(loc_rnd_gen[site].rnd_get_unif(0,1)<flip_prob)
		    {
		      //flip it
		      spins[p]=nse;
		      for(int nu=0;nu<2*nmu;nu++)
			if(spins[neighs[p][nu]]==nse) glb_par_link++;
			else glb_par_link--;
		      
		      //increase cluster size
		      cluster_size++;
		      
		      //put it in the list of "to be seen"
		      next_cluster[pc]=p;
		      pc++;
		    }
		}
	      
	    }
	}
    }
  while(pc>0);
  
  //adjust up count
  if(nse==1) glb_up_spins+=cluster_size;
  else glb_up_spins-=cluster_size;
  
  return cluster_size;
}

//print spins
void system_t::print_spins()
{
  cout<<endl;
  for(int x=0;x<L;x++)
    {
      for(int y=0;y<L;y++) cout<<spins[x+L*y]<<" ";
      cout<<endl;
    }
  cout<<endl;
}

//destructor
system_t::~system_t()
{
  delete [] spins;
  delete [] next_cluster;
  delete [] curr_cluster;
  delete [] loc_rnd_gen;
}
