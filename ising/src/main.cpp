#include <iostream>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include "global_variables.hpp"
#include "system.hpp"
#include "utils.hpp"

using namespace std;

//initialize the program
void init(const char *path)
{
  //scan the file
  FILE *fin=open_file(path,"r");
  if(fscanf(fin,"L %d\n",&L)!=1) crash("reading L");
  int seed;
  if(fscanf(fin,"seed %d\n",&seed)!=1) crash("reading seed");
  if(fscanf(fin,"beta %lg\n",&beta)!=1) crash("reading beta");
  if(fscanf(fin,"nconfs %d\n",&nconfs)!=1) crash("reading nconfs");
  if(fscanf(fin,"nsweeps %d\n",&nsweeps)!=1) crash("reading nsweeps");
  fclose(fin);
    
  //compute nsites
  nsites=L;
  for(int mu=1;mu<nmu;mu++) nsites*=L;
  
  //set the geometry
  neighs=new neighs_t[nsites];
#pragma omp parallel for
  for(int site=0;site<nsites;site++)
    {
      coords_t coords;
      coords_of_site(coords,site);
      
      //set neighbors
      for(int mu=0;mu<nmu;mu++)
        {
          int cmu=coords[mu];

          //forward neighbors
          coords[mu]=(cmu+1)%L;
          neighs[site][2*mu+0]=site_of_coords(coords);
          coords[mu]=(cmu+L-1)%L;
          neighs[site][2*mu+1]=site_of_coords(coords);
          coords[mu]=cmu;
        }
    }
  
  //generate the systems and initialize them
  for(int icopy=0;icopy<ncopies;icopy++) systems[icopy]=new system_t(seed+icopy);
  
  //define the probability
  flip_prob=1-exp(-2*beta);
}

//run the simulation
void run()
{
  //reset the number of totally flipped
  int total_flipped[ncopies];
  for(int icopy=0;icopy<ncopies;icopy++) total_flipped[icopy]=0;
  
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      double magnetizz[ncopies];
      //sweep
#pragma omp parallel for
      for(int icopy=0;icopy<ncopies;icopy++)
	{
	  int nsweeps=0;
	  do
	    {
	      total_flipped[icopy]+=systems[icopy]->change_single_cluster();
	      nsweeps++;
	    }
	  while(nsweeps<nsites && total_flipped[icopy]<iconf*nsites);
	  
	  //check that we did a number of sweeps smaller than the total number of sites
	  if(nsweeps>nsites) crash("something went wront, nsweeps: %d, nsites: %d",nsweeps,nsites);
	  
	  //compute energy and magnetization
	  //double energy=-glb_par_link*2;
	  magnetizz[icopy]=(systems[icopy]->glb_up_spins*2-nsites)/(double)(nsites);
	  //cout<<nsweeps/(double)nsites<<" ";
	}
      
      //print magnetizzation
      for(int icopy=0;icopy<ncopies;icopy++) cout<<magnetizz[icopy]<<" ";
      cout<<endl;
    }
}

//free used arrays
void close()
{
  delete [] neighs;
  for(int icopy=0;icopy<ncopies;icopy++) delete systems[icopy];
}

int main(int narg,char**arg)
{
  if(narg<2) crash("use %s input",arg[0]);
  
  init(arg[1]);
  run();
  close();
  
  return 0;
}
