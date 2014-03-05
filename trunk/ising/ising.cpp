#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <stdarg.h>
#include <bitset>
#include <omp.h>

using namespace std;

#define ncopies 8
#define nmu 2
#define FW 0
#define BW 1
#define RAN2_NTAB 32

typedef int coords_t[nmu];
typedef int neighs_t[2*nmu];
//typedef bitset<ncopies> spin_t;
typedef bool spin_t[ncopies];

//structure for the random generator
struct rnd_gen_t
{
  int idum;
  int idum2;
  int iv[RAN2_NTAB];
  int iy;
};

int nsweeps,nconfs;                                   //number of sweeps and tot confs to be generated
int L,nsites;                                         //lattice size
spin_t *spins;                                        //configuration
int *next_cluster[ncopies],*curr_cluster[ncopies];    //cluster at the end of the step and at the beginning
neighs_t *neighs;                                     //neighbors
int glb_par_link[ncopies],glb_up_spins[ncopies];      //store the total number of parallel link and up spins
double beta,flip_prob;                                //beta and probablity to flip single spin
rnd_gen_t glb_rnd_gen[ncopies],*loc_rnd_gen;          //random number generators
int total_flipped[ncopies];                           //total number of flipped per copy

/////////////////////////////////////////// utils //////////////////////////////////////

//print spins
void print_spins(int icopy)
{
  cout<<endl;
  for(int x=0;x<L;x++)
    {
      for(int y=0;y<L;y++) cout<<spins[x+L*y][icopy]<<" ";
      cout<<endl;
    }
  cout<<endl;
}

//crash reporting the expanded error message
void crash(const char *templ,...)
{
  //expand error message
  char mess[1024];
  va_list ap;
  va_start(ap,templ);
  vsprintf(mess,templ,ap);
  va_end(ap);
  
  fprintf(stderr,"ERROR: \"%s\".\n",mess);
  exit(1);
}

//open a file
FILE *open_file(const char* path,const char* mod)
{
  FILE *out=fopen(path,mod);
  if(out==NULL) crash("opening file '%s' in mode '%s'",path,mod);
  
  return out;
}

///////////////////////////////////////// randoms /////////////////////////////////////

double rnd_get_unif(rnd_gen_t *gen,double min,double max);
  
//initialize a random number generator
void start_rnd_gen(rnd_gen_t *out,int seed)
{
  const int im1=2147483563,ia1=40014;
  const int iq1=53668,ir1=12211;
  int j,k;
    
  //initialization
  out->idum=seed;
  out->idum=std::max(out->idum+1,1);
  out->idum2=out->idum;
  for(j=RAN2_NTAB+7;j>=0;j--)
    {
      k=out->idum/iq1;
      out->idum=ia1*(out->idum-k*iq1)-k*ir1;
      if(out->idum<0) out->idum+=im1;
      if(j<RAN2_NTAB) out->iv[j]=out->idum;
    }
  out->iy=out->iv[0];
}

//initialize the global random generator
void start_glb_rnd_gen(int seed)
{
  for(int icopy=0;icopy<ncopies;icopy++) start_rnd_gen(glb_rnd_gen+icopy,seed+icopy);
}

//standard ran2 from numerical recipes
double rnd_get_unif(rnd_gen_t *gen,double min,double max)
{
  const int im1=2147483563,im2=2147483399,imm1=im1-1,ia1=40014,ia2=40692;
  const int iq1=53668,iq2=52774,ir1=12211,ir2=3791,ndiv=1+imm1/RAN2_NTAB;
  const double am=1.0/im1,eps=1.2e-7,rnmx=1-eps;
  int j,k;
  double out;
    
  k=gen->idum/iq1;
  gen->idum=ia1*(gen->idum-k*iq1)-k*ir1;
  if(gen->idum<0) gen->idum+=im1;
    
  k=gen->idum2/iq2;
  gen->idum2=ia2*(gen->idum2-k*iq2)-k*ir2;
  if(gen->idum2<0) gen->idum2+=im2;
    
  j=gen->iy/ndiv;
  gen->iy=gen->iv[j]-gen->idum2;
  gen->iv[j]=gen->idum;
  if(gen->iy<0) gen->iy+=imm1;
  
  out=std::min(am*gen->iy,rnmx);
    
  return out*(max-min)+min;
}
//return a numer between 0 and 1
int rnd_get_pm_one(rnd_gen_t *gen)
{
  double r=rnd_get_unif(gen,0,1);
  if(r>0.5) return 1;
  else return -1;
}

//initialize the grid of local random number generator
void start_loc_rnd_gen(int seed)
{
  //starting global random generator
  start_glb_rnd_gen(seed);
    
  //allocate the grid of random generator, one for site and copy
  loc_rnd_gen=new rnd_gen_t[nsites*ncopies];
  for(int icopy=0;icopy<ncopies;icopy++)
    {
      int internal_seed=(int)rnd_get_unif(glb_rnd_gen+icopy,0,RAND_MAX);
      for(int site=0;site<nsites;site++)
	start_rnd_gen(&(loc_rnd_gen[nsites*icopy+site]),internal_seed+nsites*icopy+site);
    }
  //printf("Grid of local random generators initialized with internal seed: %d\n",internal_seed);
}

//////////////////////////////////////// geometry /////////////////////////////////////

//get coords of site
void coords_of_site(coords_t coords,int site)
{
  for(int mu=0;mu<nmu;mu++)
    {
      coords[mu]=site%L;
      site/=L;
    }
}

//return site of coords
int site_of_coords(coords_t x)
{
  int site=0;
  for(int mu=0;mu<nmu;mu++) site=site*L+x[mu];
  return site;
}

//////////////////////////////////// main program ////////////////////////////

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
  
  //start the local random generator
  start_loc_rnd_gen(seed);
  
  //allocate array
  spins=new spin_t[nsites];
  for(int icopy=0;icopy<ncopies;icopy++)
    {
      next_cluster[icopy]=new int[nsites];
      curr_cluster[icopy]=new int[nsites];
    }
  neighs=new neighs_t[nsites];
  
  //set the geometry
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
  
  //reset the number of parallel and up spins
  for(int icopy=0;icopy<ncopies;icopy++)
    {
      glb_par_link[icopy]=0;
      glb_up_spins[icopy]=0;
    }
  
  //generate the configuration
  for(int site=0;site<nsites;site++)
    for(int icopy=0;icopy<ncopies;icopy++)
      spins[site][icopy]=(icopy%2)?0:(bool)(rnd_get_unif(loc_rnd_gen+site*ncopies+icopy,0,1)>=0.5);
  
  //count the "up" sites and number of parallel sites
  for(int site=0;site<nsites;site++)
    for(int icopy=0;icopy<ncopies;icopy++)
      {
	glb_up_spins[icopy]+=spins[site][icopy];
	for(int mu=0;mu<nmu;mu++) if(spins[neighs[site][mu]][icopy]==spins[site][icopy]) glb_par_link[icopy]++;
      }
  
  //define the probability
  flip_prob=1-exp(-2*beta);
}

//change a single cluster using Wolf algorithm
int change_single_cluster(int icopy)
{
  //cluster size, for future reference
  int cluster_size=1;
  
  //take the site and flip it
  int site=(int)rnd_get_unif(glb_rnd_gen+icopy,0,nsites);
  bool nse=!spins[site][icopy];
  spins[site][icopy]=nse;
  
  //update parallel sites surrounding
  for(int mu=0;mu<2*nmu;mu++)
    if(spins[neighs[site][mu]][icopy]==nse) glb_par_link[icopy]++;
    else glb_par_link[icopy]--;

  //set in the cluster
  next_cluster[icopy][0]=site;
  int pc=1;

  do
    {
      //swap current cluster with next
      swap(next_cluster[icopy],curr_cluster[icopy]);
      int gc=pc;

      pc=0;

      for(int g=0;g<gc;g++)
	{
	  site=curr_cluster[icopy][g];

	  //add neighbors
	  for(int mu=0;mu<2*nmu;mu++)
	    {
	      int p=neighs[site][mu];

	      //if p is antiparallel, consider it
	      if(spins[p][icopy]!=nse)
		{
		  //try to add it
		  if(rnd_get_unif(loc_rnd_gen+site*ncopies+icopy,0,1)<flip_prob)
		    {
		      //flip it
		      spins[p][icopy]=nse;
		      for(int nu=0;nu<2*nmu;nu++)
			if(spins[neighs[p][nu]][icopy]==nse) glb_par_link[icopy]++;
			else glb_par_link[icopy]--;

		      //increase cluster size
		      cluster_size++;
		      
		      //put it in the list of "to be seen"
		      next_cluster[icopy][pc]=p;
		      pc++;
		    }
		}

	    }
	}
    }
  while(pc>0);
  
  //adjust up count
  if(nse==1) glb_up_spins[icopy]+=cluster_size;
  else glb_up_spins[icopy]-=cluster_size;
  
  return cluster_size;
}

//run the simulation
void run()
{
  //reset the number of totally flipped
  for(int icopy=0;icopy<ncopies;icopy++) total_flipped[icopy]=0;

  for(int iconf=0;iconf<nconfs;iconf++)
    {
      double magnetizz[ncopies];
      //sweep
#pragma omp parallel for private(nsweeps)
      for(int icopy=0;icopy<ncopies;icopy++)
	{
	  int nsweeps=0;
	  do
	    {
	      total_flipped[icopy]+=change_single_cluster(icopy);
	      nsweeps++;
	    }
	  while(nsweeps<nsites && total_flipped[icopy]<iconf*nsites);
	  
	  //check that we did a number of sweeps smaller than the total number of sites
	  if(nsweeps>nsites) crash("something went wront, nsweeps: %d, nsites: %d",nsweeps,nsites);
	  
	  //compute energy and magnetization
	  //double energy=-glb_par_link*2;
	  magnetizz[icopy]=((glb_up_spins[icopy]*2)-nsites)/(double)(nsites);
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
  delete [] spins;
  for(int icopy=0;icopy<ncopies;icopy++)
    {
      delete [] next_cluster[icopy];
      delete [] curr_cluster[icopy];
    }
  delete [] neighs;
  delete [] loc_rnd_gen;
}

int main(int narg,char**arg)
{
  if(narg<2) crash("use %s input",arg[0]);
  
  init(arg[1]);
  run();
  close();
  
  return 0;
}
