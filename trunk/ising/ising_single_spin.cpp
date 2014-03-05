#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

using namespace std;

#define nmu 2
#define FW 0
#define BW 1
#define RAN2_NTAB 32

typedef int coords_t[nmu];
typedef int neighs_t[2*nmu];

//The structure for the random generator
struct rnd_gen_t
{
  int idum;
  int idum2;
  int iv[RAN2_NTAB];
  int iy;
};

int nsweeps,nconfs;                 //number of sweeps and tot confs to be generated
int L,nsites;                       //lattice size
bool *spins;                        //configuration
int *next_cluster,*curr_cluster;    //cluster at the end of the step and at the beginning
neighs_t *neighs;                   //neighbors
int glb_par_link,glb_up_spins;      //store the total number of parallel link and up spins
double beta,flip_prob;              //beta and probablity to flip single spin
rnd_gen_t glb_rnd_gen,*loc_rnd_gen; //random number generators

/////////////////////////////////////////// utils //////////////////////////////////////

//print spins
void print_spins()
{
  cout<<endl;
  for(int x=0;x<L;x++)
    {
      for(int y=0;y<L;y++) cout<<spins[x+L*y]<<" ";
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
  start_rnd_gen(&(glb_rnd_gen),seed);
  //printf("Global random generators initialized with seed: %d\n",seed);
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
  //Generate the true seed
  start_glb_rnd_gen(seed);
  int internal_seed=(int)rnd_get_unif(&glb_rnd_gen,0,RAND_MAX);
    
  //allocate the grid of random generator, one for site
  loc_rnd_gen=new rnd_gen_t[nsites];
  for(int site=0;site<nsites;site++) start_rnd_gen(&(loc_rnd_gen[site]),internal_seed+site);
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
  spins=new bool[nsites];
  next_cluster=new int[nsites];
  curr_cluster=new int[nsites];
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
  glb_par_link=0;
  glb_up_spins=0;
  
  //generate the configuration
  for(int site=0;site<nsites;site++) spins[site]=(rnd_get_unif(loc_rnd_gen+site,0,1)>=0.5);
  
  //count the "up" sites and number of parallel sites
  for(int site=0;site<nsites;site++)
    {
      glb_up_spins+=spins[site];
      for(int mu=0;mu<nmu;mu++) if(spins[neighs[site][mu]]==spins[site]) glb_par_link++;
    }
  
  //define the probability
  flip_prob=1-exp(-2*beta);
}

//change a single cluster using Wolf algorithm
void change_single_cluster()
{
  //cluster size, for future reference
  int cluster_size=1;
  
  //take the site and flip it
  int site=(int)rnd_get_unif(&glb_rnd_gen,0,nsites);
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
		  if(rnd_get_unif(loc_rnd_gen+site,0,1)<flip_prob)
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
}

//run the simulation
void run()
{
  for(int iconf=0;iconf<nconfs;iconf++)
    {
      //sweep
      for(int isweep=0;isweep<nsweeps;isweep++) change_single_cluster();
      
      //compute energy and magnetization
      //double energy=-glb_par_link*2;
      double magnetizz=((glb_up_spins*2)-nsites)/(double)(nsites);
      cout<<magnetizz<<endl;
    }
}

//free used arrays
void close()
{
  delete [] spins;
  delete [] next_cluster;
  delete [] curr_cluster;
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
