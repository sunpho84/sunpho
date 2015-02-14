#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/new_types_definitions.hpp"
#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

namespace bissa
{  
  double rnd_get_unif(rnd_gen *gen,double min,double max);
  
  //initialize a random number generator
  void start_rnd_gen(rnd_gen *out,int seed)
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
  
  //print all the entries of the random generator into a string
  void convert_rnd_gen_to_text(char *text,rnd_gen *gen)
  {
    int *ptr=(int*)gen;
    
    //init output text
    sprintf(text,"%d",ptr[0]);
    
    //append all the elements
    for(int i=1;i<RAN2_NTAB+3;i++)
      {
	char temp[20];
	sprintf(temp," %d",ptr[i]);
	strcat(text,temp);
      }
  }
  
  //read all the entries of the random generator from a string
  void convert_text_to_rnd_gen(rnd_gen *gen,char *text)
  {
    int *ptr=(int*)gen;
    
    for(int i=0;i<RAN2_NTAB+3;i++)
      {
	char temp[20];
	if(sscanf(text,"%s",temp)!=1) crash("while reading element %d from %s",i,text);
	text+=1+strlen(temp);
	if(sscanf(temp,"%d",ptr+i)!=1) crash("while converting to int %s",temp);
      }
  }
  
  //initialize the global random generator
  void start_glb_rnd_gen(int seed)
  {
    if(glb_rnd_gen_inited==1) crash("global random generator already initialized");
    start_rnd_gen(&(glb_rnd_gen),seed);
    
    glb_rnd_gen_inited=1;
    master_printf("Global random generators initialized with seed: %d\n",seed);
  }
  
  //init from text
  void start_glb_rnd_gen(char *text)
  {
    if(glb_rnd_gen_inited==1) crash("global random generator already initialized");
    convert_text_to_rnd_gen(&(glb_rnd_gen),text);
    
    glb_rnd_gen_inited=1;
    master_printf("Global random generators initialized from text\n");
  }
  
  //initialize the grid of local random number generator
  void start_loc_rnd_gen(int seed)
  {
    if(loc_rnd_gen_inited==1) crash("local random generator already initialized!");
    
    //check the grid to be initiaized
    if(loc_vol==0) crash("grid not initalized!");
    
    //Generate the true seed
    if(glb_rnd_gen_inited==0) start_glb_rnd_gen(seed);
    int internal_seed=(int)rnd_get_unif(&glb_rnd_gen,0,RAND_MAX);
    
    //allocate the grid of random generator, one for site
    loc_rnd_gen=bissa_malloc("Loc_rnd_gen",loc_vol,rnd_gen);
    for(int ivol=0;ivol<loc_vol;ivol++) start_rnd_gen(&(loc_rnd_gen[ivol]),internal_seed+glblx_of_loclx[ivol]);
    
    loc_rnd_gen_inited=1;
    master_printf("Grid of local random generators initialized with internal seed: %d\n",internal_seed);
  }
  
  //init from text
  void start_loc_rnd_gen(char *text)
  {
    start_glb_rnd_gen(text);
    start_loc_rnd_gen(0); //glb grid already started so does not matter
  }
  
  //stop grid of local random generators
  void stop_loc_rnd_gen()
  {
    if(loc_rnd_gen_inited==0) crash("local random generator not initialized");
    master_printf("Stopping local random generators\n");
    
    bissa_free(loc_rnd_gen);
    loc_rnd_gen_inited=0;
  }
  
  //standard ran2 from numerical recipes
  double rnd_get_unif(rnd_gen *gen,double min,double max)
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
  int rnd_get_pm_one(rnd_gen *gen)
  {
    double r=rnd_get_unif(gen,0,1);
    if(r>0.5) return 1;
    else return -1;
  }
  
  //fill a grid of vectors with numbers between 0 and 1
  THREADABLE_FUNCTION_4ARG(rnd_fill_unif_loc_vector, double*,v, int,dps, double,min, double,max)
  {
    GET_THREAD_ID();
    BISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int i=0;i<dps;i++)
	v[ivol*dps+i]=rnd_get_unif(&(loc_rnd_gen[ivol]),min,max);
    
    set_borders_invalid(v);
  }
  THREADABLE_FUNCTION_END
  
  //return a grid of +-x numbers
  THREADABLE_FUNCTION_2ARG(rnd_fill_pm_one_loc_vector, double*,v, int,nps)
  {
    GET_THREAD_ID();
    BISSA_PARALLEL_LOOP(ivol,0,loc_vol)
      for(int i=0;i<nps;i++)
	v[ivol*nps+i]=rnd_get_pm_one(&(loc_rnd_gen[ivol]));
    
    set_borders_invalid(v);
  }
  THREADABLE_FUNCTION_END
}
