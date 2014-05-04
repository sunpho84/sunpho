#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <iostream>

#include "data.hpp"
#include "geometry.hpp"
#include "random.hpp"

using namespace std;

int init_time;

//initialize the code
void init(int cond,int seed)
{
  init_time=time(0);
  
#ifdef GOOD_GENERATOR
  //init the random generators
  rd=new random_device();
  gen=new mt19937_64((*rd)());
  gen->seed(seed);
  dis=new uniform_real_distribution<double>;
#else
  const int im1=2147483563,ia1=40014;
  const int iq1=53668,ir1=12211;
  int j,k;
  
  //initialization
  gen.idum=seed;
  gen.idum=std::max(gen.idum+1,1);
  gen.idum2=gen.idum;
  for(j=RAN2_NTAB+7;j>=0;j--)
    {
      k=gen.idum/iq1;
      gen.idum=ia1*(gen.idum-k*iq1)-k*ir1;
      if(gen.idum<0) gen.idum+=im1;
      if(j<RAN2_NTAB) gen.iv[j]=gen.idum;
    }
  gen.iy=gen.iv[0];
#endif
    
  //geometry
  V=1;
  for(int mu=0;mu<NDIMS;mu++) V*=L;
  cout<<"Volume: "<<V<<endl;
  neigh_data=new int[V*NDIMS*2];
  
  //loop over sites
  for(int site=0;site<V;site++)
    {
      //get the original coordinates
      coords c;
      coords_of_site(c,site);
      
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
  
  //set the system to hot state
  init_system_to(cond);
}
