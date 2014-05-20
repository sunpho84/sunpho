#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <cmath>
#include <iostream>

#include "action.hpp"
#include "data.hpp"
#include "geometry.hpp"
#include "meta.hpp"
#include "random.hpp"

using namespace std;

//initialize the system to hot
void init_system_to_hot()
{for(int s=0;s<V;s++) phi[s]=get_random_z3();}

//initialize to cold
void init_system_to_cold()
{for(int s=0;s<V;s++) phi[s]=0;}

//switch
void init_system_to(int cond)
{
  if(cond==HOT) init_system_to_hot();
  else          init_system_to_cold();
}

//initialize the action contribution for each site
void init_act_contr_tab()
{
  for(int a=0;a<=2*NDIMS;a++)
    for(int b=0;b<=2*NDIMS;b++)
      if(a+b<=6)
	{
	  int icase=a*(2*NDIMS+1)+b;
	  int c=6-a-b;
	  int ntypes[3]={a,b,c};
	  
	  for(int d=0;d<3;d++)
	    {
	      int nequals=ntypes[d];
	      int ndiff=2*NDIMS-nequals;
	      double E=contr_act_equals*nequals+contr_act_diff*ndiff;

	      double M=contr_act_re[d];
	      
	      act_contr_tab[icase][d]=-tau*(E+kappa*M);
	    }
	}
}

//initialize the code
void init(int cond,int seed)
{
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
  for(int mu=0;mu<NDIMS;mu++) V*=L[mu];
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
	  c[mu]=(o+L[mu]-1)%L[mu];
	  neighdw(site,mu)=site_of_coords(c);
	  
	  //forward
	  c[mu]=(o+L[mu]+1)%L[mu];
	  neighup(site,mu)=site_of_coords(c);
	  
	  //restore original
	  c[mu]=o;
	}
    }
  
  //allocate data
  phi=new z3_t[V];
  
  //initialize vector for contr
  int contr_N_tab[3]={0,+1,-1};
  int contr_N0_tab[3]={1,0,0};
#ifdef LUCINI
  double contr_act_re_tab[3]={2,-1,-1};
  contr_act_equals=2;
  contr_act_diff=-1;
#else
  double contr_act_re_tab[3]={1,0,0};
  contr_act_equals=1;
  contr_act_diff=0;
#endif
  
  for(int d=0;d<3;d++)
    {
      contr_N[d]=contr_N_tab[d];
      contr_N0[d]=contr_N0_tab[d];
      contr_act_re[d]=contr_act_re_tab[d];
    }
  
  //set the system to hot state
  init_system_to(cond);
  
  //initialize lookup table
  init_act_contr_tab();
  
  //initialize the action and N
  glb_N=compute_N(phi);
  glb_N0=compute_N0(phi);
  compute_action(phi);
  data_N.resize(ntraj);
}
