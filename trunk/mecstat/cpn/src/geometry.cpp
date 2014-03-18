#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "debug.hpp"
#include "geometry.hpp"
#include "simul.hpp"
#include "utils.hpp"

//compute internal volume
int bulk_volume(coords_t L)
{
  int intvol=1,mu=0;
    do
      {
        if(L[mu]>2) intvol*=L[mu]-2;
        else intvol=0;
        
        mu++;
      }
    while(intvol!=0 && mu<NMU);
    
    return intvol;
}

//compute the bulk volume of the local lattice, given by L/R
int bulk_recip_lat_volume(coords_t R,coords_t L)
{
  coords_t X;
  for(int mu=0;mu<NMU;mu++) X[mu]=L[mu]/R[mu];
  return bulk_volume(X);
}

//compute the variance of the border
int compute_border_variance(coords_t L,coords_t P,int factorize_processor)
{
  int S2B=0,SB=0;
  for(int ib=0;ib<NMU;ib++)
    {
      int B=1;
      for(int mu=0;mu<NMU;mu++) if(mu!=ib) B*=(factorize_processor) ? L[mu]/P[mu] : P[mu];
      SB+=B;
      S2B+=B*B;
    }
  SB/=NMU;
  S2B/=NMU;
  S2B-=SB*SB;
    
  return S2B;
}

//get coords of site
void geometry_t::coords_of_site(coords_t coords,int site,coords_t sizes)
{
  for(int mu=0;mu<NMU;mu++)
    {
      coords[mu]=site%sizes[mu];
      site/=sizes[mu];
    }
}

//return site of coords
int geometry_t::site_of_coords(coords_t coords,coords_t sizes)
{
  int site=0;
  for(int mu=0;mu<NMU;mu++) site=site*sizes[mu]+coords[mu];
  return site;
}

//divide the global lattice into ranks in homogeneous way
void geometry_t::partition_ranks_homogeneously()
{
  bool something_found=false;
  
  int NR=simul->nranks;
  int LV=glb_size/NR;
  
  //minimal variance border
  int mBV=-1;
  
  //factorize the nodes
  int list_fact_LV[log2N(LV)];
  int nfact_LV=factorize(list_fact_LV,LV);
        
  //factorize the number of rank
  int list_fact_NR[log2N(NR)];
  int nfact_NR=factorize(list_fact_NR,NR);
        
  //if nfact_LV>=nfact_NR factorize the number of rank, otherwise the local volume
  //in the first case we find the best way to assign the ranks to different directions
  //in the second case we find how many sites per direction to assign to each rank
  int factorize_rank=(nfact_LV>=nfact_NR);
  int nfact=factorize_rank ? nfact_NR : nfact_LV;
  int *list_fact=factorize_rank ? list_fact_NR : list_fact_LV;
  
  //compute the number of combinations: this is given by NMU^nfact
  int ncombo=1;
  for(int ifact=0;ifact<nfact;ifact++) ncombo*=NMU;
  
  //find the partition which minimize the surface and the surface variance
  int min_surf_LV=-1;
  int icombo=0;
  for(int mu=0;mu<NMU;mu++) nranks_per_dir[mu]=-1;
        
  do
    {
      //number of ranks in each direction for current partitioning
      coords_t R;
      for(int mu=0;mu<NMU;mu++) R[mu]=1;
      
      //find the partioning corresponding to icombo
      int ifact=nfact-1;
      int valid_partitioning=1;
      int jcombo=icombo;
      do
	{
	  //find the direction: this is given by the ifact digit of icombo wrote in base NMU
	  int mu=(jcombo&(NMU-1));
	  jcombo/=NMU;
          
	  //if we are factorizing local lattice, rank factor is given by list_fact, otherwise L/list_fact
	  R[mu]*=list_fact[ifact];
          
	  //check that the total volume glb_sizes[mu] is a multiple and it is larger than the number of proc
	  valid_partitioning=(glb_sizes[mu]%R[mu]==0 && glb_sizes[mu]>=R[mu]);
	  if(valid_partitioning) ifact--;
	}
      while(valid_partitioning && ifact>=0);
            
      //validity could have changed
      if(valid_partitioning)
	{
	  //if we are factorizing reciprocal lattice, convert back to rank grid
	  for(int mu=0;mu<NMU;mu++)
	    if(!factorize_rank)
	      R[mu]=glb_sizes[mu]/R[mu];
	  
	  //compute the surface=loc_vol-bulk_volume
	  int BV=bulk_recip_lat_volume(R,glb_sizes);
	  int surf_LV=LV-BV;
          
	  //look if this is the new minimal surface
	  int new_minimal=0;
	  //if it is the minimal surface (or first valid combo) copy it and compute the border size
	  if(surf_LV<min_surf_LV||min_surf_LV==-1)
	    {
	      new_minimal=1;
	      mBV=compute_border_variance(glb_sizes,R,factorize_rank);
	    }
	  //if it is equal to previous found surface, consider borders variance
	  if(surf_LV==min_surf_LV)
	    {
	      int BV=compute_border_variance(glb_sizes,R,factorize_rank);
	      //if borders are more homogeneus consider this grid
	      if(BV<mBV)
		{
		  mBV=BV;
		  new_minimal=1;
		}
	    }
	  
	  //save it as new minimal
	  if(new_minimal)
	    {
	      min_surf_LV=surf_LV;
	      for(int mu=0;mu<NMU;mu++) nranks_per_dir[mu]=R[mu];
	      something_found=1;
	    }
	  
	  icombo++;
	}
      //skip all remaining factorization using the same structure
      else
	{
	  int skip=1;
	  for(int jfact=0;jfact<ifact-1;jfact++) skip*=NMU;
	  icombo+=skip;
	}
    }
  while(icombo<ncombo);
  
  if(!something_found) CRASH("no valid partitioning found");
}

//initialize the geometry
void geometry_t::start(coords_t ext_glb_sizes)
{
  //copy glb sizes and compute total volume
  glb_size=1;
  for(int mu=0;mu<NMU;mu++)
    {
      glb_sizes[mu]=ext_glb_sizes[mu];
      glb_size*=glb_sizes[mu];
    }
  
  //output global size and number of ranks
  MASTER_PRINTF("Global lattice grid:\t%d",glb_sizes[0]);
  for(int mu=1;mu<NMU;mu++) MASTER_PRINTF("x%d",glb_sizes[mu]);
  MASTER_PRINTF(" = %d\n",glb_size);
  MASTER_PRINTF("Number of running ranks: %d\n",simul->nranks);
    
  //check wheter each ranks will have the same number of sites
  homogeneous_partitioning=(glb_size%simul->nranks==0);
  
  //check basical things, such as that we are not using more ranks than sites
  if(simul->nranks>glb_size) CRASH("too many ranks!");
  
  //on the base of the results, partition
  if(homogeneous_partitioning) partition_ranks_homogeneously();
  else CRASH("not implemented yet");
  MASTER_PRINTF("Creating ranks grid:\t%d",nranks_per_dir[0]);
  for(int mu=1;mu<NMU;mu++) MASTER_PRINTF("x%d",nranks_per_dir[mu]);
  MASTER_PRINTF("\n");
}
