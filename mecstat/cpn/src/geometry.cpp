#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <math.h>

#include "debug.hpp"
#include "geometry.hpp"
#include "simul.hpp"
#include "utils.hpp"

//return the size of a particular rank
coords_t geometry_t::loc_sizes_of_rank_coords(coords_t in_coords)
{
  coords_t finishing(ndims);
  coords_t starting=glb_coords_of_origin_of_rank_coords(rank_coords);
  for(size_t dim=0;dim<ndims;dim++)
    {
      finishing[dim]=std::min(loc_comm_sizes[dim]*(rank_coords[dim]+1),glb_sizes[dim]);
      if(rank_coords[dim]+1==nranks_per_dir[dim]) finishing[dim]=glb_sizes[dim];
    }
  
  return finishing-starting;
}

//return the rank hosting site of passed cords, and local index
void geometry_t::rank_and_loc_site_of_glb_coords(int &dest_rank,int &dest_site,coords_t dest_glb_coords)
{
  //find rank coords
  coords_t dest_rank_coords=dest_glb_coords/loc_comm_sizes;
  for(size_t dim=0;dim<ndims;dim++) dest_rank_coords[dim]=std::min(dest_rank_coords[dim],nranks_per_dir[dim]-1);
  coords_t dest_rank_origin=glb_coords_of_origin_of_rank_coords(dest_rank_coords);
  
  //wrap coords to ask mpi what's the rank
  int dest_rank_coords_wrap[ndims];
  for(size_t dim=0;dim<ndims;dim++) dest_rank_coords_wrap[dim]=dest_rank_coords[dim];
  MPI_Cart_rank(cart_comm,dest_rank_coords_wrap,&dest_rank);
  
  //get coords in the dest rank, and its size
  coords_t dest_loc_coords=dest_glb_coords-dest_rank_origin;
  coords_t dest_loc_sizes=loc_sizes_of_rank_coords(dest_rank_coords);
  
  //get site in the dest rank
  dest_site=site_of_coords(dest_loc_coords,dest_loc_sizes);
}

//compute internal volume
int geometry_t::bulk_volume(coords_t L)
{
  int intvol=1;
  size_t dim=0;
    do
      {
        if(L[dim]>2) intvol*=L[dim]-2;
        else intvol=0;
        
        dim++;
      }
    while(intvol!=0 && dim<ndims);
    
    return intvol;
}

//compute the bulk volume of the local lattice, given by L/R
int geometry_t::bulk_recip_lat_volume(coords_t R,coords_t L)
{
  coords_t X=L/R;
  return bulk_volume(X);
}

//compute the variance of the border
int geometry_t::compute_border_variance(coords_t L,coords_t P,int factorize_processor)
{
  int S2B=0,SB=0;
  for(size_t ib=0;ib<ndims;ib++)
    {
      int B=1;
      for(size_t dim=0;dim<ndims;dim++) if(dim!=ib) B*=(factorize_processor) ? L[dim]/P[dim] : P[dim];
      SB+=B;
      S2B+=B*B;
    }
  SB/=ndims;
  S2B/=ndims;
  S2B-=SB*SB;
    
  return S2B;
}

//get coords of site
coords_t geometry_t::coords_of_site(int site,coords_t sizes)
{
  coords_t coords(ndims);
  size_t dim=ndims;
  
  do
    {
      dim--;
      coords[dim]=site%sizes[dim];
      site/=sizes[dim];
    }
  while(dim!=0);
  
  return coords;
}

//return site of coords
int geometry_t::site_of_coords(coords_t coords,coords_t sizes)
{
  int site=0;
  for(size_t dim=0;dim<ndims;dim++) site=site*sizes[dim]+coords[dim];
  return site;
}

//divide the global lattice into ranks in homogeneous way
void geometry_t::partition_ranks_homogeneously()
{
  bool something_found=false;
  
  int NR=simul->nranks;
  int LV=nglb_sites/NR;
  
  if(NR==1) for(size_t dim=0;dim<ndims;dim++) nranks_per_dir[dim]=1;
  else
    {
      //minimal variance border
      int mBV=-1;
      
      //factorize the volume
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
      
      //compute the number of combinations: this is given by ndims^nfact
      int ncombo=1;
      for(int ifact=0;ifact<nfact;ifact++) ncombo*=ndims;
      
      //find the partition which minimize the surface and the surface variance
      int min_surf_LV=-1;
      int icombo=0;
      for(size_t dim=0;dim<ndims;dim++) nranks_per_dir[dim]=-1;
      
      do
	{
	  //number of ranks in each direction for current partitioning
	  coords_t R(ndims);
	  for(size_t dim=0;dim<ndims;dim++) R[dim]=1;
	  
	  //find the partioning corresponding to icombo
	  int ifact=nfact-1;
	  int valid_partitioning=1;
	  int jcombo=icombo;
	  do
	    {
	      //find the direction: this is given by the ifact digit of icombo wrote in base ndims
	      size_t dim=(jcombo&(ndims-1));
	      jcombo/=ndims;
	      
	      MASTER_PRINTF("ifact: %d, dim: %d\n",ifact,dim);
	      
	      //if we are factorizing local lattice, rank factor is given by list_fact, otherwise L/list_fact
	      R[dim]*=list_fact[ifact];
	      
	      //check that the total volume glb_sizes[dim] is a multiple and it is larger than the number of proc
	      valid_partitioning=(glb_sizes[dim]%R[dim]==0 && glb_sizes[dim]>=R[dim]);
	      if(valid_partitioning) ifact--;
	    }
	  while(valid_partitioning && ifact>=0);
	  
	  //validity could have changed
	  if(valid_partitioning)
	    {
	      //if we are factorizing reciprocal lattice, convert back to rank grid
	      for(size_t dim=0;dim<ndims;dim++)
		if(!factorize_rank)
		  R[dim]=glb_sizes[dim]/R[dim];
	      
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
		  for(size_t dim=0;dim<ndims;dim++) nranks_per_dir[dim]=R[dim];
		  something_found=1;
		}
	      
	      icombo++;
	    }
	  //skip all remaining factorization using the same structure
	  else
	    {
	      int skip=1;
	      for(int jfact=0;jfact<ifact-1;jfact++) skip*=ndims;
	      icombo+=skip;
	    }
	}
      while(icombo<ncombo);
      
      if(!something_found) CRASH_SOFTLY("no valid partitioning found");
    }
}

//do not rely on homogeneous local sizes
void geometry_t::partition_ranks_inhomogeneously()
{
  //factorize the number of ranks
  int list_fact_nranks[log2N(simul->nranks)];
  int nfact_nranks=factorize(list_fact_nranks,simul->nranks);
  
  //compute the number of combinations: this is given by ndims^nfact
  int ncombo=1;
  for(int ifact=0;ifact<nfact_nranks;ifact++) ncombo*=ndims;
  
  //find the partition which minimize the volume variance
  int min_vol_var=-1;
  int icombo=0;
  for(size_t dim=0;dim<ndims;dim++) nranks_per_dir[dim]=-1;
  
  do
    {
      //number of ranks in each direction for current partitioning
      coords_t nranks_per_dir_old=nranks_per_dir;
      for(size_t dim=0;dim<ndims;dim++) nranks_per_dir[dim]=1;
      
      //find the partioning corresponding to icombo
      int jcombo=icombo;
      for(int ifact=nfact_nranks-1;ifact>=0;ifact--)
	{
	  //find the direction: this is given by the ifact digit of icombo wrote in base ndims
	  size_t dim=(jcombo&(ndims-1));
	  jcombo/=ndims;
	  
	  //rank factor is given by list_fact
	  nranks_per_dir[dim]*=list_fact_nranks[ifact];
	}
      
      //find common sizes
      for(size_t dim=0;dim<ndims;dim++) loc_comm_sizes[dim]=(int)((double)glb_sizes[dim]/nranks_per_dir[dim]+0.5);
      
      //compute variance
      int v_sum=0,v_sum2=0;
      for(int irank=0;irank<simul->nranks;irank++)
	{
	  coords_t rank_coords=coords_of_site(irank,nranks_per_dir);
	  coords_t rank_sizes=loc_sizes_of_rank_coords(rank_coords);
	  int v=rank_sizes.total_product();
	  v_sum+=v;
	  v_sum+=v*v;
	  
	  //MASTER_PRINTF("  icombo %d, rank %d:\t%d",icombo,irank,rank_sizes[0]);
	  //for(size_t dim=1;dim<ndims;dim++) MASTER_PRINTF("x%d",rank_sizes[dim]);
	  //MASTER_PRINTF(" vol: %d\n",v);
	}
      v_sum/=simul->nranks;
      v_sum2/=simul->nranks;
      v_sum2-=v_sum*v_sum;
      
      //if new minimal volume variance
      int new_minimal_vol_var=0;
      //if it is the minimal volume variance (or first valid combo) copy it and compute the border size
      if(v_sum2<min_vol_var||min_vol_var==-1)
	{
	  new_minimal_vol_var=1;
	  min_vol_var=v_sum2;
	}
      
      //convert to previous one if it not a new minimal
      if(!new_minimal_vol_var) nranks_per_dir=nranks_per_dir_old;
      
      icombo++;
    }
  while(icombo<ncombo);
}

//initialize the geometry
void geometry_t::start(coords_t ext_glb_sizes)
{
  GET_THREAD_ID();
  
  THREAD_BARRIER();
  if(IS_MASTER_THREAD)
    {
      //resize all the coords
      nranks_per_dir.resize(ndims);
      rank_coords.resize(ndims);
      glb_sizes.resize(ndims);
      loc_comm_sizes.resize(ndims);
      loc_sizes.resize(ndims);
      glb_coords_of_loc_origin.resize(ndims);
      
      //copy glb sizes and compute total volume
      nglb_sites=1;
      glb_sizes=ext_glb_sizes;
      nglb_sites=glb_sizes.total_product();
      
      //set not having inited loc_rnd_gen
      loc_rnd_gen_inited=false;
      
      //output global size and number of ranks
      MASTER_PRINTF("Global lattice grid:\t%d",glb_sizes[0]);
      for(size_t dim=1;dim<ndims;dim++) MASTER_PRINTF("x%d",glb_sizes[dim]);
      MASTER_PRINTF(" = %d\n",nglb_sites);
      MASTER_PRINTF("Number of running ranks: %d\n",simul->nranks);
      
      //check that we are not using more ranks than sites
      if(simul->nranks>nglb_sites) CRASH_SOFTLY("too many ranks!");
      
      //check wheter each ranks will have the same number of sites
      homogeneous_partitioning=(nglb_sites%simul->nranks==0);
      
      //on the base of the results, partition
      if(homogeneous_partitioning) partition_ranks_homogeneously();
      else partition_ranks_inhomogeneously();
      MASTER_PRINTF("Creating ranks grid:\t%d",nranks_per_dir[0]);
      for(size_t dim=1;dim<ndims;dim++) MASTER_PRINTF("x%d",nranks_per_dir[dim]);
      MASTER_PRINTF("\n");
      
      //init the mpi communicator
      int periods[ndims];
      int nranks_per_dir_conv[ndims];
      for(size_t dim=0;dim<ndims;dim++)
	{
	  periods[dim]=1;
	  nranks_per_dir_conv[dim]=nranks_per_dir[dim];
	}
      MPI_Cart_create(MPI_COMM_WORLD,ndims,nranks_per_dir_conv,periods,1,&cart_comm);
      
      //takes rank and coord of local rank
      MPI_Comm_rank(cart_comm,&cart_rank);
      int rank_coords_conv[ndims];
      MPI_Cart_coords(cart_comm,cart_rank,ndims,rank_coords_conv);
      for(size_t dim=0;dim<ndims;dim++) rank_coords[dim]=rank_coords_conv[dim];
      
      //understand the local sizes
      nloc_sites=1;
      for(size_t dim=0;dim<ndims;dim++) loc_comm_sizes[dim]=(int)((double)glb_sizes[dim]/nranks_per_dir[dim]+0.5);
      glb_coords_of_loc_origin=glb_coords_of_origin_of_rank_coords(rank_coords);
      loc_sizes=loc_sizes_of_rank_coords(rank_coords);
      nloc_sites*=loc_sizes.total_product();
      
      //print info on 
      if(homogeneous_partitioning) MASTER_PRINTF("Local volume for each node: %d",loc_sizes[0]);
      else                         MASTER_PRINTF("Common nodes local volume: %d",loc_sizes[0]);
      for(size_t dim=1;dim<ndims;dim++) MASTER_PRINTF("x%d ",loc_sizes[dim]);
      MASTER_PRINTF("\n");
      
      //write non-homogeneous partitioning
      if(!homogeneous_partitioning)
	{
	  for(int irank=0;irank<simul->nranks;irank++)
	    {
	      if(simul->rank_id==irank)
		if(loc_sizes!=loc_comm_sizes)
		  {
		    printf(" inhomogeneous rank %d: %d",irank,loc_sizes[0]);
		    for(size_t dim=1;dim<ndims;dim++) printf("x%d ",loc_sizes[dim]);
		    printf("= %d\n",loc_sizes.total_product());
		  }
	      MPI_Barrier(MPI_COMM_WORLD);
	    }
	}
      
      //fill lookup tables
      for(int loc_site=0;loc_site<nloc_sites;loc_site++)
	{
	  loc_coords_of_loc_site_table.push_back(coords_of_site(loc_site,loc_sizes));
	  glb_coords_of_loc_site_table.push_back(loc_coords_of_loc_site(loc_site)+glb_coords_of_loc_origin);
	  glb_site_of_loc_site_table.push_back(glb_site_of_loc_coords(loc_coords_of_loc_site_table[loc_site]));
	}
    }
  THREAD_BARRIER();
  
  //start no-neighbors
  per_site_neighs_t no_neighbors_per_site;
  no_neighbors=NEW_BLOCKING("no_neighs") neighs_t(this,&no_neighbors_per_site);
  
  //start per-site first neighbors
  per_site_neighs_t first_neighbors_per_site;
  coords_t site(ndims);
  for(size_t dim=0;dim<ndims;dim++) site[dim]=0;
  for(size_t dim=0;dim<ndims;dim++) //dimension loop
    for(int du=-1;du<=1;du+=2) //down then up
      {
	site[dim]=du;
	first_neighbors_per_site.add_neighbor(site);
	site[dim]=0;
      }
  
  //initialize first neighbors
  first_neighbors=NEW_BLOCKING("first_neighs") neighs_t(this,&first_neighbors_per_site);
}

//destructor
geometry_t::~geometry_t()
{
  DELETE_NON_BLOCKING(no_neighbors);
  DELETE_NON_BLOCKING(first_neighbors);
}

//print all the infos
void geometry_t::print()
{
  printf("Ndims: %d\n",(int)ndims);
  printf("Cart_rank: %d\n",cart_rank);
  printf("Nranks_per_dir: %d",nranks_per_dir[0]);
  for(size_t dim=1;dim<ndims;dim++) printf("x%d ",nranks_per_dir[dim]);
  printf("\n");
  printf("Rank_coords: %d",rank_coords[0]);
  for(size_t dim=1;dim<ndims;dim++) printf("x%d ",rank_coords[dim]);
  printf("\n");
  printf("Homogeneous_partitioning: %d\n",homogeneous_partitioning);
  printf("Glb_sizes: %d",glb_sizes[0]);
  for(size_t dim=1;dim<ndims;dim++) printf("x%d ",glb_sizes[dim]);
  printf("\n");
  printf("Loc_sizes: %d",loc_sizes[0]);
  for(size_t dim=1;dim<ndims;dim++) printf("x%d ",loc_sizes[dim]);
  printf("\n");
  printf("Loc_comm_sizes: %d",loc_comm_sizes[0]);
  for(size_t dim=1;dim<ndims;dim++) printf("x%d ",loc_comm_sizes[dim]);
  printf("\n");
  printf("Nglb_sites: %d\n",nglb_sites);
  printf("Nloc_sites: %d\n",nloc_sites);
  printf("Glb_coords_of_loc_origin: %d",glb_coords_of_loc_origin[0]);
  for(size_t dim=1;dim<ndims;dim++) printf("x%d ",glb_coords_of_loc_origin[dim]);
  printf("\n");
}

//initialize the local number generators
void geometry_t::init_loc_rnd_gen()
{
  GET_THREAD_ID();
  
  if(IS_MASTER_THREAD)
    {
      //resize to proper volume
      loc_rnd_gen.resize(nloc_sites);
      
      //get the seed for the whole grid
      int seed=simul->glb_rnd_gen.get_unif(0,RAND_MAX);
      
      //start the grid
      for(int site=0;site<nloc_sites;site++)
	loc_rnd_gen[site].init(seed+glb_site_of_loc_site(site));
    }
  THREAD_BARRIER();
}
