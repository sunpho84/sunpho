#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <list>

#include "debug.hpp"
#include "geometry.hpp"
#include "simul.hpp"
#include "threads.hpp"
#include "neighs.hpp"

//hold temporarily info
typedef std::map<int,int> sites_to_send_t;

//mark the connections of the site
void neighs_t::mark_all_neighbors(int loc_site,coords_t glb_site_coords,per_site_neighs_t *per_site_neighs,
					site_list_per_rank_t &outer_sites_per_rank,bool recursive)
{
  //find glb_coords of neighbors
  for(size_t inei=0;inei<nneighs_per_site;inei++)
    {
      //convert to rank and site
      int dest_rank,dest_site;
      coords_t neigh_coords=glb_site_coords+(*per_site_neighs)[inei];
      geometry->rank_and_loc_site_of_rel_glb_coords(dest_rank,dest_site,neigh_coords);
      
      //if local, point directly to it
      if(dest_rank==geometry->cart_rank) neighs[loc_site*nneighs_per_site+inei]=dest_site;
      else
	{
	  //search rank list, if present continue
	  site_list_per_rank_t::iterator true_dest_rank_ptr=outer_sites_per_rank.find(dest_rank);
	  if(true_dest_rank_ptr!=outer_sites_per_rank.end())
	    {
	      //search dest site id, if present continue
	      std::map<int,int>::iterator true_dest_site_ptr=true_dest_rank_ptr->second.find(dest_site);
	      if(true_dest_site_ptr!=true_dest_rank_ptr->second.end())
		{
		  //mark it
		  int true_dest_site=true_dest_site_ptr->second;
		  neighs[loc_site*nneighs_per_site+inei]=true_dest_site;
		  
		  //connect again
		  if(recursive) mark_all_neighbors(true_dest_site,neigh_coords,per_site_neighs,outer_sites_per_rank,false);
		}
	      else neighs[loc_site*nneighs_per_site+inei]=-1;
	    }
	  else neighs[loc_site*nneighs_per_site+inei]=-1;
	}
    }
}

//construct using per-site mask
neighs_t::neighs_t(geometry_t *geometry,per_site_neighs_t *per_site_neighs) : geometry(geometry),nneighs_per_site(per_site_neighs->size())
{
  GET_THREAD_ID();
  
  THREAD_BARRIER();
  if(IS_MASTER_THREAD)
    {
      //loop on all local sites and mark all required neighbours
      site_list_per_rank_t outer_sites_per_rank;
      for(int loc_site=0;loc_site<geometry->nloc_sites;loc_site++)
	{
	  //find glb_coords of neighbors
	  coords_t glb_site_coords=geometry->glb_coords_of_loc_site(loc_site);
	  
	  for(per_site_neighs_t::iterator it=per_site_neighs->begin();it!=per_site_neighs->end();it++)
	    {
	      //convert to rank and site
	      int dest_rank,dest_site;
	      geometry->rank_and_loc_site_of_rel_glb_coords(dest_rank,dest_site,glb_site_coords+(*it));
	      
	      //add it IF not local
	      if(dest_rank!=geometry->cart_rank) outer_sites_per_rank[dest_rank][dest_site]=0;
	    }
	}
      
      //count outer sites, so indexing them
      nouter_sites=0;
      ntotal_sites=geometry->nloc_sites;
      for(site_list_per_rank_t::iterator outer_sites=outer_sites_per_rank.begin();
	  outer_sites!=outer_sites_per_rank.end();outer_sites++)
	for(std::map<int,int>::iterator site=outer_sites->second.begin();site!=outer_sites->second.end();site++)
	  {
	    site->second=ntotal_sites;
	    ntotal_sites++;
	    nouter_sites++;
	  }
      
      //allocate neighbors and connect all sites
      neighs=NEW_NON_BLOCKING("neighs") int[ntotal_sites*nneighs_per_site];
      for(int loc_site=0;loc_site<geometry->nloc_sites;loc_site++)
	mark_all_neighbors(loc_site,geometry->glb_coords_of_loc_site(loc_site),per_site_neighs,outer_sites_per_rank,true);
      
      //create the list of ranks to ask to
      nranks_to_ask=outer_sites_per_rank.size();
      ranks_to_ask=NEW_NON_BLOCKING("ranks_to_ask") rank_to_ask_t[nranks_to_ask];
      site_list_per_rank_t::iterator sites_per_rank_map=outer_sites_per_rank.begin();
      for(int irank=0;irank<nranks_to_ask;irank++)
	{
	  ranks_to_ask[irank].rank=sites_per_rank_map->first;
	  ranks_to_ask[irank].size=sites_per_rank_map->second.size();
	  if(irank==0) ranks_to_ask[irank].dest=0;
	  else         ranks_to_ask[irank].dest=ranks_to_ask[irank-1].dest+ranks_to_ask[irank-1].size;
	  
	  sites_per_rank_map++;
	}
      
      //communicate to each rank how many elements to ask
      sites_to_send_t *temp_nsites_to_send_to=new sites_to_send_t;
      for(int delta_rank=1;delta_rank<simul->nranks;delta_rank++)
	{
	  int dest_rank=(geometry->cart_rank+simul->nranks+delta_rank)%simul->nranks;
	  int recv_rank=(geometry->cart_rank+simul->nranks-delta_rank)%simul->nranks;
	  
	  int nto_ask=0,nasking;
	  int jrank=0;
	  if(nranks_to_ask>0)
	    do
	      {
		if(ranks_to_ask[jrank].rank==dest_rank) nto_ask=ranks_to_ask[jrank].size;
		jrank++;
	      }
	    while(nto_ask==0 && jrank<nranks_to_ask);
	  
	  //transfer the info and mark
	  MPI_Sendrecv(&nto_ask,1,MPI_INT,dest_rank,0, &nasking,1,MPI_INT,recv_rank,0,  MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  //MASTER_PRINTF("Told to rank %d to send %d, told by rank %d to send %d\n",dest_rank,nto_ask,recv_rank,nasking);
	  if(nasking!=0) (*temp_nsites_to_send_to)[recv_rank]=nasking;
	}
      
      //convert to store this into proper list
      nsites_to_send=0;
      nranks_asking=temp_nsites_to_send_to->size();
      ranks_asking=NEW_NON_BLOCKING("ranks_asking") rank_asking_t[nranks_asking];
      std::map<int,int>::iterator nsites_to_send_to_ptr=temp_nsites_to_send_to->begin();
      for(int irank=0;irank<nranks_asking;irank++)
	{
	  ranks_asking[irank].rank=nsites_to_send_to_ptr->first;
	  ranks_asking[irank].size=nsites_to_send_to_ptr->second;
	  if(irank==0) ranks_asking[irank].dest=0;
	  else ranks_asking[irank].dest=ranks_asking[irank-1].dest+ranks_asking[irank-1].size;
	  nsites_to_send+=ranks_asking[irank].size;
	  //SHOUT("%d/%d, %d size %d dest %d",irank,nranks_asking,ranks_asking[irank].rank,
	  //ranks_asking[irank].size,ranks_asking[irank].dest);
	  
	  nsites_to_send_to_ptr++;
	}
      
      //allocate the list of sites to send, and receive it
      list_sending=NEW_NON_BLOCKING("list_sending") int[nsites_to_send];
      MPI_Request requests[nranks_asking+nranks_to_ask];
      for(int irank=0;irank<nranks_asking;irank++)
	MPI_Irecv(list_sending+ranks_asking[irank].dest,ranks_asking[irank].size,MPI_INT,
		  ranks_asking[irank].rank,ranks_asking[irank].rank,geometry->cart_comm,requests+irank);

      //send the lists
      sites_per_rank_map=outer_sites_per_rank.begin();
      int *temp_list[nranks_to_ask];
      for(int irank=0;irank<nranks_to_ask;irank++)
	{
	  //prepare the list
	  temp_list[irank]=NEW_NON_BLOCKING("temp_list") int[ranks_to_ask[irank].size];
	  std::map<int,int>::iterator sites_list=sites_per_rank_map->second.begin();
	  for(int site=0;site<ranks_to_ask[irank].size;site++)
	    {
	      temp_list[irank][site]=sites_list->first;
	      //MASTER_PRINTF("Asking rank %d site [%d/%d] %d\n",ranks_to_ask[irank].rank,site,sites_per_rank_map->second.size(),temp_list[irank][site]);
	      sites_list++;
	    }
	  sites_per_rank_map++;
	  
	  //send it
	  MPI_Isend(temp_list[irank],ranks_to_ask[irank].size,MPI_INT,
		    ranks_to_ask[irank].rank,geometry->cart_rank,geometry->cart_comm,requests+nranks_asking+irank);
	}
      
      //wait and free temp list
      MPI_Waitall(nranks_asking+nranks_to_ask,requests,MPI_STATUS_IGNORE);
      for(int irank=0;irank<nranks_to_ask;irank++) DELETE_NON_BLOCKING(temp_list[irank]);
      
      //check
      {
	int *recv=NEW_NON_BLOCKING("recv") int[nouter_sites];
	int *send=NEW_NON_BLOCKING("send") int[nsites_to_send];
	for(size_t site=0;site<nsites_to_send;site++)
	  send[site]=list_sending[site];
	MPI_Request requests[nranks_to_ask+nranks_asking];
	for(int irank=0;irank<nranks_to_ask;irank++)
	  MPI_Irecv(recv+ranks_to_ask[irank].dest,ranks_to_ask[irank].size,MPI_INT,
		    ranks_to_ask[irank].rank,geometry->cart_rank,geometry->cart_comm,requests+irank);
	for(int irank=0;irank<nranks_asking;irank++)
	  MPI_Isend(send+ranks_asking[irank].dest,ranks_asking[irank].size,MPI_INT,
		    ranks_asking[irank].rank,ranks_asking[irank].rank,geometry->cart_comm,requests+nranks_to_ask+irank);
	MPI_Waitall(nranks_asking+nranks_to_ask,requests,MPI_STATUS_IGNORE);

	//final check
	sites_per_rank_map=outer_sites_per_rank.begin();
	int isite=0;
	for(int irank=0;irank<nranks_to_ask;irank++)
	  {
	    //prepare the list
	    std::map<int,int>::iterator sites_list=sites_per_rank_map->second.begin();
	    for(int site=0;site<ranks_to_ask[irank].size;site++)
	      {
		//MASTER_PRINTF("rank %d site[%d] to ask: %d, obtained %d\n",irank,site,sites_list->first,recv[isite]);
		if(sites_list->first!=recv[isite])
		  CRASH_SOFTLY("site %d from rank %d failed: %d!=%d",site,ranks_to_ask[irank].rank,
			       sites_list->first,recv[isite]);
		sites_list++;
		isite++;
	      }
	    sites_per_rank_map++;
	  }
	
	DELETE_NON_BLOCKING(recv);
	DELETE_NON_BLOCKING(send);
      }
      
      //for(int irank=0;irank<nranks_asking;irank++)
      //MASTER_PRINTF("rank %d dest: %d\n",irank,ranks_asking[irank].dest);
      
      //int site_i=0;
      //for(int irank=0;irank<nranks_asking;irank++)
      //for(int site=0;site<ranks_asking[irank].size;site++)
      //{
      //MASTER_PRINTF("site %d sending: %d/%d to %d[%d]\n",site_i,list_sending[site_i],geometry->nloc_sites,
      //ranks_asking[irank].rank,irank);
      //site_i++;
      //}
      
      //remove temporary list of sites to send
      delete temp_nsites_to_send_to;
    }
  THREAD_BARRIER();
}

//deallocate
neighs_t::~neighs_t()
{
  DELETE_NON_BLOCKING(neighs);
  DELETE_NON_BLOCKING(ranks_to_ask);
  DELETE_NON_BLOCKING(ranks_asking);
  DELETE_NON_BLOCKING(list_sending);
}
