#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <list>

#include "geometry.hpp"
#include "simul.hpp"
#include "total_neighs.hpp"

//mark the connections of the site
void total_neighs_t::mark_all_neighbors(int loc_site,coords_t glb_site_coords,per_site_neighs_t &per_site_neighs,
					site_list_per_rank_t &outer_sites_per_rank,bool recursive)
{
  //find glb_coords of neighbors
  for(size_t inei=0;inei<nneighs_per_site;inei++)
    {
      //convert to rank and site
      int dest_rank,dest_site;
      coords_t neigh_coords=glb_site_coords+per_site_neighs[inei];
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
total_neighs_t::total_neighs_t(geometry_t *geometry,per_site_neighs_t &per_site_neighs) : geometry(geometry),nneighs_per_site(per_site_neighs.size())
{
  //loop on all local sites and mark all required neighbours
  site_list_per_rank_t outer_sites_per_rank;
  for(int loc_site=0;loc_site<geometry->nloc_sites;loc_site++)
    {
      //find glb_coords of neighbors
      coords_t glb_site_coords=geometry->glb_coords_of_loc_site(loc_site);
      for(per_site_neighs_t::iterator it=per_site_neighs.begin();it!=per_site_neighs.end();it++)
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
  neighs=ALLOCATE("neighs",ntotal_sites*nneighs_per_site,int);
  for(int loc_site=0;loc_site<geometry->nloc_sites;loc_site++)
    mark_all_neighbors(loc_site,geometry->glb_coords_of_loc_site(loc_site),per_site_neighs,outer_sites_per_rank,true);
  
  //create the list of ranks to ask to
  nranks_to_ask=outer_sites_per_rank.size();
  ranks_to_ask=ALLOCATE("ranks_to_ask",nranks_to_ask,rank_id_size_dest_start_t);
  site_list_per_rank_t::iterator sites_per_rank=outer_sites_per_rank.begin();
  for(int irank=0;irank<nranks_to_ask;irank++)
    {
      ranks_to_ask[irank].rank=sites_per_rank->first;
      ranks_to_ask[irank].size=sites_per_rank->second.size();
      ranks_to_ask[irank].dest=(irank==0)?geometry->nloc_sites:ranks_to_ask[irank-1].dest+ranks_to_ask[irank].size;
      
      sites_per_rank++;
    }
  
  
}

//deallocate
total_neighs_t::~total_neighs_t()
{
  FREE(neighs);
  FREE(ranks_to_ask);
}
