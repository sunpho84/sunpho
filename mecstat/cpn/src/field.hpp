#ifndef _FIELDS_HPP
#define _FIELDS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>

#include "geometry.hpp"
#include "reduce.hpp"
#include "simul.hpp"

#define PARALLEL_FOR_SITES_OF_FIELD(SITE,FIELD) PARALLEL_FOR(SITE,0,FIELD.neighs_ptr->geometry->nloc_sites)

template <class T> class field_t
{
public:
  T *data;                                     //store the internal data
  neighs_t *neighs_ptr;                        //neighboring connections
  bool outer_sites_synced;                     //store the info on whether outer sites are synced
  
  field_t():data(NULL) {}
  field_t(const char *name,neighs_t *ext_neighs_ptr){create(name,ext_neighs_ptr);}
  void create(const char *name,neighs_t *ext_neighs_ptr);
  ~field_t(){destroy();}
  void destroy();
  
  void operator=(const field_t<T> &in)
  {
    GET_THREAD_ID();
    if(this!=&in) PARALLEL_FOR_SITES_OF_FIELD(site,in) (*this)[site]=in[site];
    THREAD_BARRIER();
  }
  
  int get_neigh(int i,int dir){return (*neighs_ptr)[i][dir];}   //return directly the neighbor
  
  T &operator[](int iel){return data[iel];}                //accede to internal
  const T &operator[](int iel) const{return data[iel];}    //accede to internal with const

  void sync_outer_sites();                                 //sync outer sites
  void mark_touched(bool isb=true)                         //mark as touched
  {if(isb)THREAD_BARRIER();outer_sites_synced=false;if(isb)THREAD_BARRIER();}
  void mark_touched_nonblocking(){mark_touched(false);}
  
  T reduce();          //perform the total summ
  double get_norm2();  //return the "norm2"
  void normalize();    //normalize
private:
  field_t(field_t<T> &in) {CRASH_SOFTLY("copy constructor invocated");}
};

//initialize
template<class T> void field_t<T>::create(const char *name,neighs_t *ext_neighs_ptr)
{
  //mark unsynced external sites
  outer_sites_synced=false;

  //allocate the buffer
  neighs_ptr=ext_neighs_ptr;
  data=NEW_ARRAY_BLOCKING(name,T,neighs_ptr->ntotal_sites);
  
  //check that communication buffer is large enough to allow borders comm
  size_t needed_size=sizeof(T)*neighs_ptr->nouter_sites;
  if(simul->comm_buff_size<needed_size)
    {
      MASTER_PRINTF("simul->comm_buff_size: %d, required: %d\n",(int)simul->comm_buff_size,(int)needed_size);
      if(simul->comm_buff_size!=0)
	{
	  DELETE_BLOCKING(simul->comm_buff);
	  MASTER_PRINTF("deleting old comm_buff\n");
	}
      simul->comm_buff=NEW_ARRAY_BLOCKING("comm_buff",char,needed_size);
      simul->comm_buff_size=needed_size;
    }
}

//gather sites from neighbors, send them
template <class T> void field_t<T>::sync_outer_sites()
{
  GET_THREAD_ID();

  //communicate only if needed
  if(outer_sites_synced==false)
    {
      //prepare communication buffer
      PARALLEL_FOR(site,0,neighs_ptr->nsites_to_send)
	((T*)simul->comm_buff)[site]=data[neighs_ptr->list_sending[site]];
      THREAD_BARRIER();
      
      //good moment to mark as synced
      outer_sites_synced=true;
      
      //send and receive
      if(IS_MASTER_THREAD)
	{
	  MPI_Request requests[neighs_ptr->nranks_to_ask+neighs_ptr->nranks_asking];
	  for(int irank=0;irank<neighs_ptr->nranks_to_ask;irank++) //data is of type T
	    MPI_Irecv(data+neighs_ptr->geometry->nloc_sites+neighs_ptr->ranks_to_ask[irank].dest,
		      neighs_ptr->ranks_to_ask[irank].size*sizeof(T),MPI_CHAR,neighs_ptr->ranks_to_ask[irank].rank,
		      neighs_ptr->geometry->cart_rank,neighs_ptr->geometry->cart_comm,requests+irank);
	  for(int irank=0;irank<neighs_ptr->nranks_asking;irank++) //comm_buff is of type char
	    MPI_Isend(simul->comm_buff+neighs_ptr->ranks_asking[irank].dest*sizeof(T),
		      neighs_ptr->ranks_asking[irank].size*sizeof(T),
		      MPI_CHAR,neighs_ptr->ranks_asking[irank].rank,neighs_ptr->ranks_asking[irank].rank,
		      neighs_ptr->geometry->cart_comm,requests+neighs_ptr->nranks_to_ask+irank);
	  MPI_Waitall(neighs_ptr->nranks_asking+neighs_ptr->nranks_to_ask,requests,MPI_STATUS_IGNORE);
	}
      THREAD_BARRIER();
    }
}

//destroy the field
template <class T> void field_t<T>::destroy()
{
  if(data) DELETE_BLOCKING(data);
}

//reduce across all the nodes
template <class T> T field_t<T>::reduce()
{
  GET_THREAD_ID();
  
  //reset local residue
  T loc_thread_res=0;
  
  //sum local thread contribution
  PARALLEL_FOR(site,0,neighs_ptr->geometry->nloc_sites)
    loc_thread_res+=data[site];
  
  return rank_threads_reduce(loc_thread_res);
}

//normalize all sites
template <class T> void field_t<T>::normalize()
{
  GET_THREAD_ID();  
  PARALLEL_FOR(site,0,neighs_ptr->geometry->nloc_sites) data[site].normalize();
  THREAD_BARRIER();
}

//get the squared norm
template <class T> double field_t<T>::get_norm2()
{
  GET_THREAD_ID();
  
  //reset local residue
  double loc_thread_res=0;
  
  //sum local thread contribution
  PARALLEL_FOR(site,0,neighs_ptr->geometry->nloc_sites)
    loc_thread_res+=data[site].get_norm2();
  return rank_threads_reduce(loc_thread_res);
}

#endif
