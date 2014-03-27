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
  T *data;              //store the internal data
  field_t(const char *name,neighs_t *neighs);  //constructor
  ~field_t();
  
  int get_neigh(int i,int dir){return (*neighs_ptr)[i][dir];}   //return directly the neighbor
  neighs_t *neighs_ptr;                                         //neighboring connections
  
  T &operator[](int iel){return data[iel];}    //accede to internal
  void sync_outer_sites();                     //sync outer sites
  T reduce();
private:
  field_t();
};

//constructor
template <class T> field_t<T>::field_t(const char *name,neighs_t *neighs_ptr) : neighs_ptr(neighs_ptr)
{
  //allocate the buffer
  data=NEW_BLOCKING(name) T[neighs_ptr->ntotal_sites];
  
  //check that communication buffer is large enough to allow borders comm
  size_t needed_size=sizeof(T)*neighs_ptr->nouter_sites;
  if(simul->comm_buff_size<needed_size)
    {
      if(simul->comm_buff_size!=0) DELETE_BLOCKING(simul->comm_buff);
      simul->comm_buff=NEW_BLOCKING("comm_buff") char[needed_size];
    }
}

template <class T> void field_t<T>::sync_outer_sites()
{
  GET_THREAD_ID();
  
  //prepare communication buffer
  PARALLEL_FOR(site,0,neighs_ptr->nsites_to_send)
    ((T*)simul->comm_buff)[site]=data[neighs_ptr->list_sending[site]];
  THREAD_BARRIER();
  
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

//destroy the field
template <class T> field_t<T>::~field_t()
{
  DELETE_BLOCKING(data);
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

#endif
