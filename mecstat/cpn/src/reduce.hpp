#ifndef _REDUCE_HPP
#define _REDUCE_HPP

#include "threads.hpp"

//reduce a number across all the threads
template <class T> T threads_reduce(T &in)
{
  GET_THREAD_ID();
  
  //upload the address to look at
  THREAD_BARRIER();
  simul->thread_res_arr[thread_id]=(char*)&in;
  THREAD_BARRIER();
  T *out=(T*)(simul->thread_res_arr[0]);
  
  //reduce across all the threads
  if(IS_MASTER_THREAD)
    for(int thread_jd=1;thread_jd<simul->nthreads;thread_jd++)
      (*out)+=(*((T*)(simul->thread_res_arr[thread_jd])));
  THREAD_BARRIER();
    
  return *out;
}

//find the correct datatype
inline MPI_Datatype find_MPI_type(int *in){return MPI_INT;}
inline MPI_Datatype find_MPI_type(double *in){return MPI_DOUBLE;}

//reduce a number across all the nodes
template <class T> T rank_threads_reduce(T &in)
{
  GET_THREAD_ID();
  
  //upload the address to look at
  THREAD_BARRIER();
  simul->thread_res_arr[thread_id]=(char*)&in;
  THREAD_BARRIER();
  T *out=(T*)(simul->thread_res_arr[0]);
  
  if(IS_MASTER_THREAD)
    {
      //reduce across all the threads
      T temp=0;
      for(int thread_jd=0;thread_jd<simul->nthreads;thread_jd++)
	temp+=(*((T*)(simul->thread_res_arr[thread_jd])));
      
      //reduce across all the nodes
      MPI_Allreduce(&temp,out,1,find_MPI_type(&in),MPI_SUM,MPI_COMM_WORLD);
    }
  THREAD_BARRIER();
    
  return *out;
}

#endif
