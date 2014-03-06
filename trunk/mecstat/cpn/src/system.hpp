#ifndef _SYSTEM_HPP
#define _SYSTEM_HPP

#include <mpi.h>
#include <mpi.h>

#include "thread_macros.hpp"

extern bool system_started;

//system
struct system_t
{
  int rank,nranks;
  
  //basic MPI initialization
  void init_MPI_thread(int narg,char **arg);
  void get_MPI_nranks(){MPI_Comm_size(MPI_COMM_WORLD,&nranks);} //get nranks
  void get_MPI_rank(){MPI_Comm_rank(MPI_COMM_WORLD,&rank);} //get rank
  
  //debug
  void print_backtrace_list();
  void abort(int arg);
  
  system_t(int narg,char **arg);
  ~system_t();
private:
  system_t();
};

extern system_t *system;

#endif
