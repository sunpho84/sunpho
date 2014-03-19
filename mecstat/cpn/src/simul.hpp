#ifndef _SIMUL_HPP
#define _SIMUL_HPP

#include <mpi.h>
#include <mpi.h>

#include "threads.hpp"
#include "vectors.hpp"

extern bool simul_started;

#define IS_MASTER_RANK (simul->rank_id==0)

//simul
struct simul_t
{
  int rank_id,nranks;    //rank and number of ranks
  int nthreads;          //number of threads
  
  int verbosity_lv;      //debug variables
  bool is_little_endian; //endianness flag
  
  void *returned_malloc_ptr;   //returned global pointer
  vectors_t *vectors;          //vectors manager

#ifdef THREAD_DEBUG
  char glb_barr_file[1024];
  int glb_barr_line;
#endif

  //basic MPI initialization
  void init_MPI_thread(int narg,char **arg);
  void get_MPI_nranks(){MPI_Comm_size(MPI_COMM_WORLD,&nranks);} //get nranks
  void get_MPI_rank(){MPI_Comm_rank(MPI_COMM_WORLD,&rank_id);} //get rank
  
  //only master rank and thread print
  int master_fprintf(FILE *stream,const char *format,...)
  {
    GET_THREAD_ID();
    int ret=0;
    
    if(rank_id==0 && IS_MASTER_THREAD)
      {
        va_list ap;
        va_start(ap,format);
        ret=vfprintf(stream,format,ap);
        va_end(ap);
      }
    
    return ret;
  }
  void print_backtrace_list();
  void abort(int arg);
  
  //start and close
  void start(int narg,char **arg,void(*main_function)(int narg,char **arg));
  void close();
  
  simul_t(int narg,char **arg,void(*main_function)(int narg,char **arg));
  ~simul_t();
private:
  simul_t();
};

extern simul_t *simul;

#endif
