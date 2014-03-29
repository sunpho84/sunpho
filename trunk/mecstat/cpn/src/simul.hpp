#ifndef _SIMUL_HPP
#define _SIMUL_HPP

#include <mpi.h>
#include <mpi.h>

#include "per_site_neighs.hpp"
#include "random.hpp"
#include "threads.hpp"

class vectors_t;

extern bool simul_started;

#define IS_MASTER_RANK (simul->rank_id==0)

//simul
struct simul_t
{
  int rank_id,nranks;    //rank and number of ranks
  int nthreads;          //number of threads
  char **thread_res_arr; //array to make reductions

  int verbosity_lv;      //debug variables
  bool is_little_endian; //endianness flag
  
  void *returned_malloc_ptr;   //returned global pointer
  vectors_t *vectors;          //vectors manager
  
  char *comm_buff;      //communication buffer
  int comm_buff_size;   //buffer size
  
  per_site_neighs_t *first_neighbors_per_site;   //implement first neighbors connections

  rnd_gen_t glb_rnd_gen;   //random generator
  
#ifdef THREAD_DEBUG
  //used for thread barrier debugging
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
  void init_glb_rnd_gen(int seed);
  void close();
  
  simul_t(int narg,char **arg,void(*main_function)(int narg,char **arg));
  ~simul_t();
private:
  simul_t();
};

extern simul_t *simul;

#include "vectors.hpp"

#endif
