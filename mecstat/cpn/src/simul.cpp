#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <execinfo.h>
#include <stdlib.h>

#include "debug.hpp"
#include "simul.hpp"
#include "svnversion.hpp"
#include "threads.hpp"

bool simul_started;
simul_t *simul;

//basic MPI initialization
void simul_t::init_MPI_thread(int narg,char **arg)
{
  if(simul_started) CRASH("simul already started");
    
  int provided;
  MPI_Init_thread(&narg,&arg,MPI_THREAD_SERIALIZED,&provided);
  simul_started=true;
}

//start the simul
void simul_t::start(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
  //init base things
  init_MPI_thread(narg,arg);

  //this must be done before everything otherwise rank non properly working  
  //get the number of rank and the id of the local one
  get_MPI_nranks();
  get_MPI_rank();
  printf("%d\n",rank);
  //associate sigsegv with proper handle
  signal(SIGSEGV,signal_handler);
  signal(SIGFPE,signal_handler);
  signal(SIGXCPU,signal_handler);
  
  //print SVN version and configuration and compilation time
  MASTER_PRINTF("Initializing nissa, version: %s\n",SVN_VERSION);
  MASTER_PRINTF("Configured at %s with flags: %s\n",CONFIG_TIME,CONFIG_FLAGS);
  MASTER_PRINTF("Compiled at %s of %s\n",__TIME__,__DATE__);
  

  //verb
}

//write the list of called routines
void simul_t::print_backtrace_list()
{
  void *callstack[128];
  int frames=backtrace(callstack,128);
  char **strs=backtrace_symbols(callstack,frames);
    
  //only master rank, but not master thread
  if(rank==0)
    {
      printf("Backtracing...\n");
      for(int i=0;i<frames;i++) printf("%s\n",strs[i]);
    }
    
  free(strs);
}

//close
void simul_t::close()
{
  MPI_Barrier(MPI_COMM_WORLD);
  MASTER_PRINTF("   Ciao!\n\n");
  MPI_Finalize();
}

//abort
void simul_t::abort(int err)
{
  GET_THREAD_ID();
  printf("thread %d on rank %d aborting\n",thread_id,rank);
  MPI_Abort(MPI_COMM_WORLD,0);
}

//initialize MPI and threads
simul_t::simul_t(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
  start(narg,arg,main_function);
  
  close();
}
