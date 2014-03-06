#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <execinfo.h>
#include <stdlib.h>

#include "debug.hpp"
#include "system.hpp"

bool system_started;
system_t *system;

//basic MPI initialization
void system_t::init_MPI_thread(int narg,char **arg)
{
  if(system_started) CRASH("system already started");
    
  int provided;
  MPI_Init_thread(&narg,&arg,MPI_THREAD_SERIALIZED,&provided);
  system_started=true;
}

//start the system
void system_t::start(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
  //init base things
  init_MPI_thread(narg,arg);

  //this must be done before everything otherwise rank non properly working  
  //get the number of rank and the id of the local one
  get_MPI_nranks();
  get_MPI_rank();
  
  //associate sigsegv with proper handle
  signal(SIGSEGV,signal_handler);
  signal(SIGFPE,signal_handler);
  signal(SIGXCPU,signal_handler);
  
  //print SVN version and configuration and compilation time
  master_printf("Initializing nissa, version: %s\n",SVN_VERSION);
  master_printf("Configured at %s with flags: %s\n",CONFIG_TIME,CONFIG_FLAGS);
  master_printf("Compiled at %s of %s\n",__TIME__,__DATE__);

}

//write the list of called routines
void system_t::print_backtrace_list()
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

void system_t::print_backtrace_list()
{
}

//initialize MPI and threads
system_t::system_t(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
  start(narg,arg,main_function);
  
  close();
}
