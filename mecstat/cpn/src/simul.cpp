#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

#include "debug.hpp"
#include "endianness.hpp"
#include "simul.hpp"
#include "svnversion.hpp"
#include "threads.hpp"
#include "utils.hpp"

bool simul_started;
simul_t *simul;

//basic MPI initialization
void simul_t::init_MPI_thread(int narg,char **arg)
{
  int provided;
  MPI_Init_thread(&narg,&arg,MPI_THREAD_SERIALIZED,&provided);
  simul_started=true;
}

//start the simul
void simul_t::start(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
  //set verbosity_lv
  verbosity_lv=2;
  
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
  MASTER_PRINTF("Initializing mecstat, version: %s\n",SVN_VERSION);
  MASTER_PRINTF("Configured at %s with flags: %s\n",CONFIG_TIME,CONFIG_FLAGS);
  MASTER_PRINTF("Compiled at %s of %s\n",__TIME__,__DATE__);
  
  //initialize the first vector
  vectors=new vectors_t();
  
  //set no communication buffer
  comm_buff=NULL;
  comm_buff_size=0;
  
  //check endianness
  is_little_endian=get_little_endianness();
  if(is_little_endian) MASTER_PRINTF("System endianness: little (ordinary machine)\n");
  else MASTER_PRINTF("System endianness: big (BG, etc)\n");
  
  //get number of threads
  #pragma omp parallel
  {
    nthreads=omp_get_num_threads();
  }
  MASTER_PRINTF("Using %u threads\n",nthreads);
    
  //init the thread array
  thread_res_arr=NEW_ARRAY_NON_BLOCKING("thread_res_arr",char*,nthreads);
  
  //some init missing
  
  //now start the threads
  #pragma omp parallel
  {
    //start internal main
    main_function(narg,arg);
  }
}

//only master rank and thread print
int simul_t::master_fprintf(FILE *stream,const char *format,...)
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

//write the list of called routines
void simul_t::print_backtrace_list()
{
  void *callstack[128];
  int frames=backtrace(callstack,128);
  char **strs=backtrace_symbols(callstack,frames);
    
  //only master rank, but not master thread
  if(IS_MASTER_RANK)
    {
      printf("Backtracing...\n");
      for(int i=0;i<frames;i++) printf("%s\n",strs[i]);
    }
    
  free(strs);
}

//close
void simul_t::close()
{
  //delete the thread buffer
  DELETE_NON_BLOCKING(thread_res_arr);
  
  //delete communicator buff
  DELETE_NON_BLOCKING(comm_buff);  

  MPI_Barrier(MPI_COMM_WORLD);
  //print information over the maximum amount of memory used
  MASTER_PRINTF("Maximal memory used during the run: %d bytes (",vectors->max_required_memory);
  if(IS_MASTER_RANK) fprintf_friendly_filesize(stdout,vectors->max_required_memory);
  MASTER_PRINTF(") per rank\n\n");
    
  //check wether there are still allocated vectors
  if(vectors->main_vector->next!=NULL && IS_MASTER_RANK)
    {
      printf("Warning, there are still allocated vectors:\n");
      vectors->print_all_contents();
      printf("For a total of %d bytes\n",vectors->total_memory_usage());
    }
  
  //final message
  MPI_Barrier(MPI_COMM_WORLD);
  MASTER_PRINTF("   Ciao!\n\n");
  MPI_Finalize();
}

//init the global random generator
void simul_t::init_glb_rnd_gen(int seed)
{
  GET_THREAD_ID();
  if(IS_MASTER_THREAD) glb_rnd_gen.init(seed);
  THREAD_BARRIER();
}

//abort
void simul_t::abort(int err)
{
  GET_THREAD_ID();
  printf("thread %d on rank %d aborting\n",thread_id,rank_id);
  MPI_Abort(MPI_COMM_WORLD,0);
}

//initialize MPI and threads
simul_t::simul_t(int narg,char **arg,void(*main_function)(int narg,char **arg))
{
  //check immediately that it is not already started
  if(simul_started) CRASH_SOFTLY("simulation already started");
  
  simul=this; //because must assigned before starting
  start(narg,arg,main_function);
  
  close();
}

//destructor
simul_t::~simul_t()
{
}
