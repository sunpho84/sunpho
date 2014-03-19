#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>
#include <stdlib.h>

#include "debug.hpp"
//#include "random.hpp"
#include "threads.hpp"

#if THREAD_DEBUG>=2

//wait previously delayed threads
void wait_for_delayed_threads()
{
  GET_THREAD_ID();
  
  if(delayed_thread_barrier[thread_id]==0)
    {
      if(IS_MASTER_RANK && VERBOSITY_LV3) printf("thread %d waiting for delayed threads\n",THREAD_ID);
      thread_barrier_internal();
    }
} 

//select a new state for the delay
void select_new_delay_pattern()
{
  GET_THREAD_ID();
  
  //extract a random switch:
  // if 0, we exec immediately
  // if 1, we postpone to the end of the barrier
  
  enum delay_pattern{DELAY_RANDOMLY,DELAY_SLAVES};
  const delay_pattern picked=DELAY_SLAVES;
  
  switch(picked)
    {
    case DELAY_RANDOMLY:
      //delayed_thread_barrier[thread_id]=(int)rnd_get_unif(delay_rnd_gen+thread_id,0,2);
      CRASH("unspoorted");
      break;
    case DELAY_SLAVES:
      delayed_thread_barrier[thread_id]=!IS_MASTER_THREAD;
      break;
    default:
      CRASH("Unknown delay pattern %d",picked);
    }
}

//delay marked threads
void delay_marked_threads()
{
  GET_THREAD_ID();
  
  if(delayed_thread_barrier[thread_id]==1)
    {
      if(IS_MASTER_THREAD && VERBOSITY_LV3) printf("thread %d will delay its execution\n",thread_id);
      thread_barrier_internal();
    }
}
#endif

#if THREAD_DEBUG>=1
//check that the local barrier correspond to global one
void check_barrier(const char *barr_file,int barr_line)
{
  GET_THREAD_ID();
  
  if(VERBOSITY_LV3)
    printf("thread %d rank %d barrier call on line %d of file %s)\n",thread_id,simul->rank_id,barr_line,barr_file);
  
  //copy the barrier id to the global ref
  if(IS_MASTER_THREAD)
    {
      strcpy(simul->glb_barr_file,barr_file);
      simul->glb_barr_line=barr_line;
      cache_flush();
    }
  thread_barrier_internal();
  
  //check
  if(!IS_MASTER_THREAD)
    if(simul->glb_barr_line!=barr_line||strcmp(simul->glb_barr_file,barr_file))
      CRASH("Thread %d found barrier on line %d of file %s when master thread invoked it at line %d of file %s)",
	    thread_id,barr_line,barr_file,simul->glb_barr_line,simul->glb_barr_file);
}
#endif

//thread barrier without line number
#if THREAD_DEBUG>=1  
void thread_barrier_with_check(const char *barr_file,int barr_line)
#else
void thread_barrier_without_check()
#endif
{
  //if something was delayed, make it advance
#if THREAD_DEBUG>=2
  wait_for_delayed_threads();
#endif
  
  //true barrier
  thread_barrier_internal();
  
  //check coherence of called barrier
#if THREAD_DEBUG>=1
  check_barrier(barr_file,barr_line);
#endif
  
  //selct new delay pattern and delay
#if THREAD_DEBUG>=2
  select_new_delay_pattern();
  delay_marked_threads();
#endif
}
