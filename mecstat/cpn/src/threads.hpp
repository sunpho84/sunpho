#ifndef _THREAD_HPP
#define _THREAD_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>
#include <stdint.h>

//barriers
void thread_barrier_internal();
#ifdef THREAD_DEBUG
  #define THREAD_BARRIER_FORCE() thread_barrier_internal()
  #define THREAD_BARRIER()       thread_barrier_with_check(__FILE__,__LINE__)
  void thread_barrier_with_check(const char *file,int line);
#else
 #define THREAD_BARRIER_FORCE() thread_barrier_internal()
 #define THREAD_BARRIER()       thread_barrier_without_check()
 void thread_barrier_without_check();
#endif

//parallel for
#define CHUNK_WORKLOAD(START,CHUNK_LOAD,END,EXT_START,EXT_END,CHUNK_ID,NCHUNKS) \
  int WORKLOAD=EXT_END-EXT_START,                                       \
    CHUNK_LOAD=(WORKLOAD+NCHUNKS-1)/NCHUNKS,				\
    START=EXT_START+CHUNK_ID*CHUNK_LOAD,				\
    END=START+CHUNK_LOAD< EXT_END ? START+CHUNK_LOAD : EXT_END
#define CHUNK_FOR(INDEX,EXT_START,EXT_END,CHUNK_ID,NCHUNKS)		\
  for(CHUNK_WORKLOAD(START,CHUNK_LOAD,END,EXT_START,EXT_END,CHUNK_ID,NCHUNKS),INDEX=START;INDEX<END;INDEX++)
#define PARALLEL_FOR(INDEX,START,END)				\
  CHUNK_FOR(INDEX,START,END,thread_id,simul->nthreads)

#include "debug.hpp"

#define GET_THREAD_ID() uint32_t thread_id=omp_get_thread_num()
#define IS_MASTER_THREAD (thread_id==0)

#include "simul.hpp"

//flush the cache
inline void cache_flush()
{
 #pragma omp flush
}

//barrier
inline void thread_barrier_internal()
{
  #pragma omp barrier
}

#endif
