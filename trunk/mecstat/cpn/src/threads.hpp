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

//cast a pointer from master thread to all the others
#define CAST_PTR_FROM_MASTER_THREAD(A)					\
  (IS_MASTER_THREAD?                                                    \
  (THREAD_BARRIER(),simul->returned_malloc_ptr=A,THREAD_BARRIER(),(typeof(A))simul->returned_malloc_ptr): \
   (THREAD_BARRIER(),THREAD_BARRIER(),(typeof(A))simul->returned_malloc_ptr))

#endif
