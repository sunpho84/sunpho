#ifndef _THREAD_HPP
#define _THREAD_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>
#include <stdint.h>

#define GET_THREAD_ID() uint32_t thread_id=omp_get_thread_num()
#define IS_MASTER_THREAD (thread_id==0)

//flush the cache
inline void cache_flush()
{
 #pragma omp flush
}

//barriers
void thread_barrier_internal();
#ifdef THREAD_DEBUG
 #define THREAD_BARRIER_FORCE() thread_barrier_internal()
 #define THREAD_BARRIER()       thread_barrier_with_check(__FILE__,__LINE__)
 void thread_barrier_with_check(__FILE__,__LINE__);
#else
 #define THREAD_BARRIER_FORCE() thread_barrier_internal()
 #define THREAD_BARRIER()       thread_barrier_without_check()
 void thread_barrier_without_check();
#endif

#endif
