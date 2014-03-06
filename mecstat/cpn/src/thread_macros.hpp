#ifndef _THREAD_MACROS_HPP
#define _THREAD_MACROS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>
#include <stdint.h>

#define GET_THREAD_ID() uint32_t thread_id=omp_get_thread_num()
#define IS_MASTER_THREAD (thread_id==0)

#endif
