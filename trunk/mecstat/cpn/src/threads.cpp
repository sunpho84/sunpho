#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <omp.h>

void thread_barrier_without_check()
{
#pragma omp barrier
}
