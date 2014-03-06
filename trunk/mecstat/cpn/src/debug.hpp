#ifndef _DEBUG_HPP
#define _DEBUG_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>
#include <stdio.h>

#define CRASH(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)

void internal_crash(int line,const char *file,const char *templ,...);

#endif
