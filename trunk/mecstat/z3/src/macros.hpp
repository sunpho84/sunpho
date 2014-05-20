#ifndef _MACROS_HPP
#define _MACROS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define RAN2_NTAB 32

#define NDIMS 3

#define CRASH(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)

#define HOT 1
#define COLD 0

#endif
