#ifndef _MACROS_HPP
#define _MACROS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#define OMELYAN_LAMBDA 0.1931833

#define RAN2_NTAB 32

#define NDIMS 2

#define CRASH(...) internal_crash(__LINE__,__FILE__,__VA_ARGS__)

#define TINY 1.e-13
#define HALFWAY_PRECISION 1.e-8

#define SKIP_TEST true
#define DO_TEST false

#define DRAW_EACH 100

#define PREC 12

#endif
