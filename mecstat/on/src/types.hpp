#ifndef _TYPES_HPP
#define _TYPES_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "macros.hpp"

using namespace std;

typedef int coords[NDIMS];
enum start_cond_t{COLD,HOT,LOAD};

#endif
