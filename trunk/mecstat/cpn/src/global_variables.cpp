#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef ONLY_INSTANTIATION
 #define EXTERN extern
#else
 #define EXTERN
#endif

#include "system.hpp"

#define nmu 2
#define FW 0
#define BW 1
#define RAN2_NTAB 32

typedef int coords_t[nmu];
typedef int neighs_t[2*nmu];
typedef bool spin_t;

EXTERN system_t *systems;                           //system copies
EXTERN int L,nsites;                                //lattice size
EXTERN neighs_t *neighs;                            //neighbors
