#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef ONLY_INSTANTIATION
 #define EXTERN extern
#else
 #define EXTERN
#endif

#include "system.hpp"

#define ncopies 8
#define nmu 2
#define FW 0
#define BW 1
#define RAN2_NTAB 32

typedef int coords_t[nmu];
typedef int neighs_t[2*nmu];
typedef bool spin_t;

EXTERN system_t *systems[ncopies];                           //system copies
EXTERN int nsweeps,nconfs;                                   //number of sweeps and tot confs to be generated
EXTERN int L,nsites;                                         //lattice size
EXTERN neighs_t *neighs;                                     //neighbors
EXTERN double beta,flip_prob;                                //beta and probablity to flip single spin
