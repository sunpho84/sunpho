#ifndef _INIT_HPP
#define _INIT_HPP

#ifndef EXTERN_INIT
 #define EXTERN_INIT extern
#endif

enum start_cond_t{COLD,HOT,LOAD};

struct read_pars_t
{
  int seed;
  start_cond_t start_cond;
  int nterm,nsweep,nmicro,use_hmc;
};

void read_input(read_pars_t &read_pars,const char *path);
void init(int &base_isweep,read_pars_t &read_pars);

EXTERN_INIT int init_time;

#endif
