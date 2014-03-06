#ifndef _FIELDS_HPP
#define _FIELDS_HPP

#include "random.hpp"

//fields
struct fields_t
{
  rnd_gen_t *loc_rnd_gen;                 //random number generators
  
  fields_t(pars_t *pars);
  ~fields_t();
};

#endif
