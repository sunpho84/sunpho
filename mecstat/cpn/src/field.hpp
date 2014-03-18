#ifndef _FIELDS_HPP
#define _FIELDS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>

#include "geometry.hpp"
#include "local_mask.hpp"

//structure for a field
template <class T> struct field_t
{
  T *data;
  local_mask_t *mask;
  field_t(local_mask_t *mask): mask(mask) {};
private:
  field_t();
};

#endif
