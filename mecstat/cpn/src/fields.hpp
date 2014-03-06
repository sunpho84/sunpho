#ifndef _FIELDS_HPP
#define _FIELDS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>

#include "geometry.hpp"

//structure for a field
struct field_t
{
  geometry_t *geometry;
  field_t(geometry_t *geometry): geometry(geometry) {};
private:
  field_t();
};

#endif
