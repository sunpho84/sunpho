#ifndef _FIELDS_HPP
#define _FIELDS_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <mpi.h>

#include "geometry.hpp"

template <class T> class field_t
{
public:
  T *data;
private:
  field_t();
};

#endif
