#ifndef _BISSA_HPP
#define _BISSA_HPP

//including config.hpp
#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "base/close.hpp"
#include "base/debug.hpp"
#include "base/global_variables.hpp"
#include "base/init.hpp"
#include "base/linalgs.hpp"
#include "base/macros.hpp"
#include "base/random.hpp"
#include "base/thread_macros.hpp"
#include "base/vectors.hpp"

#include "communicate/communicate.hpp"

#include "geometry/geometry_eo.hpp"
#include "geometry/geometry_lx.hpp"

#include "io/checksum.hpp"
#include "io/endianness.hpp"
#include "io/ILDG_File.hpp"
#include "io/input.hpp"
#include "io/reader.hpp"
#include "io/writer.hpp"

#include "new_types/new_types_definitions.hpp"

#include "operations/remap_vector.hpp"

#include "routines/ios.hpp"
#include "routines/math_routines.hpp"
#include "routines/mpi_routines.hpp"
#ifdef USE_THREADS
 #include "routines/thread.hpp"
#endif

#endif
