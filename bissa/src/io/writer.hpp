#ifndef _writer_hpp
#define _writer_hpp

#include <mpi.h>

namespace bissa
{
  void write_double_vector(ILDG_File &file,double *data,size_t nreals_per_site,size_t nbits,const char *header_message,ILDG_message *mess=NULL);
}

#endif
