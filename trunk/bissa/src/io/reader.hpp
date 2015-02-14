#ifndef _READER_H
#define _READER_H

namespace bissa
{
  void read_real_vector(double *out,const char *path,const char *record_name,uint64_t nreals_per_site,ILDG_message *mess=NULL);
}

#endif
