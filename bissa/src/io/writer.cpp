#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/global_variables.hpp"
#include "base/debug.hpp"
#include "base/linalgs.hpp"
#include "base/vectors.hpp"
#include "geometry/geometry_lx.hpp"
#include "new_types/new_types_definitions.hpp"
#include "routines/ios.hpp"

#include "checksum.hpp"
#include "endianness.hpp"
#include "ILDG_File.hpp"

namespace bissa
{
  //Write a vector of double, in 32 or 64 bits according to the argument
  void write_double_vector(ILDG_File &file,double *data,size_t nreals_per_site,size_t nbits,const char *header_message,ILDG_message *mess=NULL)
  {
    if(nbits!=32 && nbits!=64) crash("Error, asking %u precision, use instead 32 or 64\n",nbits);
    
    //take initial time
    double time=-take_time();
    
    //write all the messages
    if(mess!=NULL) ILDG_File_write_all_messages(file,mess);
    
    //compute float or double site
    size_t nreals_loc=nreals_per_site*loc_vol;
    size_t nbytes_per_real=nbits/8;
    size_t nbytes_per_site=nreals_per_site*nbytes_per_real;
    
    //buffer to reorder data in ILDG format and change endianness
    char *buffer=bissa_malloc("buffer",nreals_loc*nbytes_per_real,char);
    
    //possibly reduce to 32 bit
    if(nbits==64) parallel_memcpy(buffer,data,nreals_loc*nbytes_per_real);
    else doubles_to_floats_same_endianness((float*)buffer,data,nreals_loc);
    
    //compute the checksum
    checksum check;
    checksum_compute_bissa_data(check,buffer,nbytes_per_site,nbytes_per_real*8);
    
    //change endianness if needed
    if(little_endian)
      {
	if(nbits==64) change_endianness((double*)buffer,(double*)buffer,nreals_loc);
	else change_endianness((float*)buffer,(float*)buffer,nreals_loc);
      }
    
    //write
    ILDG_File_write_ildg_data_all(file,buffer,nbytes_per_site,header_message);
    
    //append the checksum
    ILDG_File_write_checksum(file,check);
    
    //delete the swapped data
    bissa_free(buffer);
    
    //take final time
    time+=take_time();
    verbosity_lv2_master_printf("Time elapsed in writing: %f s\n",time);
  }
}
