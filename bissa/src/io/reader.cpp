#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <string.h>

#include "base/debug.hpp"
#include "base/global_variables.hpp"
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
  //////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //read a real vector
  void read_real_vector(double *out,const char *path,const char *record_name,uint64_t nreals_per_site,ILDG_message *mess=NULL)
  {
    master_printf("Reading vector: %s\n",path);
    
    //Take inital time
    double start_time=take_time();
    
    //open file
    ILDG_File file=ILDG_File_open_for_read(path);
    
    //search the record
    ILDG_header header;
    int found=ILDG_File_search_record(header,file,record_name,mess);
    if(!found) crash("Error, record %s not found.\n",record_name);
    
    //check the size of the data block
    int loc_nreals_tot=nreals_per_site*loc_vol;
    uint64_t nbytes=header.data_length;
    uint64_t nbytes_per_site_read=nbytes/glb_vol;
    if(nbytes_per_site_read>nreals_per_site*sizeof(double))
      crash("Opsss! The file contain %d bytes per site and it is supposed to contain not more than %d!",
	    nbytes_per_site_read,nreals_per_site*sizeof(double));
    
    //read
    ILDG_File_read_ildg_data_all(out,file,header);
    
    //check read size
    uint64_t nbytes_per_site_float=nreals_per_site*sizeof(float);
    uint64_t nbytes_per_site_double=nreals_per_site*sizeof(double);
    
    //read the checksum
    checksum read_check={0,0};
    ILDG_File_read_checksum(read_check,file);
    
    //close the file
    ILDG_File_close(file);
    
    //check precision
    int single_double_flag=-1;
    const char single_double_str[2][10]={"single","double"};
    if(nbytes_per_site_read==nbytes_per_site_float) single_double_flag=0;
    if(nbytes_per_site_read==nbytes_per_site_double) single_double_flag=1;
    if(single_double_flag==-1)
      crash("Opsss! The file contain %d bytes per site and it is supposed to contain: %d (single) or %d (double)",
	    nbytes_per_site_read,nbytes_per_site_float,nbytes_per_site_double);    
    verbosity_lv3_master_printf("Vector is stored in %s precision\n",single_double_str[single_double_flag]);
    
    //change endianess
    if(little_endian)
      {
	if(single_double_flag==0) change_endianness((float*)out,(float*)out,loc_nreals_tot);
	else change_endianness(out,out,loc_nreals_tot);
      }
    
    //check the checksum
    if(read_check[0]!=0||read_check[1]!=0)
      {
	master_printf("Checksums read:      %#010x %#010x\n",read_check[0],read_check[1]);
	
	//compute checksum
	checksum comp_check;
	checksum_compute_bissa_data(comp_check,out,nbytes_per_site_read,(single_double_flag+1)*32);
	
	//print the comparison between checksums
	master_printf("Checksums computed:  %#010x %#010x\n",comp_check[0],comp_check[1]);
	if((read_check[0]!=comp_check[0])||(read_check[1]!=comp_check[1]))
	  master_printf("Warning, checksums do not agree!\n");
      }
    else master_printf("Data checksum not found.\n");
    
    //cast to double if needed
    if(single_double_flag==0) floats_to_doubles_same_endianness(out,(float*)out,loc_nreals_tot);

    set_borders_invalid(out);    
    verbosity_lv2_master_printf("Total time elapsed including possible conversion: %f s\n",take_time()-start_time);
  }
}
