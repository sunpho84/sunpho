#ifndef _NEW_TYPES_DEFINITIONS_H
#define _NEW_TYPES_DEFINITIONS_H

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#ifdef USE_MPI
 #include <mpi.h>
#endif

#include <math.h>
#include <stdint.h>
#include <sstream>
#include <map>
#include <vector>
#include <string>

#include "base/macros.hpp"

namespace bissa
{  
  ///////////////// New types ///////////////////
  
  typedef int coords[NDIM];
  typedef double momentum_t[NDIM];
  
  typedef int intpair[2];
  typedef uint32_t checksum[2];
  
#ifdef USE_MPI
#ifdef USE_MPI_IO
  typedef MPI_Offset ILDG_Offset;
  typedef MPI_File ILDG_File;
#else
  typedef off_t ILDG_Offset;
  typedef FILE* ILDG_File;
#endif
#endif
  
  //Random types
  enum rnd_t{RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_UNIF,RND_Z2,RND_Z4,RND_GAUSS};
  //Starting condition for a gauge conf
  enum start_conf_cond_t{UNSPEC_START_COND,HOT_START_COND,COLD_START_COND};
  
  ///////////////////// New structures ////////////////////
  
  //The structure for the random generator
  struct rnd_gen
  {
    int idum;
    int idum2;
    int iv[RAN2_NTAB];
    int iy;
  };
  
  //bissa vector
  struct bissa_vect
  {
    int nel;
    int size_per_el;
    
    char tag[BISSA_VECT_STRING_LENGTH];
    char type[BISSA_VECT_STRING_LENGTH];
    
    char file[BISSA_VECT_STRING_LENGTH];
    int line;
    
    bissa_vect *prev;
    bissa_vect *next;
    
    uint32_t flag;
    
    //padding to keep memory alignment
    char pad[BISSA_VECT_ALIGNMENT-(3*sizeof(int)+2*sizeof(void*)+3*BISSA_VECT_STRING_LENGTH+sizeof(uint32_t))%BISSA_VECT_ALIGNMENT];
  };

  //all to all communicators initializing structure
  struct all_to_all_gathering_list_t : std::map<int,int>
  {int add_conf_link_for_paths(coords g,int mu);};
  struct all_to_all_scattering_list_t : std::vector<std::pair<int,int> > {};
  struct temp_build_t
  {
    int *nper_rank_to_temp,*nper_rank_fr_temp;
    int *out_buf_cur_per_rank,*in_buf_cur_per_rank;
    std::map<int,int> rank_to_map_list_ranks_to,rank_fr_map_list_ranks_fr;
    temp_build_t();
    ~temp_build_t();
  };

  //all to all communicators
  struct all_to_all_comm_t
  {
    int nel_out,nel_in;
    int nranks_fr,*list_ranks_fr,*in_buf_dest,*nper_rank_fr,*in_buf_off_per_rank;
    int nranks_to,*list_ranks_to,*out_buf_source,*nper_rank_to,*out_buf_off_per_rank;
    
    all_to_all_comm_t(all_to_all_gathering_list_t &gl);
    all_to_all_comm_t(all_to_all_scattering_list_t &sl);
    ~all_to_all_comm_t();
    void communicate(void *out,void *in,int bps,void *buf_out=NULL,void *buf_in=NULL,int tag=-1);

    void setup_knowing_where_to_send(all_to_all_scattering_list_t &sl);
    void setup_knowing_what_to_ask(all_to_all_gathering_list_t &gl);
    void setup_nper_rank_other_temp(int *nper_rank_other_temp,int *nper_rank_temp);
    void common_setup_part1(temp_build_t &build);
    void common_setup_part2(int nel_note,int *&buf_note,int nranks_note,int *list_ranks_note,int *buf_note_off_per_rank,int *nper_rank_note,int *buf_expl,int nranks_expl,int *list_ranks_expl,int *buf_expl_off_per_rank,int *nper_rank_expl);
    all_to_all_comm_t() {};
  };
}
  
#include "base/vectors.hpp"

namespace bissa
{  
  //ILDG header
  struct ILDG_header
  {
    uint32_t magic_no;
    uint16_t version;
    uint16_t mbme_flag;
    uint64_t data_length;
    char type[128];
  };
  
  //store messages
  struct ILDG_message
  {
    bool is_last;
    char *data;
    char *name;
    uint64_t data_length;
    ILDG_message *next;
  };
  
  //ILDG file view
  struct ILDG_File_view
  {
#ifdef USE_MPI
#ifdef USE_MPI_IO
    MPI_Datatype etype;
    MPI_Datatype ftype;
    MPI_Offset view_pos;
    MPI_Offset pos;
#endif
#endif
    char format[100];
  };
  
  //storable vector
  ILDG_message* ILDG_string_message_append_to_last(ILDG_message *mess,const char *name,const char *data);
  template<class T> struct storable_vector_t : std::vector<T>
  {
    //append to last message
    ILDG_message *append_to_message_with_name(ILDG_message &mess,const char *name)
    {
      std::ostringstream os;
      os.precision(16);
      for(typename std::vector<T>::iterator it=this->begin();it!=this->end();it++) os<<*it<<" ";
      return ILDG_string_message_append_to_last(&mess,name,os.str().c_str());
    }
    //convert from a text message
    void convert_from_message(ILDG_message &mess)
    {
      std::istringstream is(mess.data);
      T temp;
      while(is>>temp) this->push_back(temp);
    }
  };
  
#ifdef USE_MPI
  //out and in buffer
  struct comm_t
  {
    //destinations and source ranks
    int send_rank[2*NDIM],recv_rank[2*NDIM];
    //requests and message
    MPI_Request requests[4*NDIM];
    int nrequest,imessage;
    
    //communication in progress
    int comm_in_prog;
    //local size
    size_t nbytes_per_site;
    //size of the message
    uint64_t tot_mess_size;
    //offsets
    int send_offset[2*NDIM],message_length[2*NDIM],recv_offset[2*NDIM];
    
    //constructor
    bool initialized;
    comm_t(){initialized=false;}
  };
#endif
}
#endif
  
