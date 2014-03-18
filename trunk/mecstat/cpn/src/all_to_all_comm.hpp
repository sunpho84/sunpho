#ifndef _ALL_TO_ALL_COMM_HPP
#define _ALL_TO_ALL_COMM_HPP

#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include <map>
#include <vector>

#include "simul.hpp"
#include "simul.hpp"

struct all_to_all_gathering_list_t : std::map<int,int> {};
struct all_to_all_scattering_list_t : std::vector<std::pair<int,int> > {};
struct temp_all_to_all_build_t
{
  int *nper_rank_to_temp,*nper_rank_fr_temp;
  int *out_buf_cur_per_rank,*in_buf_cur_per_rank;
  std::map<int,int> rank_to_map_list_ranks_to,rank_fr_map_list_ranks_fr;
  temp_all_to_all_build_t()
  {
    nper_rank_to_temp=ALLOCATE("nper_rank_to_temp",simul->nranks,int);
    nper_rank_fr_temp=ALLOCATE("nper_rank_fr_temp",simul->nranks,int);
  }
  ~temp_all_to_all_build_t()
  {
    FREE(nper_rank_to_temp);
    FREE(nper_rank_fr_temp);
    FREE(out_buf_cur_per_rank);
    FREE(in_buf_cur_per_rank);
  }
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
  void common_setup_part1(temp_all_to_all_build_t &build);
  void common_setup_part2(int nel_note,int *&buf_note,int nranks_note,int *list_ranks_note,int *buf_note_off_per_rank,int *nper_rank_note,int *buf_expl,int nranks_expl,int *list_ranks_expl,int *buf_expl_off_per_rank,int *nper_rank_expl);
  all_to_all_comm_t() {};
};

#endif
