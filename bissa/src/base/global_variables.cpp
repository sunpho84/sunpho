#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

#include "new_types/new_types_definitions.hpp"

#ifdef ONLY_INSTANTIATION
 #define EXTERN extern
#else
 #define EXTERN
#endif

namespace bissa
{
  //nomenclature: 
  //-glb is relative to the global grid
  //-loc to the local one
  EXTERN int glb_size[4],glb_vol,glb_spat_vol,glb_volh;
  EXTERN int loc_size[4],loc_vol,loc_spat_vol,loc_volh;
  EXTERN int bulk_vol,non_bw_surf_vol,non_fw_surf_vol;
  EXTERN int surf_vol,bw_surf_vol,fw_surf_vol;
  EXTERN int vsurf_vol,vsurf_volh;
  EXTERN int vdir_bord_vol,vdir_bord_volh;
  EXTERN double glb_vol2,loc_vol2;
  //-lx is lexicografic
  EXTERN coords *glb_coord_of_loclx;
  EXTERN coords *loc_coord_of_loclx;
  EXTERN int *glblx_of_loclx;
  EXTERN int *glblx_of_bordlx;
  EXTERN int *loclx_of_bordlx;
  EXTERN int *surflx_of_bordlx;
  EXTERN int *Wsklx_of_loclx;
  EXTERN int *loclx_of_Wsklx;
  //EXTERN int *Wsklx_hopping_matrix_output_pointer;
  //EXTERN int *Wsklx_hopping_matrix_final_output;
  EXTERN int *loclx_of_bulklx;
  EXTERN int *loclx_of_surflx;
  EXTERN int *loclx_of_non_bw_surflx;
  EXTERN int *loclx_of_non_fw_surflx;
  EXTERN int *loclx_of_bw_surflx;
  EXTERN int *loclx_of_fw_surflx;
  EXTERN int lx_geom_inited;
  //-eo is even-odd
  EXTERN int *loclx_parity;
  EXTERN int *loceo_of_loclx;
  EXTERN int *loclx_of_loceo[2];
  EXTERN int *surfeo_of_bordeo[2];
  EXTERN coords *loceo_neighup[2];
  EXTERN coords *loceo_neighdw[2];
  EXTERN int eo_geom_inited;
  EXTERN int use_eo_geom;
  
  //neighbours of local volume + borders
  EXTERN coords *loclx_neighdw,*loclx_neighup;
  EXTERN coords *loclx_neigh[2];
  
//timings
  EXTERN double tot_time;
  EXTERN int issued_cg_warning; //hacking
#ifdef BENCH
 #ifdef ONLY_INSTANTIATION
   EXTERN double tot_comm_time;
 #else
   EXTERN double tot_comm_time=0;
 #endif
#endif

  //bissa_config parameters
  EXTERN int verb_call;
  EXTERN int verbosity_lv;
  EXTERN int warn_if_not_disallocated;
  EXTERN int use_async_communications;
  EXTERN int warn_if_not_communicated;
  EXTERN coords fix_nranks;
  EXTERN int use_128_bit_precision;
  EXTERN int vnode_paral_dir;
  
  //size of the border
  EXTERN int bord_vol,bord_volh;
  //size along various dir
  EXTERN int bord_dir_vol[4],bord_offset[4];
  EXTERN int bord_offset_eo[2][8]; //eo, 8 dirs

#ifdef USE_MPI
  EXTERN int start_lx_bord_send_up[4],start_lx_bord_rece_up[4];
  EXTERN int start_lx_bord_send_dw[4],start_lx_bord_rece_dw[4];
  EXTERN int start_eo_bord_send_up[4],start_eo_bord_rece_up[4];
  EXTERN int start_eo_bord_send_dw[4],start_eo_bord_rece_dw[4];
  
  //volume, plan and line communicator
  EXTERN MPI_Comm cart_comm;
  EXTERN MPI_Comm plan_comm[4];
  EXTERN MPI_Comm line_comm[4];
#endif
  //ranks
  EXTERN int rank,nranks,cart_rank;
  EXTERN coords rank_coord;
  EXTERN coords rank_neigh[2],rank_neighdw,rank_neighup;
  EXTERN coords plan_rank,line_rank,line_coord_rank;
  EXTERN coords nrank_dir;
  EXTERN int grid_inited;
  EXTERN int nparal_dir;
  EXTERN coords paral_dir;
  
/////////////////////////////////////////////// threads //////////////////////////////////////////
#ifdef USE_THREADS

 #ifdef THREAD_DEBUG
   EXTERN int glb_barr_line;
   EXTERN char glb_barr_file[1024];
  #if THREAD_DEBUG >=2 
    EXTERN rnd_gen *delay_rnd_gen;
    EXTERN int *delayed_thread_barrier;
  #endif
 #endif

 #ifndef ONLY_INSTANTIATION
  bool thread_pool_locked=true;
  unsigned int nthreads=1;
 #else
  EXTERN bool thread_pool_locked;
  EXTERN unsigned int nthreads;
 #endif
  
  EXTERN void *broadcast_ptr;
  EXTERN double *glb_double_reduction_buf;

  EXTERN void(*threaded_function_ptr)();
#endif
  
  //endianness
  EXTERN int little_endian;

  //global input file handle
  EXTERN FILE *input_global;
  
  //vectors
  EXTERN int max_required_memory;
  EXTERN int required_memory;
  EXTERN void *main_arr;
  EXTERN bissa_vect main_vect;
  EXTERN bissa_vect *last_vect;
  EXTERN void *return_malloc_ptr;
  
  //random generator stuff
  EXTERN rnd_gen glb_rnd_gen;
  EXTERN int glb_rnd_gen_inited;
  EXTERN rnd_gen *loc_rnd_gen;
  EXTERN int loc_rnd_gen_inited;
  EXTERN enum rnd_t rnd_type_map[6]
#ifndef ONLY_INSTANTIATION
  ={RND_ALL_PLUS_ONE,RND_ALL_MINUS_ONE,RND_Z2,RND_Z2,RND_Z4,RND_GAUSS}
#endif
    ;
  
  /////////////////////////////////////////// buffered comm ///////////////////////////////////
  
  EXTERN int ncomm_allocated;
  EXTERN int comm_in_prog;
  
  //buffers
  EXTERN size_t recv_buf_size,send_buf_size;
  EXTERN char *recv_buf,*send_buf;
  
  //communicators
#ifdef USE_MPI
  EXTERN comm_t lx_double_comm,eo_double_comm;
#endif
  
}

#undef EXTERN
